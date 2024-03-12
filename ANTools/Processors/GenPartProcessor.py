
import coffea.processor
import coffea.nanoevents.methods.candidate
import awkward as ak
from collections import defaultdict
import ANTools
import uproot
import os

class GenPartSelector(coffea.processor.ProcessorABC):
    
    def __init__(self, chunksize=75000, ):
        self.chunksize = chunksize
        self.nphocut = {"origin": 1}
        self.calocut = {"origin": 10}
        self.trkcut = {"origin": 5}
        self.hoecut = {"origin": 0.5}
        
        
    def process(self, events):
        import warnings
        warnings.filterwarnings('ignore')
        dataset = events.metadata["dataset"]
        filename = events.metadata["filename"]
        cutflow = defaultdict(int)
        events = ANTools.add_SC_eta(events)       
        
        if len(events) == 0:
            return {'cutflow': {dataset: cutflow},}
        recopho_fields = [name for name in events.fields if name[:7] == "Photon_"]
        # genpho_fields = [name for name in events.fields if name[:14] == "CurvedGenPart_"]
        
        recopho = ak.zip({bname[7:]: events[bname] for bname in recopho_fields}, 
                        with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
            
        genPart = ak.zip({
            "pt": events.GenPart_pt,
            "eta": events.GenPart_eta,
            "phi": events.GenPart_phi,
            "mass": events.GenPart_mass,
            "pdgId": events.GenPart_pdgId,
            "status": events.GenPart_status,
            "statusFlags": events.GenPart_statusFlags,
            # "charge": events.GenPart_charge,
            "IdxMother": events.GenPart_genPartIdxMother
        }, with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        curvedgenPart = ak.zip({
            "pt": events.CurvedGenPart_pt,
            "eta": events.CurvedGenPart_eta,
            "phi": events.CurvedGenPart_phi,
            "mass": events.CurvedGenPart_mass,
            "pdgId": events.CurvedGenPart_pdgId,
            "status": events.CurvedGenPart_status,
            "statusFlags": events.CurvedGenPart_statusFlags,
            "charge": events.CurvedGenPart_charge,
            "IdxMother": events.CurvedGenPart_genPartIdxMother
        }, with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        # do the matching privately
        reco_gen = ANTools.matcher(ids = [11, 22], dr = 0.1, dpt = 10000)
        genPart["DgenPartIdx"] = reco_gen.process(genPart, recopho)

        # pick up good lepton
        genlep = genPart[(abs(genPart.pdgId) == 11) | (abs(genPart.pdgId) == 13)]
        leading_lep, trailing_lep = ANTools.SelectLepton(genlep)

        cutflow["NEvents"] += ak.num(recopho, axis=0)
        # channel determination
        mask_is_Z_to_tautau = ak.any(genPart[genPart[(abs(genPart.pdgId) == 15)].IdxMother].pdgId == 23, axis=1)
        Z_ele = (abs(genPart.pdgId) == 11) & (genPart[genPart.IdxMother].pdgId == 23) 
        Z_mu = (abs(genPart.pdgId) == 13) & (genPart[genPart.IdxMother].pdgId == 23) 
        mask_is_Z_to_ee = ak.where(mask_is_Z_to_tautau, False, ak.any(Z_ele, axis=1))
        mask_is_Z_to_mumu = ak.where(mask_is_Z_to_tautau, False, ak.any(Z_mu, axis=1))
        mask_has_recopho = ak.num(recopho, axis=1) != 0


        cutflow["NEvents_Ztoee"] += ak.sum(mask_is_Z_to_ee, axis=0)
        cutflow["NEvents_Ztomumu"] += ak.sum(mask_is_Z_to_mumu, axis=0)
        cutflow["NEvents_ZtoTauTau"] += ak.sum(mask_is_Z_to_tautau, axis=0)
        
        # add the isolation variables
        # photon electron seeds isolation
        seeds = ANTools.GetElePhoSeedsIso(genPart, curvedgenPart)
        genpho = ANTools.GetFrixioneIso(genPart)
        # tag the mother pid and GMPID
        seeds = ANTools.tag_relationId(seeds, curvedgenPart)
        genpho = ANTools.tag_relationId(genpho, genPart)
 
        # categorizations 
        genpho = ANTools.tag_category(genpho, genPart)
        # combine gen-photon and seeds photon

        mask_has_pho = ak.num(recopho, axis=1) != 0
        for channel, mask in zip([11, 13, 15], [mask_is_Z_to_ee, mask_is_Z_to_mumu, mask_is_Z_to_tautau]):
            cutflow[f"NEvent_hasphoton_{channel}"] += ak.sum((ak.num(genpho, axis=1) != 0) & mask, axis=0)
            cutflow[f"NEvent_hasrecopho_{channel}"] += ak.sum(mask_has_pho & mask, axis=0)
            cutflow[f"NEvent_hasphoton_recopho_{channel}"] += ak.sum((ak.num(genpho, axis=1) != 0) & mask_has_pho & mask, axis=0)
            cutflow[f"NEvent_nophoton_hasrecopho_{channel}"] += ak.sum((ak.num(genpho, axis=1) == 0) & mask_has_pho & mask, axis=0)
            has_matched_ele = ak.any((abs(genPart.pdgId) == 11) & (genPart.DgenPartIdx != -999), axis=1) 
            cutflow[f"NEvent_nophoton_hasrecopho_from_ele_{channel}"] += ak.sum((ak.num(genpho, axis=1) == 0) & mask_has_pho & mask & has_matched_ele, axis=0)
            cutflow[f"NEvent_nophoton_norecopho_{channel}"] += ak.sum((ak.num(genpho, axis=1) == 0) & (mask_has_pho == False) & mask, axis=0)
            
        match_seed = ak.argcartesian([seeds, genpho], axis=1)
        idxs, idxgenpho = ak.unzip(match_seed[
            (seeds[match_seed["0"]].pt == genpho[match_seed["1"]].pt) & 
            (seeds[match_seed["0"]].eta == genpho[match_seed["1"]].eta) & 
            (seeds[match_seed["0"]].phi == genpho[match_seed["1"]].phi)])
        seeds = seeds[idxs]
        genpho = genpho[idxgenpho]
        # mask no match gen-photon
        # combine gen & reco photon
        newpho = {bname: seeds[bname] for bname in seeds.fields if (bname in ["pt", "eta", "phi"]) == False}
        newpho.update({bname: genpho[bname] for bname in genpho.fields})
        newpho = ak.zip(newpho, with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        mask_has_pho = ak.num(newpho, axis=1) != 0
        for channel, mask in zip([11, 13, 15], [mask_is_Z_to_ee, mask_is_Z_to_mumu, mask_is_Z_to_tautau]):
            for has_genpho in [True, False]:
                outdict = {
                    "newpho": newpho[mask & (mask_has_pho == has_genpho)],
                    "recopho": recopho[mask & (mask_has_pho == has_genpho)],
                    "lep1": leading_lep[mask & (mask_has_pho == has_genpho)],
                    "lep2": trailing_lep[mask & (mask_has_pho == has_genpho)]
                }
                # store the root tree
                filename = filename.split("/")[-1]
                postfix = "with_genpho" if has_genpho == True else "without_genpho"
                os.system(f"mkdir -p /home/JYChen/photon_nanoAOD/miniTree/{dataset}")
                os.system(f"mkdir -p /home/JYChen/photon_nanoAOD/miniTree/{dataset}/{postfix}")
                with uproot.recreate(
                    f"/home/JYChen/photon_nanoAOD/miniTree/{dataset}/{postfix}/{filename.replace('.root','')}_{channel}.root"
                    ) as fout:
                    fout["Events"] = outdict
            

               
        # cutflow["NEvent"] += len(events)
        return cutflow
    
    def postprocess(self, accumulater):
        pass
        