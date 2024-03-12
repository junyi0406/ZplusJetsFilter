
import coffea.processor
import coffea.nanoevents.methods.candidate
import awkward as ak
from collections import defaultdict
import ANTools
import uproot
import os
import numpy as np
class ZgSelectionCheck(coffea.processor.ProcessorABC):
    
    def __init__(self, chunksize=75000, ):
        self.chunksize = chunksize
    def set_searchChannel(self, channel: int):
        self.channel = channel
        
    def set_EMcut(self, hoe, calo, trk):
        self.hoecut = hoe
        self.calocut = calo
        self.trkcut = trk   
        
    def process(self, events, ):
        import warnings
        warnings.filterwarnings('ignore')
        dataset = events.metadata["dataset"]
        filename = events.metadata["filename"]
        cutflow = defaultdict(int)
        events = ANTools.add_SC_eta(events)       
        eventIDX = ak.local_index(events, axis=0)
        
        if len(events) == 0:
            return {'cutflow': {dataset: cutflow},}
        genPart = ak.zip({bname[8:]: events[bname] for bname in events.fields if bname[:7] == "GenPart"}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        curvedgenPart = ak.zip({bname[14:]: events[bname] for bname in events.fields if bname[:13] == "CurvedGenPart"}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        recopho = ak.zip({bname[7:]: events[bname] for bname in events.fields if bname[:6] == "Photon"}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        # do the matching privately
        # reco_gen = ANTools.matcher(ids = [11, 22], dr = 0.1, dpt = 10000)
        # genPart["DgenPartIdx"] = reco_gen.process(genPart, recopho)
        # add the isolation variables
        # photon electron seeds isolation
        seeds = ANTools.GetElePhoSeedsIso(genPart, curvedgenPart)
        genpho = ANTools.GetFrixioneIso(genPart)
        # tag the mother pid and GMPID
        genpho = ANTools.tag_relationId(genpho, genPart)
        # categorizations 
        genpho = ANTools.tag_category(genpho, genPart) 
        # match_seed = ak.argcartesian([seeds, genpho], axis=1)
        # idxs, idxgenpho = ak.unzip(match_seed[
        #     (seeds[match_seed["0"]].pt == genpho[match_seed["1"]].pt) & 
        #     (seeds[match_seed["0"]].eta == genpho[match_seed["1"]].eta) & 
        #     (seeds[match_seed["0"]].phi == genpho[match_seed["1"]].phi)])
        # seeds = seeds[idxs]
        # genpho = genpho[idxgenpho]
        # newpho = {bname: seeds[bname] for bname in seeds.fields if (bname in ["pt", "eta", "phi"]) == False}
        # newpho.update({bname: genpho[bname] for bname in genpho.fields})
        # newpho = ak.zip(newpho, with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        sort_genpho_pt = ak.argsort(genpho.pt, ascending=False, axis=1)
        sort_genpho_pt = ak.singletons(ak.firsts(sort_genpho_pt))
        genpho = genpho[sort_genpho_pt]
        
        # HZg selections
        # selector and Filter initialization
        lepSelector = ANTools.LeptonSelector()
        phoSelector = ANTools.PhotonSelector()
        genFilter = ANTools.genZplusJetFilter()
        # print(f"selections_{filename.replace('.root','')}")
        # select required particles
        channelIDX, lep1, lep2 = lepSelector.select_lepton(self.channel, eventIDX, events)
        mask_channel = ZgSelectionCheck.compare_selectedIDX(eventIDX, channelIDX)
        prunedIDX, lep1, lep2, phos = phoSelector.select_pho(channelIDX, events[mask_channel], (lep1, lep2))

        # print(f"filter_{filename.replace('.root','')}")
        # gen filter selections
        genFilter.set_EMcut(hoe = self.hoecut, calo = self.calocut, trk = self.trkcut)
        EMfilter_prunedIDX = genFilter.filter(eventIDX, seeds)
        dcrfilter_prunedIDX = genFilter.double_counting_removal(eventIDX, genpho)
        
        cutflow["total_event"] += ak.num(eventIDX, axis=0)
        cutflow[f"{self.channel}_event"] += ak.num(channelIDX, axis=0)

        
        # print(f"compare_{filename.replace('.root','')}")
        mask_selections = ZgSelectionCheck.compare_selectedIDX(eventIDX, prunedIDX)
        mask_DCRfilter = ZgSelectionCheck.compare_selectedIDX(eventIDX, dcrfilter_prunedIDX)
        mask_EMfilter = ZgSelectionCheck.compare_selectedIDX(eventIDX, EMfilter_prunedIDX)
        if (ak.sum(mask_selections, axis=0) == 0 | (ak.sum(mask_DCRfilter, axis=0) == 0)  | (ak.sum(mask_EMfilter, axis=0) == 0)):
            return {'cutflow': {dataset: cutflow},}
  
        phos["isPassSelection"] =  ak.broadcast_arrays(phos.pt, mask_selections[mask_selections])[1]
        phos["isPassDCRFilter"] = ak.broadcast_arrays(phos.pt, mask_DCRfilter[mask_selections])[1]
        phos["isPassEMFilter"] = ak.broadcast_arrays(phos.pt, mask_EMfilter[mask_selections])[1]
        
        cutflow["pass_selections"] += ak.num(prunedIDX, axis=0)

        cutflow["Zg_EM1"] += ak.sum(mask_channel & mask_selections & mask_EMfilter, axis=0)
        cutflow["Zg_EM2"] += ak.sum(mask_channel & (mask_selections == False) & mask_EMfilter, axis=0)
        cutflow["Zg_EM3"] += ak.sum(mask_channel & mask_selections & (mask_EMfilter == False), axis=0)
        cutflow["Zg_EM4"] += ak.sum(mask_channel & (mask_selections == False) & (mask_EMfilter == False), axis=0)
        
        cutflow["Zg_dupli1"] += ak.sum(mask_channel & mask_selections & mask_DCRfilter, axis=0)
        cutflow["Zg_dupli2"] += ak.sum(mask_channel & (mask_selections == False) & mask_DCRfilter, axis=0)
        cutflow["Zg_dupli3"] += ak.sum(mask_channel & mask_selections & (mask_DCRfilter == False), axis=0)
        cutflow["Zg_dupli4"] += ak.sum(mask_channel & (mask_selections == False) & (mask_DCRfilter == False), axis=0)
        
        cutflow["Zg_dupli_EM1"] += ak.sum(mask_channel & mask_selections & (mask_DCRfilter & mask_EMfilter), axis=0)
        cutflow["Zg_dupli_EM2"] += ak.sum(mask_channel & (mask_selections == False) & (mask_DCRfilter & mask_EMfilter), axis=0)
        cutflow["Zg_dupli_EM3"] += ak.sum(mask_channel & mask_selections & ((mask_DCRfilter & mask_EMfilter) == False), axis=0)
        cutflow["Zg_dupli_EM4"] += ak.sum(mask_channel & (mask_selections == False) & ((mask_DCRfilter & mask_EMfilter) == False), axis=0)



        # add the reco-level variables
        outdict = {
            "genpho": genpho[mask_selections], 
            "seeds": seeds[mask_selections], 
            "GenPart": genPart[mask_selections], 
            "CurvedGenPart": curvedgenPart[mask_selections], 
            "lep1": lep1, 
            "lep2": lep2, 
            "pho": phos
        }
        
        # store the root tree
        filename = filename.split("/")[-1]
        os.system(f"mkdir -p /home/JYChen/photon_nanoAOD/miniTree/{dataset}/{self.channel}")
        print(f"/home/JYChen/photon_nanoAOD/miniTree/{dataset}/{self.channel}/{filename.replace('.root','')}_{int(events.metadata['entrystop']/self.chunksize)}.root")
        with uproot.recreate(
            f"/home/JYChen/photon_nanoAOD/miniTree/{dataset}/{self.channel}/{filename.replace('.root','')}_{int(events.metadata['entrystop']/self.chunksize)}.root"
            ) as fout:
            output_root = {}
            # particle kinematics
            for bname in outdict.keys():
                if not outdict[bname].fields:
                    output_root[bname] = ak.packed(ak.without_parameters(outdict[bname]))
                else:
                    b_nest = {}
                    for n in outdict[bname].fields:
                        # print(type(ak.packed(ak.without_parameters(outdict[bname][n]))))
                        # print(ak.packed(ak.without_parameters(outdict[bname][n])))
                        # print(ak.type(ak.packed(ak.without_parameters(outdict[bname][n]))) )
                        if bname in ["lep1", "lep2", "pho"]:
                            b_nest[n] = ak.packed(ak.without_parameters(outdict[bname][n])).to_list()
                        else:
                            b_nest[n] = ak.packed(ak.without_parameters(outdict[bname][n]))
                    output_root[bname] = ak.zip(b_nest)
            fout["Events"] = output_root   


        return {'cutflow': {dataset: cutflow},}
    
    def compare_selectedIDX(origin, IDX1):
        # return the same shape like original idx array
        mask_IDX1 =  ak.any((ak.broadcast_arrays(origin, [IDX1])[1]) == origin, axis=1)
        return mask_IDX1
    
    def postprocess(self, accumulater):
        pass
        
