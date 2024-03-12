import numpy as np
import awkward as ak
import coffea.processor
import coffea.hist as hist
from collections import defaultdict
import coffea.nanoevents.methods.candidate

class QAE(coffea.processor.ProcessorABC):
    
    def __init__(self, chunksize=75000, process_name=""):
        self.chunksize = chunksize
        self.calocuts = [i*2 + 2 for i in range(8)]
        self.trkcuts = [i*1.5 + 0.5 for i in range(8)]
        self.hoecuts = [i*0.1+0.2 for i in range(8)]
        
    def pre_selection(pt, eta, isEB, isEE):
        # photon pt-eta pre-selection
        con_reco_pt = pt > 10 
        con_reco_eta = abs(eta) < 2.5
        con_not_gap = isEB | isEE
        return con_reco_pt & con_reco_eta & con_not_gap

    def double_counting_removal(gen_photons):
        # double counting removal
        con_pt15 = gen_photons.pt > 15
        con_PromptOrHardProc = ((gen_photons.statusFlags & 1) == 1) | ((gen_photons.statusFlags & 256) == 256)
        con_isZgammaPho = gen_photons.IsoFrixione & con_pt15 & con_PromptOrHardProc
        return ak.any(con_isZgammaPho, axis=1)
    
    def EMFilter(gen_photons, hoecut, calocut, trkcut):
        # EM filter        
        con_genHoverE = (gen_photons.IsoHadrEt / gen_photons.IsoETInCn) < hoecut
        con_genCalIso = gen_photons.IsoCalIso < calocut
        con_genTrkIso = gen_photons.IsoTkIsoE < trkcut
        passEM = con_genHoverE & con_genCalIso & con_genTrkIso
        return ak.all(passEM == False, axis=1) 
    
    def concatenate_all(**particles_collections):
        new_particles = ak.firsts(list(particles_collections.values())[0])
        for header, particles in list(particles_collections.items())[1:]:
            field_names = particles.fields
            for field in field_names:
                new_particles[header+"_"+field] = ak.fill_none(ak.firsts(particles[field]), np.float32(-999.))  
        return new_particles        
        
        
                
    def process(self, events):
        import warnings
        warnings.filterwarnings('ignore')
        dataset = events.metadata["dataset"]
        dataset, channel = dataset.split("/")
        filename = events.metadata["filename"]
        cutflow = defaultdict(int)
        genpho_fields = [name for name in events.fields if name[:7] == "newpho_" ]
        recopho_fields = [name for name in events.fields if name[:8] == "recopho_"]
        lep_fields = [name[5:] for name in events.fields if name[:5] == "lep1_"]


        reco_photon = ak.zip( {bname[8:]: events[bname] for bname in recopho_fields}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        gen_photon = ak.zip( {bname[7:]: events[bname] for bname in genpho_fields}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        lep_lead = ak.zip( {bname: events["lep1_"+bname] for bname in lep_fields}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        lep_trail = ak.zip( {bname: events["lep2_"+bname] for bname in lep_fields}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)



        mask_dcr = QAE.double_counting_removal(gen_photon)
        EMcutDict = defaultdict(int)
        for calocut in self.calocuts:
            for trkcut in self.trkcuts:
                for hoecut in self.hoecuts:
                    EMcutDict[f"calocut_{calocut:.1f}_trkcut_{trkcut:.1f}_hoecut_{hoecut:.1f}"] = QAE.EMFilter(gen_photon, hoecut, calocut, trkcut)
        
        gen_photon = gen_photon[gen_photon.DgenPartIdx != -999]
        reco_photon = reco_photon[gen_photon.DgenPartIdx]
        # pre-selection
        mask_reco_selection = QAE.pre_selection(reco_photon.pt, reco_photon.eta, reco_photon.isScEtaEB, reco_photon.isScEtaEE)
        gen_photon = gen_photon[mask_reco_selection]
        reco_photon = reco_photon[mask_reco_selection]
        
        sort_reco_pt = ak.argsort(reco_photon.pt, ascending=False, axis=1)
        sort_reco_pt = ak.singletons(ak.firsts(sort_reco_pt))
        gen_photon = gen_photon[sort_reco_pt]
        reco_photon = reco_photon[sort_reco_pt]

        particles_collection = QAE.concatenate_all(gen = gen_photon, reco = reco_photon, lep1 = lep_lead, lep2 = lep_trail)

        
        
        # print(filename)
        import os
        filename = filename.split("/")[-1]
        file_directory = f"/home/JYChen/photon_nanoAOD/miniTree/{dataset}"
        os.system(f"mkdir -p {file_directory}")
        os.system(f"mkdir -p {file_directory}/{channel}")
        
        import uproot 
        with uproot.recreate(f"{file_directory}/{channel}/{filename}") as fout:
            fout["photons_origin"] = ak.flatten(ak.singletons(particles_collection))
            fout["photons_dupli"] = ak.flatten(ak.singletons(particles_collection[mask_dcr]))
            cutflow["photons_origin"] = ak.sum(particles_collection.Category == 11122, axis=0) 
            cutflow["photons_dupli"] = ak.sum(particles_collection.Category[mask_dcr] == 11122, axis=0) 
            for names, EMmasks in EMcutDict.items():
                fout[names] = ak.flatten(ak.singletons(particles_collection[mask_dcr & EMmasks]))
                cutflow[names] = ak.sum(particles_collection.Category[mask_dcr & EMmasks] == 11122, axis=0) 
                
        return cutflow
    

        
    def postprocess(self, accumulater):
        pass
        