
import ANTools
import awkward as ak
import coffea.processor
from collections import defaultdict
import coffea.nanoevents.methods.candidate

class datasetProducer(coffea.processor.ProcessorABC):
    
    def __init__(self, chunksize=75000, process_name=""):
        self.chunksize = chunksize
        self.cutflow = defaultdict(int)
        self.process_name = process_name
    
    def process(self, events, ):
        import warnings
        warnings.filterwarnings('ignore')
        dataset_name = events.metadata["dataset"]
        filename = events.metadata["filename"]
        # eval(f"process_{self.process_name}(events, dataset, filename)")
        # return self.cutflow

    # def process_photon_dataset(self, events, dataset_name, filename):
        import os
        import pandas as pd
        seeds_fields = [name for name in events.fields if name[:6] == "seeds_" ]
        genpho_fields = [name for name in events.fields if name[:7] == "genpho_" ]
        recopho_fields = [name for name in events.fields if name[:7] == "Photon_"]
        seeds = ak.zip( {bname[6:]: events[bname] for bname in seeds_fields}, 
                        with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        genpho = ak.zip( {bname[7:]: events[bname] for bname in genpho_fields}, 
                        with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        recopho = ak.zip( {bname[7:]: events[bname] for bname in recopho_fields}, 
                        with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        # match seed and genpho into same index
        # match_seed = ak.argcartesian([seeds, genpho], axis=1)
        # idxs, idxgenpho = ak.unzip(match_seed[
        #     (seeds[match_seed["0"]].pt == genpho[match_seed["1"]].pt) & 
        #     (seeds[match_seed["0"]].eta == genpho[match_seed["1"]].eta) & 
        #     (seeds[match_seed["0"]].phi == genpho[match_seed["1"]].phi)])
        # seeds = seeds[idxs]
        # genpho = genpho[idxgenpho]
        
        mask_notMatched = genpho.DgenPartIdx != -999
        seeds = seeds[mask_notMatched]
        genpho = genpho[mask_notMatched]
        # sort the recopho together
        idx_sort = ak.argsort(genpho.DgenPartIdx, axis=1)
        seeds = seeds[idx_sort]
        genpho = genpho[idx_sort]
        recopho = recopho[genpho.DgenPartIdx]
        
        # flatten to pandas
        seeds = ak.to_pandas(ak.flatten(seeds))
        genpho = ak.to_pandas(ak.flatten(genpho))
        recopho = ak.to_pandas(ak.flatten(recopho))
        # rename the kinematics
        seeds.rename(columns={name: "seeds_" + name for name in seeds.columns}, inplace=True)
        genpho.rename(columns={name: "gen_" + name for name in genpho.columns}, inplace=True)
        recopho.rename(columns={name: "Photon_" + name for name in recopho.columns}, inplace=True)
        dataset = pd.concat([seeds, genpho, recopho], axis=1)
        # dataset["seeds_IsoHoverE"] = dataset["seeds_IsoHadrEt"] / dataset["seeds_IsoETInCn"]
        # store the dataset 
        filename = filename.split("/")[-1]
        file_directory = f"/home/JYChen/photon_nanoAOD/miniTree/{dataset_name}"
        os.system(f"mkdir -p {file_directory}")
    #     self.save_as_root(dataset, file_dir=file_directory, filename = f"photons_{filename}")

    # def save_as_root(data, file_dir, filename):
        import uproot 
        with uproot.recreate(f"{file_directory}/{filename}") as fout:
            fout["photons"] = dataset
        return self.cutflow
    def postprocess(self, accumulater):
        pass