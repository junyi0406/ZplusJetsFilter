
import coffea.processor
import coffea.nanoevents.methods.candidate
import awkward as ak
from collections import defaultdict
import ANTools
import uproot
import os

class MatchInspector(coffea.processor.ProcessorABC):
    
    def __init__(self, chunksize=75000,):
        self.chunksize = chunksize
  
    
    def process(self, events):
        import warnings
        warnings.filterwarnings('ignore')
        dataset = events.metadata["dataset"]
        filename = events.metadata["filename"]
        cutflow = defaultdict(int)

        genPart = ak.zip({bname[8:]: events[bname] for bname in events.fields if bname[:8] == "GenPart_" },
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior) 
        recopho = ak.zip({bname[7:]: events[bname] for bname in events.fields if bname[:7] == "Photon_"}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        
        cutflow["NEvent"] +=  ak.num(genPart, axis=0)
        # reco_gen = ANTools.matcher(ids = [22], dr = 0.1, dpt = 10000)
        # genPart["DgenPartIdx"] = reco_gen.process(genPart, recopho)
        combination = ak.argcartesian([genPart, recopho], axis=1)
        conID = abs(genPart[combination["0"]].pdgId) == 22
        combination = combination[conID]
        cutflow["NEvent_before_pre_selection"] +=  ak.num(combination[ak.num(combination, axis=1) > 0], axis=0)
        pre_selection = (genPart[combination["0"]].pt > 10) & (abs(genPart[combination["0"]].eta) < 3) & (abs(recopho[combination["1"]].eta) < 2.5)
        combination = combination[pre_selection]
        cutflow["NEvent_after_pre_selection"] +=  ak.num(combination[ak.num(combination, axis=1) > 0], axis=0)
        
        for i in range(1, 2):
            dr = i * 0.1
            condr = genPart[combination["0"]].delta_r(recopho[combination["1"]]) < dr
            cutflow[f"NEvent_with_dR_0.{i}"] +=  ak.num(combination[ak.num(combination[condr], axis=1) > 0], axis=0)
            
        
        return cutflow
    
    def postprocess(self, accumulator):
        pass