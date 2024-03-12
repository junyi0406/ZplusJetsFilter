import numpy as np
import awkward as ak
import coffea.processor
import coffea.hist as hist
from collections import defaultdict
import coffea.nanoevents.methods.candidate

class IsoCutOptimizer(coffea.processor.ProcessorABC):
    
    def __init__(self, chunksize=75000, channel = -99):
        self.chunksize = chunksize
        self.channel = channel
        self.calocuts = [i*2 + 2 for i in range(8)]
        self.trkcuts = [i*1.5 + 0.5 for i in range(8)]
        self.hoecuts = [i*0.1+0.2 for i in range(8)]

    def process(self, events):
        import warnings
        warnings.filterwarnings('ignore')
        dataset = events.metadata["dataset"]
        filename = events.metadata["filename"]
        cutflow = defaultdict(int)
        genpho_fields = [name for name in events.fields if name[:7] == "newpho_" ]
        # recopho_fields = [name for name in events.fields if name[:7] == "Photon_"]


        gen_photon = ak.zip( {bname[7:]: events[bname] for bname in genpho_fields}, 
            with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        if self.channel == 11:
            gen_photon = gen_photon[gen_photon.channel_ele]
        elif self.channel == 13:
            gen_photon = gen_photon[gen_photon.channel_mu]
        elif self.channel == 15:
            gen_photon = gen_photon[gen_photon.channel_tautau]
        else:
            print("error")
        mask_has_pho = ak.num(gen_photon, axis=1) != 0
        mask_isNotZgamma = ak.all(gen_photon.isZgamma == 0, axis=1)
        
        gen_photon = gen_photon[mask_has_pho & mask_isNotZgamma]
        sort_gen_pt = ak.argsort(gen_photon.pt, ascending=False, axis=1)
        gen_photon = gen_photon[sort_gen_pt]
        for calocut in self.calocuts:
            for trkcut in self.trkcuts:
                for hoecut in self.hoecuts:
                    con_genHoverE = (gen_photon.IsoHadrEt / gen_photon.IsoETInCn) < hoecut
                    con_genCalIso = gen_photon.IsoCalIso < calocut
                    con_genTrkIso = gen_photon.IsoTkIsoE < trkcut
                    regist_gen_photon_is_faked = ak.firsts(gen_photon.Category == 11122)
                    mask_EMKeep = ak.any(con_genHoverE & con_genCalIso & con_genTrkIso, axis=1)
                    mask_twoEMKeep = ak.sum(con_genHoverE & con_genCalIso & con_genTrkIso, axis=1) > 1
                    
                    cutflow[f"calocut_{calocut:.1f}_trkcut_{trkcut:.1f}_hoecut_{hoecut:.1f}"] += ak.sum(regist_gen_photon_is_faked & mask_EMKeep, axis=0)
                    cutflow[f"calocut_{calocut:.1f}_trkcut_{trkcut:.1f}_hoecut_{hoecut:.1f}_two"] += ak.sum(regist_gen_photon_is_faked & mask_twoEMKeep, axis=0)
                    cutflow[f"calocut_{calocut:.1f}_trkcut_{trkcut:.1f}_hoecut_{hoecut:.1f}_EMtotal"] += ak.sum(mask_EMKeep, axis=0)
                    cutflow[f"calocut_{calocut:.1f}_trkcut_{trkcut:.1f}_hoecut_{hoecut:.1f}_EMtotal_two"] += ak.sum(mask_twoEMKeep, axis=0)

                    
        cutflow[f"Events_isfaked"] += ak.sum(ak.firsts(gen_photon.Category == 11122), axis=0)
        cutflow[f"Events_total"] += ak.num(gen_photon, axis=0)
        return cutflow
    
    def postprocess(self, accumulater):
        pass