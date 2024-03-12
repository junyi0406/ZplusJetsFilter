import awkward as ak

class genZplusJetFilter():
    
    def __init__(self) -> None:
        pass
    def set_EMcut(self, hoe, calo, trk):
        self.hoecut = hoe
        self.calocut = calo
        self.trkcut = trk
    
    def filter(self, eventIDX, seeds):
        seeds = genZplusJetFilter.pre_selection(seeds)
        # duplicate_events = genZplusJetFilter.double_counting_removal(seeds)
        EMhealth_events = genZplusJetFilter.EMfilter(self, phos = seeds)
        return eventIDX[EMhealth_events]
    
    def pre_selection(seeds):
        # cut_pt = seeds.pt > 5
        cut_eta = abs(seeds.eta) < 3
        # cut_gap_region = (abs(seeds.eta) < 1.4442) & (abs(seeds.eta) > 1.566)
        return seeds[cut_eta ]
    
    def double_counting_removal(self, eventIDX, phos):
        # double counting removal
        con_pt15 = phos.pt > 15
        con_PromptOrHardProc = ((phos.statusFlags & 1) == 1) | ((phos.statusFlags & 256) == 256)
        con_isZgammaPho = phos.IsoFrixione & con_pt15 & con_PromptOrHardProc
        return eventIDX[ak.all(con_isZgammaPho == False, axis=1)]
        
    def EMfilter(self, phos):
        # EM filter        
        con_genHoverE = (phos.IsoHadrEt / phos.IsoETInCn) < self.hoecut
        con_genCalIso = phos.IsoCalIso < self.calocut
        con_genTrkIso = phos.IsoTkIsoE < self.trkcut
        passEM = con_genHoverE & con_genCalIso & con_genTrkIso
        return ak.sum(passEM, axis=1) > 0
    

    
    