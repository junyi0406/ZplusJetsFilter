

import awkward as ak


def SelectLepton(leptons):
    
    # assign charge
    anti_matter = leptons.pdgId < 0
    leptons["charge"] = ak.where(anti_matter, 1, -1)
    # sort by pt
    sort_lep_pt = ak.argsort(leptons.pt, ascending=False, axis=1)
    leptons = leptons[sort_lep_pt]
    leptons_pair = ak.combinations(leptons, 2, axis=1, replacement = False)
    
    leptons_selection = (
        # pdgId
        (abs(leptons_pair["0"].pdgId) == abs(leptons_pair["1"].pdgId)) 
        & ((leptons_pair["0"].charge + leptons_pair["1"].charge) == 0)
        & ak.where(abs(leptons_pair["0"].pdgId) == 11, (leptons_pair["0"].pt > 25), (leptons_pair["0"].pt > 20), axis = 1)
        & ak.where(abs(leptons_pair["1"].pdgId) == 11, (leptons_pair["1"].pt > 15), (leptons_pair["1"].pt > 10), axis = 1)
        & ak.where(abs(leptons_pair["0"].pdgId) == 11, (abs(leptons_pair["0"].eta) < 2.5), (abs(leptons_pair["0"].eta) < 2.4), axis = 1)
        & ak.where(abs(leptons_pair["1"].pdgId) == 11, (abs(leptons_pair["1"].eta) < 2.5), (abs(leptons_pair["1"].eta) < 2.4), axis = 1)
        & (leptons_pair["0"].pt > leptons_pair["1"].pt)
        & ((leptons_pair["0"] + leptons_pair["1"]).mass > 50)
    )
    leptons_pair = leptons_pair[leptons_selection]
    leptons_pair = ak.singletons(ak.firsts(leptons_pair[ak.argsort(abs((leptons_pair["0"] + leptons_pair["1"]).mass - 91.18), axis=1, ascending=True)], axis=1))
    
    # lep_leading, lep_trailing = [leptons_pair[idx] for idx in "01"]
    lep_leading, lep_trailing = ak.unzip(leptons_pair)    
    
    return (lep_leading, lep_trailing)
    
    
    
    
    

    
    