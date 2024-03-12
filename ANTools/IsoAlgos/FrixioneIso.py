
import awkward as ak
import numpy as np





def GetFrixioneIso(uncurgp: ak.Array):
    
    r0, xn, eps = 0.05, 1, 1

    con_is22 = uncurgp.pdgId == 22
    con_stable = uncurgp.status == 1
    
    genpho = uncurgp[con_is22 & con_stable]
    genstable = uncurgp[con_stable]
    combination = ak.argcartesian([genpho, genstable])
    
    # Idxgenpho, Idxgenpart = ak.unzip(combination)
    con_not_same_pt = genpho[combination["0"]].pt != genstable[combination["1"]].pt
    con_not_same_eta = genpho[combination["0"]].eta != genstable[combination["1"]].eta
    con_dr04 = genpho[combination["0"]].delta_r(genstable[combination["1"]]) < r0
    
    # it is the angular distance between 
    ang_dist = genpho[combination["0"]].delta_r(genstable[combination["1"]])
    threshold = genpho[combination["0"]].pt * (1-np.cos(genpho[combination["0"]].delta_r(genstable[combination["1"]])))/(1-np.cos(r0))
    
    ang_dist = ang_dist[con_not_same_pt & con_not_same_eta & con_dr04]
    threshold = threshold[con_not_same_pt & con_not_same_eta & con_dr04]
    combination = combination[con_not_same_pt & con_not_same_eta & con_dr04]
    


    idx_level = ak.copy(combination["0"]) 
    idx_level = 0
    idx_mask = (idx_level == combination["0"])
    database = [[0]*len(genpho)] * ak.max(ak.num(genpho, axis=1), axis=0)
    while(ak.any(idx_mask)): 
        sel_gp = genpho[combination[idx_mask]["0"]]
        sel_sgp = genstable[combination[idx_mask]["1"]]
        # isolation and comapre 
        rank_ad = ak.argsort(ang_dist[idx_mask], axis=1)
        ad_masked = ang_dist[idx_mask][rank_ad]
        th_masked = threshold[idx_mask][rank_ad]
        pt_masked = sel_sgp.pt[rank_ad]
        max_iso_level = ak.max(ak.num(ad_masked, axis=1), axis=0)
        frixiiso = [[True]*len(sel_sgp)] * max_iso_level
        for iso_level in range(max_iso_level):
            iso_mask = (ak.local_index(ad_masked)>iso_level) == False
            sub_iso = ak.sum(pt_masked[iso_mask], axis=1)
            sub_th  = ak.max(th_masked[iso_mask], axis=1)
            frixiiso[iso_level] = ak.fill_none(sub_iso < sub_th, True)
        iso = ak.all(frixiiso, axis=0)
        database[idx_level] = iso
        # increae the level
        idx_level += 1
        idx_mask = (idx_level == combination["0"])
    builder = ak.ArrayBuilder()
    for iev in range(len(genpho)):
        builder.begin_list()
        ngenpho = len(genpho[iev])
        for ipho in range(ngenpho):
            iso = database[ipho][iev]
            builder.boolean(iso)
        builder.end_list()

    IsoFrixione = builder.snapshot()
    genpho["IsoFrixione"] = IsoFrixione
    return genpho


