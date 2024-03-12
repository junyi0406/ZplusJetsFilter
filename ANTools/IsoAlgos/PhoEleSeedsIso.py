
import awkward as ak


def pickseeds(particles):
    filter_eta_min_= 0.
    filter_eta_max_= 2.6
    seedthreshold  = 5.
    
    con_pt = particles.pt > seedthreshold
    con_eta = (abs(particles.eta) < filter_eta_max_) & (abs(particles.eta) > filter_eta_min_)
    con_pdgid = (particles.pdgId == 22) #| (abs(particles.pdgId) == 11)
    con_status = particles.status == 1
    return (con_pt) & (con_eta) & (con_pdgid) & (con_status)
  
def deltaRxyAtEE(parts_seed, parts_cur):
    import numpy as np
    EECapZ = 304.5
    part0_theta = 2 * np.arctan(np.exp(-parts_seed.eta))
    part0_rxy = EECapZ * np.tan(part0_theta)
    part0_x   = np.cos(parts_seed.phi) * part0_rxy
    part0_y   = np.sin(parts_seed.phi) * part0_rxy

    part2_theta = 2 * np.arctan(np.exp(-parts_cur.eta))
    part2_rxy = EECapZ * np.tan(part2_theta)
    part2_x   = np.cos(parts_cur.phi) * part2_rxy
    part2_y   = np.sin(parts_cur.phi) * part2_rxy

    dxy = np.sqrt((part2_x-part0_x)*(part2_x-part0_x) + (part2_y-part0_y)*(part2_y-part0_y))
    return dxy
 
def build_isolation(parts, part_idx, builder, ismtk):
    s, un, cv = parts
    # sid, unid, cvid = ak.unzip(part_idx)
    
    idx_level = ak.copy(part_idx["0"]) 
    idx_level = 0
    idx_mask = (idx_level == part_idx["0"])
    database = [[0]*len(s)] * ak.max(ak.num(s, axis=1), axis=0)
    while(ak.any(idx_mask)): # check is there any seeds
        # isolation calculation in each mask layer
        sel_sp = s[part_idx[idx_mask]["0"]]
        sel_cg = cv[part_idx[idx_mask]["2"]]
        pt_reqset = ak.mask(sel_cg, ak.num(sel_cg, axis=1) != 0).pt
        if ismtk:
            iso = ak.num(pt_reqset, axis=1) > 0
        else:
            iso = ak.sum(pt_reqset, axis=1)
            
        # package isolation
        iso = ak.fill_none(iso, 0)
        # print(iso )
        database[idx_level] = iso
        # increae the level
        idx_level += 1
        idx_mask = (idx_level == part_idx["0"])
    # database (event*number of seed)
    # print(np.shape(database))
    for i in range(len(cv)): # event loop
        # shape should be the same as curved gen-particles
        builder.begin_list()
        nSeed = len(s[i])
        # looping the database
        for iS in range(nSeed):
            iso = database[iS][i]
            if ismtk:
                builder.boolean(iso)
            else:
                builder.real(iso)
        builder.end_list()
    return builder
            
    
def GetElePhoSeedsIso(genpart: ak.Array, curvedgenpart: ak.Array):
    
    EBmaxEta=1.479
    consizeEE = 15.
    consizeIso = 0.2                                     
    mask_seed = pickseeds(curvedgenpart)
    seeds = curvedgenpart[mask_seed]
    # the candidates in single event contain seeds, uncurved, curved
    part_combination = ak.argcartesian([seeds, genpart, curvedgenpart], axis=1)
    Idxseed, IdxUnc, IdxCur = ak.unzip(part_combination)
    # list all the material of the selection
    # con_softPt = curvedgenpart[IdxCur].pt > 1.
    con_TkPt = curvedgenpart[IdxCur].pt > 5.
    con_isChg = curvedgenpart[IdxCur].charge != 0 # neutral = False
    con_deta = abs(seeds[Idxseed].eta - curvedgenpart[IdxCur].eta) < 0.03 
    con_dphi = abs(seeds[Idxseed].phi - curvedgenpart[IdxCur].phi) < 0.2
    con_drxyEE   = deltaRxyAtEE(seeds[Idxseed], curvedgenpart[IdxCur]) < consizeEE
    con_is22 = curvedgenpart[IdxCur].pdgId == 22
    
    # seeds_Mparticle = curvedgenpart[seeds[Idxseed].IdxMother]
    # seeds_Mparticle_islepton = (abs(seeds_Mparticle.pdgId) == 11) | (abs(seeds_Mparticle.pdgId) == 13) | (abs(seeds_Mparticle.pdgId) == 15) 
    # con_seedsIsPrompt = (seeds[Idxseed].statusFlags & 1) == 1
    # con_islepton = (abs(curvedgenpart[IdxCur].pdgId) == 11) | (abs(curvedgenpart[IdxCur].pdgId) == 13) | (abs(curvedgenpart[IdxCur].pdgId) == 15)
    # con_isfromhardprompt = ((curvedgenpart[IdxCur].statusFlags & 256) == 256) & ((curvedgenpart[IdxCur].statusFlags & 1) == 1)
    # con_isnothardprolepton = ak.where(con_islepton & con_isfromhardprompt, False, True)
    
    con_is11 = abs(curvedgenpart[IdxCur].pdgId) == 11
    con_is211 = abs(curvedgenpart[IdxCur].pdgId) == 211
    con_is321 = abs(curvedgenpart[IdxCur].pdgId) == 321
    con_drcsIso01  = seeds[Idxseed].delta_r(genpart[IdxUnc]) < consizeIso
    con_drcsIso02  = seeds[Idxseed].delta_r(curvedgenpart[IdxCur]) < consizeIso
    con_samepart = (genpart[IdxUnc].pdgId == curvedgenpart[IdxCur].pdgId) & (genpart[IdxUnc].pt == curvedgenpart[IdxCur].pt) & (genpart[IdxUnc].eta == curvedgenpart[IdxCur].eta)
    con_isbarrel = abs(seeds[Idxseed].eta) < EBmaxEta  # for the seed particles
    con_status = (genpart[IdxUnc].status == 1) & (curvedgenpart[IdxCur].status == 1)
    # this step create the mask which is used to be applied on genParticles.
    mask_MatTrk = ak.where(con_isbarrel, 
                (con_samepart) & (con_deta) & (con_dphi) & (con_is11 | con_is211 | con_is321) & (con_TkPt) & (con_status), # condition in barrel 
                (con_samepart) & (con_drxyEE) & (con_is11 | con_is211 | con_is321) & (con_TkPt) & (con_status))             # condition in endcap
    mask_ETInCn = ak.where(con_isbarrel, 
                (con_samepart) & (con_deta) & (con_dphi) & (con_is11 | con_is22 | con_is211 | con_is321) & (con_status), 
                (con_samepart) & (con_drxyEE) & (con_is11 | con_is22 | con_is211 | con_is321) & (con_status))             
    mask_HadrEt = ak.where(con_isbarrel,
                (con_samepart) & (con_deta) & (con_dphi) & (ak.where(con_is11 | con_is22 | con_is211 | con_is321, False, True)) & (con_status),
                (con_samepart) & (con_drxyEE) & (ak.where(con_is11 | con_is22 | con_is211 | con_is321, False, True)) & (con_status))
    mask_CalIso = ak.where(con_isbarrel,
                (con_samepart) & (ak.where((con_deta) & (con_dphi), False, True)) & (((ak.where(con_isChg, False, True)) & (con_drcsIso02) & (ak.where(con_is22, False, True))) | ((con_isChg) & (con_drcsIso01))) & (con_status),
                (con_samepart) & (ak.where(con_drxyEE, False, True)) & (((ak.where(con_isChg, False, True)) & (con_drcsIso02) & (ak.where(con_is22, False, True))) | ((con_isChg) & (con_drcsIso01))) & (con_status))
    mask_TkIsoE = ak.where(con_isbarrel,
                (con_samepart) & (ak.where((con_deta) & (con_dphi), False, True)) & (con_isChg) & (con_drcsIso01) & (con_status),
                (con_samepart) & (ak.where(con_drxyEE, False, True)) & (con_isChg) & (con_drcsIso01) & (con_status))
    # remove hardprocess lepton for FSR photon
    # mask_CalIso = ak.where(seeds_Mparticle_islepton & con_seedsIsPrompt, mask_CalIso & con_isnothardprolepton, mask_CalIso)
    # mask_TkIsoE = ak.where(seeds_Mparticle_islepton & con_seedsIsPrompt, mask_TkIsoE & con_isnothardprolepton, mask_TkIsoE)
    # need to sum up the ratio of energy by the mask (old method)
    IsoMatTrk = build_isolation([seeds, genpart, curvedgenpart], part_combination[mask_MatTrk], ak.ArrayBuilder(), True).snapshot()
    IsoETInCn = build_isolation([seeds, genpart, curvedgenpart], part_combination[mask_ETInCn], ak.ArrayBuilder(), False).snapshot()
    IsoHadrEt = build_isolation([seeds, genpart, curvedgenpart], part_combination[mask_HadrEt], ak.ArrayBuilder(), False).snapshot()
    IsoCalIso = build_isolation([seeds, genpart, curvedgenpart], part_combination[mask_CalIso], ak.ArrayBuilder(), False).snapshot()
    IsoTkIsoE = build_isolation([seeds, genpart, curvedgenpart], part_combination[mask_TkIsoE], ak.ArrayBuilder(), False).snapshot()
    # need to sum up the ratio of energy by the mask (new method)
    
    
    # assign value
    seeds["MatchTrack"] = IsoMatTrk
    seeds["IsoETInCn"] = IsoETInCn # /seeds.pt
    seeds["IsoHadrEt"] = IsoHadrEt # /seeds.pt
    seeds["IsoCalIso"] = IsoCalIso # /seeds.pt
    seeds["IsoTkIsoE"] = IsoTkIsoE # /seeds.pt
    return seeds

