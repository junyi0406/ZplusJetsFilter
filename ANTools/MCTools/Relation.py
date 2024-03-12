

def findMotherCollection(target, source):
    import awkward as ak
    mIdx = target.genPartIdxMother
    mother_particles_idx = ak.where( mIdx != -1, source[mIdx].genPartIdxMother, -1)
    return mother_particles_idx

def findGMotherCollection(target, source):
    import awkward as ak
    mother_collection_idx = findMotherCollection(target, source)
    return ak.where(mother_collection_idx != -1, findMotherCollection(target[mother_collection_idx], source), -1)
    
def tag_relationId(target, source):
    import awkward as ak
    # mother_gp_idx = findMotherCollection(target, source)
    # gmother_gp_idx = findGMotherCollection(target, source)
    mother_gp_idx = target.genPartIdxMother
    gmother_gp_idx = ak.where( mother_gp_idx != -1, source[mother_gp_idx].genPartIdxMother, -1)
    # gmother_gp_idx = ak.where( mother_gp_idx != -1, source[mother_gp_idx].genPartIdxMother, -1)
    target["pdgMId"] = ak.where((mother_gp_idx != -1), source[mother_gp_idx].pdgId, -999)
    target["pdgGMId"] = ak.where(((mother_gp_idx != -1) & (gmother_gp_idx != -1)), source[gmother_gp_idx].pdgId, -999)
    target["Mstatus"] = ak.where((mother_gp_idx != -1), source[mother_gp_idx].status, -999)
    target["GMstatus"] = ak.where(((mother_gp_idx != -1) & (gmother_gp_idx != -1)), source[gmother_gp_idx].status, -999)
    target["MstatusFlags"] = ak.where(((mother_gp_idx != -1)), source[mother_gp_idx].statusFlags, -999)
    target["GMstatusFlags"] = ak.where(((mother_gp_idx != -1) & (gmother_gp_idx != -1)), source[gmother_gp_idx].statusFlags, -999)
    return target
    