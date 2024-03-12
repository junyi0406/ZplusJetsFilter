import awkward as ak

def tag_category(parts: ak.Array, sources: ak.Array):
    # pdgID criteria
    con_mpid_lep   = (abs(parts.pdgMId) == 11) | (abs(parts.pdgMId) == 13)
    con_mpid_fsrtau = abs(parts.pdgMId) == 15
    con_mpid_meson = (parts.pdgMId == 111) | (parts.pdgMId == 221) | (parts.pdgMId == 223) 
    con_mpid_charged_meson = abs(parts.pdgMId == 521) | abs(parts.pdgMId == 523) | abs(parts.pdgMId == 411) | abs(parts.pdgMId == 413) | abs(parts.pdgMId == 433)
    con_mpid_quark = abs(parts.pdgMId) < 6
    con_gmpid_tau = abs(parts.pdgGMId) == 15
    con_gmpid_not_tau = abs(parts.pdgGMId) != 15
    # con_gmpid_quark = (abs(parts.pdgGMId) < 6) | (parts.pdgGMId == 21)
    # con_nogmpid = sources[parts.IdxMother].IdxMother == -1
    # statusFlags criteria
    con_isPrompt = (parts.statusFlags & 1) == 1
    con_fromHardProcess = (parts.statusFlags & 256) == 256
    # charged mesons
    # con_mpid_charge_mesons = (abs(parts.pdgGMId) == 423) | (abs(parts.pdgGMId) == 513) | (abs(parts.pdgGMId) == 523) 
    # photon from tau
    parts["Category"] = -99
    parts["Category"] = ak.where(con_isPrompt & con_fromHardProcess & con_mpid_quark, 22, parts.Category)
    parts["Category"] = ak.where(con_isPrompt & con_mpid_lep, 1122, parts.Category)
    parts["Category"] = ak.where(con_gmpid_not_tau & con_mpid_meson , 11122, parts.Category)
    parts["Category"] = ak.where(con_gmpid_tau & con_mpid_meson , 1522, parts.Category)
    parts["Category"] = ak.where(con_isPrompt & con_mpid_fsrtau, 151522, parts.Category)
    parts["Category"] = ak.where(con_mpid_charged_meson, 43322, parts.Category)
    return parts