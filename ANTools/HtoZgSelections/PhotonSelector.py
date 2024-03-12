import os 
import numpy as np
import awkward as ak
import coffea.nanoevents.methods.candidate

class PhotonSelector():
    
    def __init__(self) -> None:
        pass
    
    def select_pho(self, eventIdx: ak.Array, events: ak.Array, selectedLep: tuple) -> ak.Array:
        # return the good photon (particle level selection)
        photons = ak.zip({bname[7:]: events[bname] for bname in events.fields if bname[:6] == "Photon"}, 
                with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        leadsLep, trailLep = selectedLep
        photons = PhotonSelector.pho_dr_cuts(leadsLep, trailLep, photons)
        phoCalibEt = photons.pt 
        isGoodPho_EB = (phoCalibEt > 15.) & photons.electronVeto & photons.isScEtaEB
        isGoodPho_EE = (phoCalibEt > 15.) & photons.electronVeto & photons.isScEtaEE
        isGoodPho = isGoodPho_EB | isGoodPho_EE
        photons = photons[isGoodPho]
        photons = PhotonSelector.select_higgs(selectedLep, photons)
        # pruned no good photon events
        atLeast_one_photon = ak.fill_none(ak.num(photons.pt, axis=1), 0) > 0

        return (eventIdx[atLeast_one_photon], leadsLep[atLeast_one_photon], trailLep[atLeast_one_photon], photons[atLeast_one_photon])
    
    def pho_dr_cuts(lep1, lep2, phos):
        llg = ak.argcartesian([lep1, lep2, phos], axis=1)
        dr_cut = (phos[llg["2"]].delta_r(lep1[llg["0"]]) > 0.4) & (phos[llg["2"]].delta_r(lep2[llg["1"]]) > 0.4)
        return phos[ak.unzip(llg[dr_cut])[2]]    
    
    def select_higgs(selectedLep: tuple, selectedPho: ak.Array):
        
        leadLep, trailLep = selectedLep
        phos = PhotonSelector.reco_higgs(leadLep, trailLep, selectedPho)
        return phos

    def reco_higgs(lep1, lep2, phos):
        
        phos["mass"] = 0.
        phos["charge"] = 0.
        Z = lep1+lep2
        higgs = ak.cartesian([Z, phos], axis=1)
        higgs_mass = (higgs["0"]+higgs["1"]).mass
        higgs_slection = (
            (higgs_mass < 180.) & (higgs_mass > 100.) &
            ((higgs_mass +higgs["0"].mass) > 185.) &
            (higgs["1"].pt/higgs_mass > 15./110.)
        )
        phos = ak.unzip(higgs[higgs_slection])[1]
        # phos = phos[prunedPhoIDX]
        return phos

