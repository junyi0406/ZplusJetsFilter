import os 
import numpy as np
import awkward as ak
import coffea.nanoevents.methods.candidate

class LeptonSelector():
    
    def __init__(self) -> None:
        pass
    
    def select_lepton(self, channel, eventIdx, events):
        events = LeptonSelector.add_SC_eta("Electron", events)
        electron = ak.zip({bname[9:]: events[bname] for bname in events.fields if bname[:8] == "Electron"}, 
                with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        muon = ak.zip({bname[5:]: events[bname] for bname in events.fields if bname[:4] == "Muon"}, 
                with_name="PtEtaPhiMCandidate", behavior=coffea.nanoevents.methods.candidate.behavior)
        if channel == 11:
            # for good electron
            eleMVAID = LeptonSelector.SelectLooseMVAIDEle(electron)
            goodEle = (electron.pt > 15) & (abs(electron.SCEta) < 2.5) & (abs(electron.dxy) < 0.5) & (abs(electron.dz) < 1.) & eleMVAID
            return LeptonSelector.selection_z_to_2l(channel, eventIdx, events, goodEle, electron)
        elif channel == 13:
            # for good muon
            muon_mass = 105.7 * 0.001
            LooseMu = LeptonSelector.SelectLooseMuon(muon)
            goodMu = (muon.pt > 10.) & LooseMu
            return LeptonSelector.selection_z_to_2l(channel, eventIdx, events, goodMu, muon)
        else:
            print("unsupported channel!")
            return False
        # additional lepton selection for lepton tag
        # isaddEl = (electron.pt > 7.) & (abs(electron.SCEta) < 2.5) & (abs(electron.dxy) < 0.5) & (abs(electron.dz) < 1.) & eleMVAID
        # isaddMu = LooseMu
                
        

    def SelectLooseMVAIDEle(electron : ak.Array):
        cut_mva_eta1 = ak.broadcast_arrays(electron.pt, [-0.145237])[1]
        cut_mva_eta2 = ak.broadcast_arrays(electron.pt, [-0.0315746])[1]
        cut_mva_eta3 = ak.broadcast_arrays(electron.pt, [-0.032173])[1]
        cut_low_mva_eta1 = ak.broadcast_arrays(electron.pt, [0.604775])[1]
        cut_low_mva_eta2 = ak.broadcast_arrays(electron.pt, [0.628743])[1]
        cut_low_mva_eta3 = ak.broadcast_arrays(electron.pt, [0.896462])[1]

        mask_high_pt = electron.pt >= 10.
        mask_eta1 = electron.SCEta < 0.8
        mask_eta2 = (electron.SCEta >= 0.8) & (electron.SCEta < 1.479)
        mask_eta3 = electron.SCEta >= 1.479
        isEleGoodMVA = ak.where(mask_high_pt & mask_eta1, electron.mvaIso_Fall17V2 > cut_mva_eta1, electron.mvaIso_Fall17V2 > cut_low_mva_eta1)
        isEleGoodMVA = ak.where(mask_high_pt & mask_eta2, electron.mvaIso_Fall17V2 > cut_mva_eta2, electron.mvaIso_Fall17V2 > cut_low_mva_eta2) | isEleGoodMVA
        isEleGoodMVA = ak.where(mask_high_pt & mask_eta3, electron.mvaIso_Fall17V2 > cut_mva_eta3, electron.mvaIso_Fall17V2 > cut_low_mva_eta3) | isEleGoodMVA        
        return isEleGoodMVA
        
    def SelectLooseMuon(muons):
        
        isLooseMuon = (abs(muons.eta) < 2.4) & (abs(muons.dxy) < 0.5) & (abs(muons.dz) < 1.) & (muons.sip3d < 4.) & (muons.isGlobal | (muons.isTracker & (muons.nStations > 0))) #& (muons. != 2) & muons.muon.pt > 5.
        isTrkHighPtMuon = (muons.isTracker) & (muons.nStations > 1) & (muons.nTrackerLayers > 5) & (abs(muons.dxy) < 0.2) & (abs(muons.dz) < 0.5) & ((muons.ptErr / muons.pt) < 0.3)# & muPixelHits > 0
        isHZZTightMuon = (isLooseMuon & (muons.pt < 200.) & muons.isPFcand) | (isLooseMuon & (muons.pt >= 200.) & (muons.isPFcand | isTrkHighPtMuon))
        # muIso03 = LeptonSelector.MuonNoneZeroIso03(muons.pfRelIso03_chg, muons.PFNeuIso03, muons.PFPhoIso03, muons.PFPUIso03)
        # isHZZIsoMuon = (muIso03/muons.pt) < 0.35
        isGoodMu = isHZZTightMuon# & isHZZIsoMuon
        return isGoodMu
        


    def add_SC_eta(particle: str, events: ak.Array):
        # photons = events.Photon
        tg_theta_over_2 = eval(f"np.exp(-events.{particle}_eta)")
        tg_theta = 2* tg_theta_over_2 / (1-tg_theta_over_2*tg_theta_over_2)
        
        # barrel 
        radius = 130
        a0 = np.arctan(events.PV_y/events.PV_x)
        a1 = np.arctan(events.PV_y/events.PV_x) + np.pi
        angle_x0_y0 = ak.where((events.PV_x > 0) , a0, a1)
        alpha = eval(f"angle_x0_y0 + (np.pi - events.{particle}_phi)")
        x2y2 = events.PV_x*events.PV_x + events.PV_y*events.PV_y
        sin_beta = np.sqrt(x2y2) / radius * np.sin(alpha)
        beta = abs(np.arcsin(sin_beta))
        gamma = np.pi/2 - alpha - beta 
        l = np.sqrt(radius*radius + x2y2 - 2*radius*np.sqrt(x2y2)*np.cos(gamma))
        z0_zSC = l / tg_theta
        tg_sctheta_EB = radius / (events.PV_z + z0_zSC)
        # endcap
        intersection_z = eval(f"ak.where(events.{particle}_eta > 0, 310, -310)")
        base = intersection_z - events.PV_z
        r = base * tg_theta
        crystalx = eval(f"events.PV_x + r*np.cos(events.{particle}_phi)")
        crystaly = eval(f"events.PV_y + r*np.sin(events.{particle}_phi)")
        tg_sctheta_EE = np.sqrt(crystalx*crystalx + crystaly*crystaly) / intersection_z

        if particle == "Electron":
            events[f"{particle}_isScEtaEB"] = eval(f"abs(events.{particle}_eta < 1.4442)")
            events[f"{particle}_isScEtaEE"] = eval(f"abs(events.{particle}_eta > 1.566) & abs(events.{particle}_eta < 2.5)")
        # combine EB and EE
        sctheta = eval(f"ak.where(events.{particle}_isScEtaEB == 1, np.arctan(tg_sctheta_EB), np.arctan(tg_sctheta_EE))")
        sctheta = ak.where(sctheta < 0, sctheta + np.pi, sctheta)
        tg_sctheta_over_2 = np.tan(sctheta/2)
        SCEta = - np.log(tg_sctheta_over_2)
        # remove gap
        SCEta = eval(f"ak.where((events.{particle}_isScEtaEB == 1) | (events.{particle}_isScEtaEE == 1), SCEta, events.{particle}_eta)")        
        events[f"{particle}_SCEta"] = SCEta
        return events
    
    # def RochesterCorrection(RoccoR *rc, int HasMC, int nMu, Long64_t ev, VecF_t muPhi, VecF_t muEta, VecF_t muPt, VecI_t muCharge, VecI_t muTrkLayers):
    #     int year = 16;
    #     ROOT::RVec<float> vec; vec.clear();
    #     Double_t SF[nMu], randd;   
    #     TRandom3 *r = new TRandom3(ev);
    #     for(int i=0; i<nMu; i++){
    #         if(HasMC) SF[i] = rc->kSmearMC(muCharge[i], muPt[i], muEta[i], muPhi[i], muTrkLayers[i], r->Rndm(), 0, 0);
    #         else SF[i] = rc->kScaleDT(muCharge[i], muPt[i], muEta[i], muPhi[i], 0, 0);
    #         vec.push_back(muPt[i]*SF[i]);
    #     }
    #     return vec;

    def MuonNoneZeroIso03(muPFChIso03, muPFNeuIso03, muPFPhoIso03, muPFPUIso03):
        iso = muPFNeuIso03 + muPFPhoIso03 - 0.5*muPFPUIso03
        return muPFChIso03 + ak.where(iso < 0, 0, iso)
    
    # Find first lepton pair with invariant mass closest to Z mass
    def reco_z_to_2l(isGoodlep, leptons):
        leptons = leptons[isGoodlep]
        LepCombinations = ak.combinations(leptons, 2, axis=1)
        # leptons-pair selection
        charge_selection = ((LepCombinations["0"].charge + LepCombinations["1"].charge) == 0)
        eta_selection = (abs(LepCombinations["0"].eta) > 1.2) & (abs(LepCombinations["1"].eta) > 1.2)
        eta_product_selection = abs(LepCombinations["0"].eta)*abs(LepCombinations["1"].eta) > 0
        delta_phi = LepCombinations["0"].phi - LepCombinations["1"].phi
        delta_phi_selection = ak.where(delta_phi < -np.pi, delta_phi + np.pi, delta_phi - np.pi) < 70.*np.pi/180
        Z_mass_selection = (LepCombinations["0"] + LepCombinations["1"]).mass > 50.
        LepCombinations = LepCombinations[charge_selection & eta_selection & eta_product_selection & delta_phi_selection & Z_mass_selection]
        # Find first lepton pair with invariant mass closest to Z mass
        lep_pair_mass = (LepCombinations["0"] + LepCombinations["1"]).mass 
        cloest_to_z_mass = ak.min(abs(lep_pair_mass - 91.18), axis=1) == abs(lep_pair_mass - 91.18)
        lep_leading, lep_trailing = ak.unzip(LepCombinations[cloest_to_z_mass])    

        return (lep_leading, lep_trailing)
    
    def selection_z_to_2l(channel, eventsIDX, events, isGoodLep, leptons):
        # HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v >> 15
        # HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v >> 42
        hlt = False
        if channel == 11:
            hlt = events.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL # HLTEleMuX >> 40
            hlt = hlt | events.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ # HLTEleMuX >> 5
        elif channel == 13:
            hlt = events.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL # HLTEleMuX >> 41
            hlt = hlt | events.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ # HLTEleMuX >> 14
        else:
            print("unsupported channel")
            os._exit(0)
        atLeast_two_lep = ak.sum(isGoodLep, axis=1) > 1
        ent_selection = atLeast_two_lep & hlt
        prunedIDX = eventsIDX[ent_selection]
        leadingLep, trailingLep = LeptonSelector.reco_z_to_2l(isGoodLep[ent_selection], leptons[ent_selection])
        
        return (prunedIDX, leadingLep, trailingLep)
       