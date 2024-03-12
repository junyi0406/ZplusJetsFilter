
import awkward as ak
import numba
import numpy as np


def add_SC_eta(events: ak.Array):
    # photons = events.Photon
    tg_theta_over_2 = np.exp(-events.Photon_eta)
    tg_theta = 2* tg_theta_over_2 / (1-tg_theta_over_2*tg_theta_over_2)
    
    # barrel 
    radius = 130
    a0 = np.arctan(events.PV_y/events.PV_x)
    a1 = np.arctan(events.PV_y/events.PV_x) + np.pi
    angle_x0_y0 = ak.where((events.PV_x > 0) , a0, a1)
    alpha = angle_x0_y0 + (np.pi - events.Photon_phi)
    x2y2 = events.PV_x*events.PV_x + events.PV_y*events.PV_y
    sin_beta = np.sqrt(x2y2) / radius * np.sin(alpha)
    beta = abs(np.arcsin(sin_beta))
    gamma = np.pi/2 - alpha - beta 
    l = np.sqrt(radius*radius + x2y2 - 2*radius*np.sqrt(x2y2)*np.cos(gamma))
    z0_zSC = l / tg_theta
    tg_sctheta_EB = radius / (events.PV_z + z0_zSC)
    # endcap
    intersection_z = ak.where(events.Photon_eta > 0, 310, -310)
    base = intersection_z - events.PV_z
    r = base * tg_theta
    crystalx = events.PV_x + r*np.cos(events.Photon_phi)
    crystaly = events.PV_y + r*np.sin(events.Photon_phi)
    tg_sctheta_EE = np.sqrt(crystalx*crystalx + crystaly*crystaly) / intersection_z

    # combine EB and EE
    sctheta = ak.where(events.Photon_isScEtaEB == 1, np.arctan(tg_sctheta_EB), np.arctan(tg_sctheta_EE))
    sctheta = ak.where(sctheta < 0, sctheta + np.pi, sctheta)
    tg_sctheta_over_2 = np.tan(sctheta/2)
    SCEta = - np.log(tg_sctheta_over_2)
    # remove gap
    SCEta = ak.where((events.Photon_isScEtaEB == 1) | (events.Photon_isScEtaEE == 1), SCEta, events.Photon_eta)        
    events["Photon_SCEta"] = SCEta
    return events
    
    