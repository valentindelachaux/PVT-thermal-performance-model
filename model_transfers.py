import math
import copy
import pandas as pd
import numpy as np
import heat_transfer as bht
import ht

import model_ht as mht

from CoolProp.CoolProp import PropsSI

# Need the absorber's fin base temperature T_B - function not used
def qp_fin(componentSpecs, var):
    lambd_abs = componentSpecs["lambd_abs"]
    k_abs = componentSpecs["k_abs"]

    L_af = componentSpecs["L_af"]

    #T_PV = var["T_PV"]
    T_B = var["T_Base_mean"]

    j = var["j"]
    b = var["b"]

    m = var["m"]
    
    q = k_abs*lambd_abs*m*((b/j)-T_B)*math.tanh(m*L_af)
    var["qp_fin"] = q

# Eq. 560.38
def q_tube_fluid(componentSpecs,stepConditions,var):
    """Calculates the heat flux from the fluid to the tube and stores it in var["q_tube_fluid"]
    
    $$
    q_{tube,fluid} = \frac{\dot{m}C_p(T_{fluid,out}-T_{fluid,in})}{L_{tube}N_{harp}}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        T_fluid_in (float): fluid inlet temperature
        var (dict): dictionary containing the variables
    
    Returns:
        None"""
    N_harp = componentSpecs["N_harp"]
    L = componentSpecs["L_tube"]
    mdot = stepConditions["mdot"]
    Cp = var["Cp"]    

    T_fluid_in = var["T_fluid_in"]
    
    T_f_out = var["T_fluid_out"]
    res = (mdot*Cp*(T_f_out-T_fluid_in))/(L*N_harp)

    var["q_tube_fluid"] = res

def q_Base_tube(componentSpecs,var):
    """Calculates the heat flux from the base to the tube and stores it in var["q_Base_tube"]
    
    $$
    q_{Base,tube} = -\theta_{Bt}q_{tube,fluid} + \kappa_{Bt}T_{fluid,out} + \epsilon_{Bt}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables

    Returns:
        None"""
    Ka_Bt = var["Ka_Bt"]
    Th_Bt = var["Th_Bt"]
    Ep_Bt = var["Ep_Bt"]

    T_fluid = var["T_fluid_mean"]

    q_tube_fluid = var["q_tube_fluid"]

    var["q_Base_tube"] = -Th_Bt*q_tube_fluid + Ka_Bt*T_fluid + Ep_Bt

def Qdot_sun_glass(componentSpecs,stepConditions,var):

    alpha = componentSpecs["alpha_g"]
    G = stepConditions["G"]

    W = componentSpecs['W']
    L_tube = componentSpecs['L_tube']

    var["Qdot_sun_glass"] = W*L_tube*alpha*G

def Qdot_sun_PV(componentSpecs,stepConditions,var):

    W=componentSpecs["W"]
    L_tube = componentSpecs["L_tube"]

    var["Qdot_sun_PV"] = W*L_tube*var["S"]

# Eq. 560.47
def Qdot_top_conv(componentSpecs,stepConditions,var):

    T_glass_m = var["T_glass"]
    T_amb = stepConditions["T_amb"]

    h_top_g = var["h_top_g"]
    W = componentSpecs["W"]
    L = componentSpecs["L_tube"]

    var["Qdot_top_conv"] = (W*L)*(T_glass_m-T_amb)*h_top_g

def Qdot_top_rad(componentSpecs,stepConditions,var):

    h_r = var["h_rad_g"]
    T_glass = var["T_glass"]
    T_sky = stepConditions["T_sky"]
    W = componentSpecs["W"]

    L = componentSpecs["L_tube"]

    var["Qdot_top_rad"] = W*L*h_r*(T_glass-T_sky)

def Qdot_PV_sky(componentSpecs,stepConditions,var):

    h_r = var["h_rad"]
    T_PV_m = var["T_PV"]
    T_sky = stepConditions["T_sky"]
    W = componentSpecs["W"]

    L = componentSpecs["L_tube"]

    var["Qdot_PV_sky"] = W*L*h_r*(T_PV_m-T_sky)

def Qdot_tube_sky(componentSpecs,stepConditions,var):
    h_r = var["h_rad_tube_sky"]
    p_tube_sky = componentSpecs["p_tube_sky"]
    T_tube = var["T_tube_mean"]
    T_sky = stepConditions["T_sky"]
    W = componentSpecs["W"]

    L = componentSpecs["L_tube"]

    var["Qdot_tube_sky"] = h_r*p_tube_sky*(T_tube-T_sky)

def Qdot_glass_PV(componentSpecs,stepConditions,var):
    
    R_g = componentSpecs["R_g"]
    T_glass_m = var["T_glass"]
    T_PV_m = var["T_PV"]
    W = componentSpecs["W"]

    L = componentSpecs["L_tube"]

    var["Qdot_glass_PV"] = (W*L)*(T_glass_m-T_PV_m)/R_g

def Qdot_PV_plate(componentSpecs,var):

    R_inter = componentSpecs["R_inter"]
    W = componentSpecs["W"]

    T_PV_m = var["T_PV"]
    T_abs_m = var["T_abs_mean"]
    L = componentSpecs["L_tube"]

    var["Qdot_PV_plate"] = (W*L)*(T_PV_m-T_abs_m)/R_inter

def power_balance_1(componentSpecs,var):
    S = var["S"]
    Q1 = var["Qdot_top_conv"]
    Q2 = var["Qdot_top_rad"]
    Q3 = var["Qdot_PV_plate"]
    W = componentSpecs["W"]
    L = componentSpecs["L_tube"]
 
    var["power_balance_1"] = (W*L)*S-Q1-Q2-Q3

def Qdot_PV_Base(componentSpecs,var):

    R_inter = componentSpecs["R_inter"]
    l_B = componentSpecs["l_B"]

    T_PV_mB = var["T_PV_Base_mean"]
    T_Base_m = var["T_Base_mean"]
    L = componentSpecs["L_tube"]

    var["Qdot_PV_Base"] = L*l_B*((T_PV_mB-T_Base_m)/R_inter)

def Qdot_PV_absfin(componentSpecs,var):
    R_inter = componentSpecs["R_inter"]
    L_af = componentSpecs["L_af"]

    T_PV_absfin_mean = var["T_PV_absfin_mean"]
    T_absfin_mean = var["T_absfin_mean"]
    L = componentSpecs["L_tube"]

    var["Qdot_PV_absfin"] = L*2*L_af*((T_PV_absfin_mean-T_absfin_mean)/R_inter)

def qp_PV_Base(componentSpecs,var):

    R_inter = componentSpecs["R_inter"]
    l_c = componentSpecs["l_c"]

    T_PV_m = var["T_PV"]
    T_Base_m = var["T_Base_mean"]
    L = componentSpecs["L_tube"]

    var["qp_PV_Base"] = l_c*((T_PV_m-T_Base_m)/R_inter)

def Qdot_absfin_back(componentSpecs,stepConditions,var):

    R_b = componentSpecs["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    L_af = componentSpecs["L_af"]

    T_absfin_m = var["T_absfin_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    var["Qdot_absfin_back"] = L*2*L_af*(T_absfin_m-T_back)/R_b

def Qdot_absfin_back_rad(componentSpecs,stepConditions,var):

    R_b = componentSpecs["R_2"] + 1/(var["h_rad_back"])
    L_af = componentSpecs["L_af"]

    T_absfin_m = var["T_absfin_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    var["Qdot_absfin_back_rad"] = L*2*L_af*(T_absfin_m-T_back)/R_b

def Qdot_absfin_back_conv(componentSpecs,stepConditions,var):
    
    R_b = componentSpecs["R_2"] + 1/(var["h_back"])
    L_af = componentSpecs["L_af"]

    T_absfin_m = var["T_absfin_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    var["Qdot_absfin_back_conv"] = L*2*L_af*(T_absfin_m-T_back)/R_b

def Qdot_absfin_tube(componentSpecs,stepConditions,var):

    T_absfin_m = var["T_absfin_mean"]
    T_tube_m = var["T_tube_mean"]
    L = componentSpecs["L_tube"]

    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    h_rad_tube_abs = var["h_rad_tube_abs"]

    var["Qdot_absfin_tube"] = L*p_ext_tube_rad*h_rad_tube_abs*(T_absfin_m-T_tube_m)

def Qdot_tube_back(componentSpecs,stepConditions,var):

    T_tube_m = var["T_tube_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    var["Qdot_tube_back"] = L*gamma_back*(T_tube_m - T_back)

def Qdot_f01(componentSpecs,stepConditions,var):

    L = componentSpecs["L_tube"]

    if componentSpecs["fin_0"]==1:
        gamma_0_int = var["gamma_0_int"]
    else:
        gamma_0_int = 0
    if componentSpecs["fin_1"]==1:
        gamma_1_int = var["gamma_1_int"]
    else:
        gamma_1_int = 0

    gamma = gamma_0_int + gamma_1_int
    
    T_tube_m = var["T_tube_mean"]
    T_back = stepConditions["T_back"]

    Q = L*gamma*(T_tube_m-T_back)

    var["Qdot_f01"] = Q

def Qdot_tube_back_wo_ins_conv(componentSpecs,stepConditions,var):
    T_tube_m = var["T_tube_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    h_back_tube = var["h_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    gamma_back = p_ext_tube*h_back_tube
    # gamma_0_int = var["gamma_0_int"]
    # gamma_1_int = var["gamma_1_int"]
    # gamma = gamma_back + gamma_0_int + gamma_1_int
    gamma = gamma_back

    var["Qdot_tube_back_conv"] = L*gamma*(T_tube_m - T_back)

def Qdot_tube_back_wo_ins_rad(componentSpecs,stepConditions,var):
    T_tube_m = var["T_tube_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    h_rad_back_tube = var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    gamma_back = p_ext_tube*h_rad_back_tube

    var["Qdot_tube_back_rad"] = L*gamma_back*(T_tube_m - T_back)

def Qdot_ins_tube_back_conv(componentSpecs,stepConditions,var):

    T_ins_tube_m = var["T_ins_tube_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    h_back_tube = var["h_back_tube"]

    var["Qdot_ins_tube_back_conv"] = L*p_ext_tube*h_back_tube*(T_ins_tube_m - T_back)

def Qdot_ins_tube_back_rad(componentSpecs,stepConditions,var):

    T_ins_tube_m = var["T_ins_tube_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    h_rad_back_tube = var["h_rad_back_tube"]

    var["Qdot_ins_tube_back_rad"] = L*p_ext_tube*h_rad_back_tube*(T_ins_tube_m - T_back)

def Qdot_ins_absfin_back_conv(componentSpecs,stepConditions,var):

    T_ins_absfin_m = var["T_ins_absfin_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]
    L_af = componentSpecs["L_af"]
    h_back = var["h_back"]

    var["Qdot_ins_absfin_back_conv"] = L*2*L_af*h_back*(T_ins_absfin_m - T_back)

def Qdot_ins_absfin_back_rad(componentSpecs,stepConditions,var):

    T_ins_absfin_m = var["T_ins_absfin_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]
    L_af = componentSpecs["L_af"]
    h_back = var["h_rad_back"]

    var["Qdot_ins_absfin_back_rad"] = L*2*L_af*h_back*(T_ins_absfin_m - T_back)

def Qdot_ins_conv(componentSpecs,stepConditions,var):

    T_ins_m = var["T_ins_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    h_back = var["h_back"]

    var["Qdot_ins_conv"] = L*h_back*(T_ins_m - T_back)

def Qdot_ins_rad(componentSpecs,stepConditions,var):

    T_ins_m = var["T_ins_mean"]
    T_back = stepConditions["T_back"]
    L = componentSpecs["L_tube"]

    h_back = var["h_rad_back"]

    var["Qdot_ins_rad"] = L*h_back*(T_ins_m - T_back)

def Qdot_tube_fluid(componentSpecs,stepConditions,var):

    h_fluid = var["h_fluid"]
    p_int_tube = componentSpecs["p_int_tube"]

    chi = 1/(h_fluid*p_int_tube)

    T_tube_m = var["T_tube_mean"]
    T_fluid_m = var["T_fluid_mean"]
    L = componentSpecs["L_tube"]

    var["Qdot_tube_fluid"]=(L/chi)*(T_tube_m-T_fluid_m)

def Qdot_Base_tube(componentSpecs,var):
    
    C_B = componentSpecs["C_B"]
    
    T_Base_m = var["T_Base_mean"]
    T_tube_m = var["T_tube_mean"]

    l_B = componentSpecs["l_B"]
    L = componentSpecs["L_tube"]

    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    h_rad_tube_abs = var["h_rad_tube_abs"]


    var["Qdot_Base_tube"] = L*(T_Base_m-T_tube_m)*(C_B+h_rad_tube_abs*p_ext_tube_rad)

def Qdot_Base_back(componentSpecs,stepConditions,var):

    R_b = componentSpecs["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    iota = componentSpecs["iota"]

    T_Base_m = var["T_Base_mean"]
    T_back = stepConditions["T_back"]

    L = componentSpecs["L_tube"]

    var["Qdot_Base_back"] = L*iota*(T_Base_m-T_back)/R_b

def qp_Base_back(componentSpecs,stepConditions,var):

    R_b = componentSpecs["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    iota = componentSpecs["iota"]

    T_Base_m = var["T_Base_mean"]
    T_back = stepConditions["T_back"]

    L = componentSpecs["L_tube"]

    var["qp_Base_back"] = iota*(T_Base_m-T_back)/R_b

def Qdot_absfin_Base(componentSpecs,var):
    q = var["qp_fin"]
    L = componentSpecs["L_tube"]

    var["Qdot_absfin_Base"] = 2*L*q

def Qdot_abs_back2(componentSpecs,var):
    var["Qdot_abs_back2"] = var["Qdot_absfin_Base"] - var["Qdot_PV_plate"] + var["Qdot_PV_Base"]

def power_balance_3(componentSpecs,var):
    Qdot_PV_Base = var["Qdot_PV_Base"]
    Qdot_absfin_Base = var["Qdot_absfin_Base"]
    Qdot_fluid = var["Qdot_tube_fluid"]
    Qdot_Base_back = var["Qdot_Base_back"]
    Qdot_fluid_back = var["Qdot_fluid_back"]

    var["power_balance_3"] = Qdot_PV_Base + Qdot_absfin_Base - (Qdot_fluid + Qdot_fluid_back) - Qdot_Base_back

def PB_3(componentSpecs,var):
    PB3 = var["qp_PV_Base"] - var["qp_Base_back"] + 2*var["qp_fin"]-var["qp_fluid"]
    var["PB_3"] = PB3
    # print(PB3)

def Qdot_fluid_back(componentSpecs,stepConditions,var):

    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]

    k = componentSpecs["k_ail"]
    h_fluid = var["h_fluid"]
    p_int_tube = componentSpecs["p_int_tube"]

    chi = 1/(h_fluid*p_int_tube)

    L = componentSpecs["L_tube"]

    if componentSpecs["fin_0"]==1:
        gamma_0_int = var["gamma_0_int"]
    else:
        gamma_0_int = 0
    if componentSpecs["fin_1"]==1:
        gamma_1_int = var["gamma_1_int"]
    else:
        gamma_1_int = 0

    R_2 = componentSpecs["R_2"]
    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)

    T_fluid_m = var["T_fluid_mean"]
    T_back = stepConditions["T_back"]

    zeta = (gamma_back)/(1+chi*(gamma_back+gamma_1_int+gamma_0_int))

    var["Qdot_fluid_back"] = L*zeta*(T_fluid_m-T_back)

def qp_f0(componentSpecs,stepConditions,var):

    T_fluid_m = var["T_fluid_mean"]
    T_back = stepConditions["T_back"]

    gamma_0_int = var["gamma_0_int"]

    var["qp_f0"] = gamma_0_int*(T_fluid_m-T_back)

def qp_f1(componentSpecs,stepConditions,var):

    T_fluid_m = var["T_fluid_mean"]
    T_back = stepConditions["T_back"]

    gamma_1_int = var["gamma_1_int"]

    var["qp_f1"] = gamma_1_int*(T_fluid_m-T_back)

def qp_f2(componentSpecs,stepConditions,var):

    T_abs_m = var["T_abs_mean"]
    T_back = stepConditions["T_back"]

    gamma_2_int = var["gamma_2_int"]

    var["qp_f2"] = gamma_2_int*(T_abs_m-T_back)

def Qdot_f2(componentSpecs,var):

    var["Qdot_f2"] = var["qp_f2"] * componentSpecs["L_tube"]
