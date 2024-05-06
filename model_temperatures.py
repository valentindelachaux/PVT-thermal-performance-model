import math
import copy
import pandas as pd
import numpy as np
import heat_transfer as bht
import ht

import model_ht as mht

from CoolProp.CoolProp import PropsSI

# Eq. 560.36
def T_fluid_out(componentSpecs,var):
    """Calculates the fluid outlet temperature and stores it in var["T_fluid_out"]
    
    $$
    T_{fluid,out} = (T_{fluid,in}+\frac{b_f}{a_f})e^{a_fL_{tube}} - \frac{b_f}{a_f}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        T_fluid_in (float): fluid inlet temperature
        var (dict): dictionary containing the variables
        
    Returns:
        None
    """

    a = var["a_f"]
    b = var["b_f"]

    L_tube = componentSpecs["L_tube"]

    T_fluid_in = var['T_fluid_in']

    res = (T_fluid_in+(b/a))*math.exp(a*L_tube) - b/a
    var["T_fluid_out"] = res

# Eq. 560.40
def T_fluid_mean(componentSpecs,var):
    """Calculates the mean fluid temperature and stores it in var["T_fluid_mean"]
    
    $$
    T_{fluid,mean} = \frac{T_{fluid,in}+\frac{b_f}{a_f}}{a_fL_{tube}}e^{a_fL_{tube}} - \frac{b_f}{a_fL_{tube}}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        T_fluid_in (float): fluid inlet temperature
        var (dict): dictionary containing the variables
    
    Returns:
        None"""

    L_tube = componentSpecs["L_tube"]

    T_fluid_in = var['T_fluid_in']

    h_back_tube = var["h_back_tube"]
    if h_back_tube == None:
        print(var["T_tube_mean"])
        h_back_tube = 3.

    a = var["a_f"]
    b = var["b_f"]

    res = ((T_fluid_in+(b/a))/(a*L_tube))*math.exp(a*L_tube) - (T_fluid_in+(b/a))/(a*L_tube) - b/a
    var["T_fluid_mean"] = res

def T_Base_mean(componentSpecs, stepConditions,var): #T_fluid has already been used for q_f_p and T_f_mean calculations
    """Calculates the mean base temperature and stores it in var["T_Base_mean"]
    
    $$
    T_{Base,mean} = \frac{q_{tube,fluid}}{C_B+h_{rad,f}p_{ext,tube,rad}} + T_{fluid,mean}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        
    Returns:
        None"""
    # C_B = componentSpecs["C_B"]

    # p_int_tube = componentSpecs["p_int_tube"]
    # h_fluid = var["h_fluid"]
    # chi = 1/(h_fluid*p_int_tube)

    # h_back = var["h_back"]+var["h_rad_back"]
    # p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    # R_2 = componentSpecs["R_2"]
    # gamma_back = p_ext_tube/(R_2+1/h_back)
    # gamma_0_int = var["gamma_0_int"]
    # gamma_1_int = var["gamma_1_int"]
    # gamma = gamma_back + gamma_0_int + gamma_1_int

    # h_rad_f = var["h_rad_f"]

    # e0 = var["e0"]
    T_fluid = var["T_fluid_mean"]
    q_tube_fluid = var["q_tube_fluid"]
    T_back = stepConditions["T_back"]
    T_sky = stepConditions["T_sky"]

    # res = (1/(C_B+h_rad_f))*(e0*q_tube_fluid - (1/chi)*T_fluid - gamma*T_back)
    # var["T_Base_mean"] = res

    # res = (1/(h_fluid*p_int_tube)+1/C_B)*q_tube_fluid + T_f_mean

    b1 = var["b1"]
    b2 = var["b2"]
    b3 = var["b3"]
    b4 = var["b4"]

    res = b1*q_tube_fluid + b2*T_fluid + b3*T_back + b4*T_sky
    var["T_Base_mean"] = res

def T_B_check(componentSpecs,stepConditions,var):

    b1 = var["b1"]
    b2 = var["b2"]
    b3 = var["b3"]
    b4 = var["b4"]

    var['T_B_check'] = b1*var['Qdot_tube_fluid']/componentSpecs['L_tube'] + b2*var['T_fluid_mean'] + b3*stepConditions['T_back'] + b4*stepConditions['T_sky']

# Eq. 560.42 -> calculate the mean fin temperature
def T_absfin_mean(componentSpecs,stepConditions,var):
    """Calculates the mean fin temperature and stores it in var["T_absfin_mean"]
    
    $$
    T_{abs,fin,mean} = \frac{b}{j} + \left(T_{Base,mean}-\frac{b}{j}\right)\frac{\tanh(mL_{af})}{mL_{af}}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        
    Returns:
        None"""
    W = componentSpecs["W"]
    L_af = componentSpecs["L_af"]

    S = var["S"]

    b = var["b"]
    j = var["j"]
    m = var["m"]

    T_B_mean = var["T_Base_mean"]

    var["T_absfin_mean"] = b/j + (T_B_mean-(b/j))*math.tanh(m*L_af)/(m*L_af)

# Eq. 560.43 -> calculate the mean absorber temperature
def T_abs_mean(componentSpecs,var):
    """Calculates the mean absorber temperature and stores it in var["T_abs_mean"]
    
    $$
    T_{abs,mean} = \frac{l_B T_{Base,mean}+2L_{af}T_{abs,fin,mean}}{W}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
    
    Returns:
        None
    """
    W = componentSpecs["W"]
    l_B = componentSpecs["l_B"]
    L_af = componentSpecs["L_af"]

    T_Base_mean = var["T_Base_mean"]
    T_absfin_mean = var["T_absfin_mean"]

    res = (l_B*T_Base_mean+(L_af*2)*T_absfin_mean)/W
    var["T_abs_mean"] = res

    # if stepConditions["compt"] >= 1:
    #     T_abs_mean_old = var["T_abs_mean"]
    #     var["T_abs_mean"] = (res+T_abs_mean_old)/2

def T_tube_mean(componentSpecs,stepConditions,var):
    """Calculates the mean tube temperature and stores it in var["T_tube_mean"]
    
    $$
    T_{tube,mean} = e_1T_{Base,mean} + e_2T_{fluid,mean} + e_3T_{back} + e_4T_{sky}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    
    Returns:
        None"""
    e1 = var["e1"]
    e2 = var["e2"]
    e3 = var["e3"]
    e4 = var["e4"]

    T_B = var["T_Base_mean"]
    T_fluid = var["T_fluid_mean"]
    T_back = stepConditions["T_back"]
    T_sky = stepConditions["T_sky"]

    var["T_tube_mean"] = e1*T_B + e2*T_fluid + e3*T_back + e4*T_sky

def T_glass_mean(componentSpecs,stepConditions,var):
    """Calculates the mean glass temperature and stores it in var["T_glass_mean"]
    
    $$
    T_{glass,mean} = \frac{1}{\alpha_g}\left(\alpha_gT_{PV} + h_{top,g}T_{amb} + h_{rad,g}T_{sky} + \frac{T_{PV}}{R_g}\right)
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        
    Returns:
        None"""
    alpha = componentSpecs["alpha_g"]

    h_top_g = var["h_top_g"]
    h_rad_g = var["h_rad_g"]
    R_g = componentSpecs["R_g"]
    G = stepConditions["G"]
    T_amb = stepConditions["T_amb"]
    T_sky = stepConditions["T_sky"]
    T_PV = var["T_PV"]

    res = (1/(h_top_g+h_rad_g+(1/R_g)))*(alpha*G + h_top_g * T_amb + h_rad_g * T_sky + (1/R_g)*T_PV)

    var["T_glass"] = res

def T_PV_mean(componentSpecs,stepConditions,var):
    """Calculates the mean PV temperature and stores it in var["T_PV_mean"]
    
    $$
    T_{PV,mean} = \frac{1}{\alpha_{PV}}\left(\alpha_{PV}T_{abs,mean} + h_{rad}T_{sky} + \frac{T_{abs,mean}}{R_{inter}}\right)
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    
    Returns:
        None"""
    
    R_inter = componentSpecs["R_inter"]
    T_sky = stepConditions["T_sky"]
    
    h_rad = var["h_rad"]
    S = var["S"]
    Fprime = var["Fp"]
    a2 = var["a2"]
    T_abs_mean = var["T_abs_mean"]

    res = Fprime*R_inter*(S-a2+h_rad*T_sky+(T_abs_mean/R_inter))

    var["T_PV0"] = var["T_PV"]
    var["T_PV"] = res

def T_PV_Base_mean(componentSpecs,stepConditions,var):

    R_inter = componentSpecs["R_inter"]
    T_sky = stepConditions["T_sky"]
    
    h_rad = var["h_rad"]
    S = var["S"]
    Fprime = var["Fp"]
    a2 = var["a2"]

    T_Base_mean = var["T_Base_mean"]

    res = Fprime*R_inter*(S-a2+h_rad*T_sky+(T_Base_mean/R_inter))

    var["T_PV_Base_mean"] = res

def T_PV_absfin_mean(componentSpecs,var):
    L_af = componentSpecs["L_af"]
    W = componentSpecs["W"]

    T_PV_Base_mean = var["T_PV_Base_mean"]
    T_PV_mean = var["T_PV"]

    var["T_PV_absfin_mean"] = (W*T_PV_mean - (W-2*L_af)*T_PV_Base_mean)/(2*L_af)

def T_ins_tube_mean(componentSpecs,var):
    R_2 = componentSpecs["R_2"]
    T_tube_mean = var["T_tube_mean"]
    Qdot_tube_back = var["Qdot_tube_back"]

    L = componentSpecs["L_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"]

    var["T_ins_tube_mean"] = T_tube_mean - R_2*Qdot_tube_back/(L*p_ext_tube)

def T_ins_absfin_mean(componentSpecs,var):
    R_2 = componentSpecs["R_2"]
    T_absfin_mean = var["T_absfin_mean"]
    Qdot_absfin_back = var["Qdot_absfin_back"]

    L = componentSpecs["L_tube"]
    L_af = componentSpecs["L_af"]

    var["T_ins_absfin_mean"] = T_absfin_mean - R_2*Qdot_absfin_back/(L*2*L_af)

def T_ins_mean(componentSpecs,var):

    T_ins_tube_mean = var["T_ins_tube_mean"]
    T_ins_absfin_mean = var["T_ins_absfin_mean"]
    l_B = componentSpecs["l_B"]
    lambd_riser_back = componentSpecs["lambd_riser_back"]
    L_af = componentSpecs["L_af"]
    W = componentSpecs["W"]

    var["T_ins_mean"] = (l_B*T_ins_tube_mean + 2*L_af*T_ins_absfin_mean)/W
