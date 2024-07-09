import math
import copy
import pandas as pd
import numpy as np
import heat_transfer as bht
import ht

import model_ht as mht
import model_temperatures as mtemp

from CoolProp.CoolProp import PropsSI

import model_transfers as mtr

mean_list = ["T_glass","T_PV","T_PV_Base_mean","T_PV_absfin_mean","T_abs_mean","T_Base_mean","T_absfin_mean","T_ins_mean","T_ins_tube_mean","T_ins_absfin_mean","T_tube_mean","T_fluid_mean","h_top_g","h_rad","h_back","h_rad_back","h_back_tube","h_rad_back_tube","h_back_fins","h_rad_tube_abs","h_fluid","X_celltemp","eta_PV","S"]
add_list = ["Qdot_sun_glass","Qdot_sun_PV","Qdot_top_conv","Qdot_top_rad","Qdot_glass_PV","Qdot_PV_sky","Qdot_PV_plate","Qdot_PV_Base","Qdot_PV_absfin","Qdot_absfin_Base","Qdot_absfin_back","Qdot_absfin_back_conv","Qdot_absfin_back_rad","Qdot_Base_tube","Qdot_Base_back","Qdot_tube_sky","Qdot_tube_fluid","Qdot_tube_back","Qdot_ins_tube_back_conv","Qdot_ins_tube_back_rad","Qdot_ins_absfin_back_conv","Qdot_ins_absfin_back_rad","Qdot_tube_back_conv","Qdot_tube_back_rad","Qdot_absfin_back","Qdot_f01"]

# Iteration solving functions

def a0(componentSpecs,stepConditions,var):
    """Calculates the a0 factor and stores it in var["a0"]
    
    $$a_0 = \frac{1}{h_{top}+h_{rad}+1/R_g}\left(\alpha_g G + h_{top} T_{amb} + h_{rad} T_{sky}\right)$$

    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    
    Returns:
        None
    """

    var["a0"] = (1/(var["h_top_g"]+var["h_rad_g"]+1/componentSpecs["R_g"]))*(componentSpecs["alpha_g"]*stepConditions["G"] + var["h_top_g"]*stepConditions["T_amb"] + var["h_rad_g"]*stepConditions["T_sky"])

def a1(componentSpecs,stepConditions,var):
    """Calculates the a1 factor and stores it in var["a1"]
    
    $$a_1 = \frac{1}{h_{top}+h_{rad}+1/R_g}\left(\frac{1}{R_g}\right)$$
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables

    Returns:
        None
    """

    var["a1"] = (1/(var["h_top_g"]+var["h_rad_g"]+1/componentSpecs["R_g"]))*(1/componentSpecs["R_g"])

def a3(componentSpecs,stepConditions,var):
    """Calculates the a3 factor and stores it in var["a3"]
    
    $$a_3 = (1-a_1) \frac{1}{R_g}$$
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    
    Returns:
        None
    """

    a1 = var["a1"]
    var["a3"] = (1-a1)*(1/componentSpecs["R_g"])

def a2(componentSpecs,stepConditions,var):
    """Calculates the a2 factor and stores it in var["a2"]
    
    $$a_2 = -\frac{a_0}{R_g}$$
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    Returns:
        None
    """

    var["a2"] = - var["a0"]/componentSpecs["R_g"]

def X_celltemp(componentSpecs,var):
    """Calculates the X_celltemp factor and stores it in var["X_celltemp"]
    
    $$X_{celltemp} = 1+\Eff_{T}(T_{PV}-T_{ref})$$
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None
    """
    Eff_T = componentSpecs["Eff_T"]
    T_ref = componentSpecs["T_ref"]

    T_PV = var["T_PV"]

    X = 1+Eff_T*(T_PV-T_ref)

    var["X_celltemp"]=X

def eta_PV(componentSpecs,stepConditions,var):
    """Calculates the PV efficiency and stores it in var["eta_PV"]
    
    $$
    \eta_{PV} = \eta_{nom}X_{celltemp}X_{rad}X_{corr}
    $$

    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    
    Returns:
        None"""
    
    eta_nom = componentSpecs["eta_nom"]
    G = stepConditions["G"]
    X_rad = componentSpecs["X_rad"]
    X_corr = componentSpecs["X_corr"]

    X_celltemp = var["X_celltemp"]

    eta = eta_nom*X_celltemp*X_rad*X_corr
    var["eta_PV"] = eta

def S(componentSpecs,stepConditions,var):
    """Calculates the net absorbed solar radiation (total absorbed - PV power production) and stores it in var["S"]
    
    $$
    S = \tau_{\theta}G(1-\eta_{PV})
    $$

    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        
    Returns:
        None"""
    
    tau_g = componentSpecs["tau_g"]
    G = stepConditions["G"]

    #T_PV = var["T_PV"]
    eta_PV = var["eta_PV"]

    S = tau_g*G*(1-eta_PV)

    var["S"] = S

def S_star(componentSpecs,stepConditions,var):
    """Calculates the $S^{*}$ factor and stores it in var["S_star"]
    
    $$S^{*} = S - a_2 + h_{rad}T_{sky}$$

    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    Returns:
        None
    """

    var["S_star"] = var["S"] - var["a2"] + var["h_rad"]*stepConditions["T_sky"]

def Fp(componentSpecs, var):
    """Calculates the Fp factor and stores it in var["Fp"]
    
    $$F_p = \frac{1}{a_3 R_{inter}+h_{rad}R_{inter}+1}$$
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters    
        var (dict): dictionary containing the variables
        
    Returns:
        None"""


    R_inter = componentSpecs["R_inter"]

    h_rad = var["h_rad"]
    a3 = var["a3"]

    Fp = 1/(a3*R_inter+h_rad*R_inter+1)
    var["Fp"] = Fp

def j(componentSpecs,var):
    """Calculates the j factor and stores it in var["j"]
    
    $$
    j = \frac{1}{R_{inter}F'}+\frac{1}{R_bF'}-\frac{1}{R_{inter}}+\frac{h_{rad,f}}{F'}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""
    
    R_inter = componentSpecs["R_inter"]
    R_b = componentSpecs["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    h_rad_tube_abs = var["h_rad_tube_abs"]

    Fprime = var["Fp"]

    j = 1/(Fprime*R_b) + 1/(R_inter*Fprime) - 1/R_inter

    if componentSpecs["fin_2"] == 1:

        gamma_int = var["gamma_2_int"]

        j += (gamma_int)/Fprime

    var["j"] = j

def m(componentSpecs, var):
    """Calculates the m factor and stores it in var["m"]
    
    $$
    m = \sqrt{\frac{F'j}{k_{abs}\lambda_{abs}}}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    lambd_abs = componentSpecs["lambd_abs"]
    k_abs = componentSpecs["k_abs"]

    Fprime = var["Fp"]

    j = var["j"]

    m = math.sqrt((Fprime*j)/(k_abs*lambd_abs))

    var["m"] = m

def b(componentSpecs,stepConditions,var):
    """Calculates the b factor and stores it in var["b"]
    
    $$
    b = S+h_{rad}T_{sky}+\frac{T_{amb}}{R_t}+\frac{T_{back}}{R_bF'}+\frac{h_{rad,f}T_{tube,mean}}{F'}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    T_back = stepConditions["T_back"]
    R_b = componentSpecs["R_2"] + 1/(var["h_back"]+var["h_rad_back"])

    h_rad = var["h_rad"]
    S = var["S"]
    Fprime = var["Fp"]

    h_rad_tube_abs = var["h_rad_tube_abs"]

    # if h_rad_tube_abs != 0:
    #     T_tube_mean = var["T_tube_mean"]
    #     b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime) + (h_rad_tube_abs*T_tube_mean)/Fprime
    # else:
    #     b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime)

    # b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime)

    S_star = var["S_star"]

    b = S_star + T_back/(R_b*Fprime)


    if componentSpecs["fin_2"]==1:
        gamma_int = var["gamma_2_int"]

        b += (gamma_int*T_back)/Fprime

    var["b"] = b

def calc_gamma(componentSpecs, var):

    h_back_tube = var["h_back_tube"]
    h_rad_back_tube = var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]

    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_rad_back = p_ext_tube_rad/(R_2+1/h_rad_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]

    gamma = gamma_back + gamma_rad_back + gamma_0_int + gamma_1_int

    return gamma

def e0(componentSpecs, var):
    """Calculates the e0 factor and stores it in var["e0"]
    
    $$
    c_0 = \frac{1}{R_{inter}}+\frac{1}{R_bF'}+\frac{h_{rad,f}}{F'}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
    
    Returns:
        None"""

    p_int_tube = componentSpecs["p_int_tube"]
    h_fluid = var["h_fluid"]
    C_B = componentSpecs["C_B"]   
    chi = 1/(h_fluid*p_int_tube)
 
    gamma = calc_gamma(componentSpecs, var)
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    h_rad_tube_abs = var["h_rad_tube_abs"]

    h_rad_tube_sky = var["h_rad_tube_sky"]
    p_tube_sky = componentSpecs["p_tube_sky"]

    # vérifier l'homogénéité
    var["e0"] = 1/chi + h_rad_tube_abs*p_ext_tube_rad + C_B + gamma + h_rad_tube_sky*p_tube_sky

def e1(componentSpecs,var):
    """Calculates the e1 factor and stores it in var["e1"]
    
    $$
    e_1 = \frac{C_B+p_{ext,tube,rad}h_{rad,f}}{c_0}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    C_B = componentSpecs["C_B"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    h_rad_tube_abs = var["h_rad_tube_abs"]

    e0 = var["e0"]

    var["e1"] = (1/e0)*(C_B+p_ext_tube_rad*h_rad_tube_abs)

def e2(componentSpecs,var):
    """Calculates the e2 factor and stores it in var["e2"]

    $$
    e_2 = \frac{1}{c_0\chi}
    $$

    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables 

    Returns:
        None"""

    p_int_tube = componentSpecs["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    e0 = var["e0"]

    var["e2"] = (1/e0)*(1/chi)

def e3(componentSpecs,var):
    """Calculates the e3 factor and stores it in var["e3"]
    
    $$
    e_3 = \frac{\gamma}{c_0}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    gamma = calc_gamma(componentSpecs, var)

    e0 = var["e0"]

    var["e3"] = (1/e0)*gamma

def e4(componentSpecs,var):
    """Calculates the e4 factor and stores it in var["e4"]
    
    $$e_4 = \frac{1}{c_0}h_{rad,tube,sky}$$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None
    """

    e0 = var["e0"]

    h_rad_tube_sky = var["h_rad_tube_sky"]
    p_tube_sky = componentSpecs["p_tube_sky"]

    var["e4"] = (1/e0)*h_rad_tube_sky*p_tube_sky

def f0(componentSpecs,var):
    """Calculates the f0 factor and stores it in var["f0"]
    
    $$
    f_0 = C_B + h_{rad,f}p_{ext,tube,rad} + \gamma
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""
    
    C_B = componentSpecs["C_B"]

    h_rad_tube_abs = var["h_rad_tube_abs"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    h_rad_tube_sky = var["h_rad_tube_sky"]
    p_tube_sky = componentSpecs["p_tube_sky"]

    gamma = calc_gamma(componentSpecs, var)

    var["f0"] = C_B + h_rad_tube_abs*p_ext_tube_rad + gamma + h_rad_tube_sky*p_tube_sky

def b1(componentSpecs, var):
    """Calculates the b1 factor and stores it in var["b1"]
    
    $$
    b_1 = \frac{1}{C_B + h_{rad,f}p_{ext,tube,rad}}\frac{1}{1-c_2}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    C_B = componentSpecs["C_B"]  
    h_rad_tube_abs = var["h_rad_tube_abs"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    e1 = var["e1"]
    f0 = var["f0"]
    # print('b1')
    # print('e1',e1)
    # print('gamma',gamma)

    var["b1"] = 1/(C_B + h_rad_tube_abs*p_ext_tube_rad - f0*e1)
                   
def b2(componentSpecs, var):
    """Calculates the b2 factor and stores it in var["b2"]
    
    $$
    b_2 = \frac{1}{c_0\chi}\frac{1}{1-c_2}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    f0 = var["f0"]
    e2 = var["e2"]
    b1 = var["b1"]

    var["b2"] = f0*e2*b1

def b3(componentSpecs, var):
    """Calculates the b3 factor and stores it in var["b3"]
    
    $$
    b_3 = -\frac{1}{C_B + h_{rad,f}p_{ext,tube,rad}}\gamma
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    C_B = componentSpecs["C_B"]  
    h_rad_tube_abs = var["h_rad_tube_abs"]

    gamma = calc_gamma(componentSpecs, var)

    e3 = var["e3"]
    f0 = var["f0"]
    # c'est vérifié c'est égal
    # var["b3"] = (-1/(C_B + h_rad_tube_abs*p_ext_tube_rad))*gamma
    # var["b3"] =  (-1/(C_B + h_rad_tube_abs*p_ext_tube_rad))*gamma

    var["b3"] = (f0*e3 - gamma)*var["b1"]

def b4(componentSpecs, var):
    """Calculates the b4 factor and stores it in var["b4"]
    
    $$b_4 = b_1(f_0e_4 - h_{rad,tube,sky})$$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
    
    Returns:
        None
    """

    b1 = var["b1"]
    C_B = componentSpecs["C_B"]
    
    gamma = calc_gamma(componentSpecs, var)

    h_rad_tube_abs = var["h_rad_tube_abs"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    h_rad_tube_sky = var["h_rad_tube_sky"]
    p_tube_sky = componentSpecs["p_tube_sky"]

    e4 = var["e4"]
    f0 = var["f0"]

    var["b4"] = b1*(f0*e4 - h_rad_tube_sky*p_tube_sky)

def d1(componentSpecs, var):
    """Calculates the d1 factor and stores it in var["d1"]
    
    $$
    d_1 = \frac{1}{c_0\chi}\frac{1}{1-c_2}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    d0 = 1/var["e0"]

    p_int_tube = componentSpecs["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    gamma = calc_gamma(componentSpecs, var)

    h_rad_tube_abs = var["h_rad_tube_abs"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    e1 = var["e1"]

    # var["d1"] = (-gamma*d0)/R_B + h_rad_tube_abs*(1-(d0/R_B))
    var["d1"] = -e1*(gamma+h_rad_tube_abs*p_ext_tube_rad) + h_rad_tube_abs*p_ext_tube_rad

def d2(componentSpecs, var):
    """Calculates the d2 factor and stores it in var["d2"]
    
    $$
    d_2 = \frac{1}{c_0\chi}\frac{1}{1-c_2}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    d0 = 1/var["e0"]

    p_int_tube = componentSpecs["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    gamma = calc_gamma(componentSpecs, var)
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    h_rad_tube_abs = var["h_rad_tube_abs"]

    # var["d2"] = (-gamma*d0)/chi - (h_rad_tube_abs*d0)/chi

    e2 = var["e2"]

    var["d2"] = -e2*(gamma+h_rad_tube_abs*p_ext_tube_rad)

def d3(componentSpecs, var):
    """Calculates the d3 factor and stores it in var["d3"]
    
    $$
    d_3 = \frac{1}{c_0\chi}\frac{1}{1-c_2}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    d0 = 1/var["e0"]

    p_int_tube = componentSpecs["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    gamma = calc_gamma(componentSpecs, var)
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    h_rad_tube_abs = var["h_rad_tube_abs"]

    # var["d3"] = -gamma**2*d0 + gamma - (h_rad_tube_abs*d0)/(1/gamma)

    e3 = var["e3"]
    var["d3"] = -e3*(gamma+h_rad_tube_abs*p_ext_tube_rad) + gamma

def d4(componentSpecs, var):
    """Calculates the d4 factor and stores it in var["d4"]
    
    $$d_4 = \frac{1}{c_0\chi}\frac{1}{1-c_2}$$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    gamma = calc_gamma(componentSpecs, var)
    p_ext_tube_rad  = componentSpecs["p_ext_tube_rad"]

    h_rad_tube_abs = var["h_rad_tube_abs"]

    # var["d3"] = -gamma**2*d0 + gamma - (h_rad_tube_abs*d0)/(1/gamma)

    e4 = var["e4"]

    var["d4"] = -e4*(gamma + h_rad_tube_abs*p_ext_tube_rad)

def KTE_Bt(componentSpecs,stepConditions,var):
    """Calculates Ka_Bt, Th_Bt, and Ep_Bt factors and stores them in var["Ka_Bt"], var["Th_Bt"], and var["Ep_Bt"]
    
    $$
    g_1 = (l_c/R_{inter})*(F'-1)*b_1
    g_2 = (l_c/R_{inter})*(F'-1)*b_2
    g_3 = (l_c/R_{inter})*(F'-1)*b_3
    g_4 = (l_c/R_{inter})*(F'-1)*b_4
    g_5 = l_c*F'*S^{*}

    h_1 = (-iota/R_b)*b_1
    h_2 = (-iota/R_b)*b_2
    h_3 = (iota/R_b)*(1-b_3)
    h_4 = (-iota/R_b)*b_4
    h_5 = 0

    i_1 = -2k_{abs}\lambda_{abs}m\tanh(mL_{af})(b_1/j)
    i_2 = -2k_{abs}\lambda_{abs}m\tanh(mL_{af})b_2
    i_3 = -2k_{abs}\lambda_{abs}m\tanh(mL_{af})b_3
    i_4 = -2k_{abs}\lambda_{abs}m\tanh(mL_{af})b_4
    i_5 = 2k_{abs}\lambda_{abs}m\tanh(mL_{af})(b/j)

    \kappa_{Bt} = g_2+h_2+i_2
    \theta_{Bt} = -(g_1+h_1+i_1)
    \epsilon_{Bt} = (g_3+h_3+i_3)*T_{back} + (g_4+h_4+i_4)*T_{sky} + (g_5+h_5+i_5)  
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables

    Returns:
        None
    """
    lambd_abs = componentSpecs["lambd_abs"]
    k_abs = componentSpecs["k_abs"]
    W = componentSpecs["W"]
    L_af = componentSpecs["L_af"]
    l_B = componentSpecs["l_B"]
    l_c = componentSpecs["l_c"]
    p_int_tube = componentSpecs["p_int_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    R_inter = componentSpecs["R_inter"]

    R_b = componentSpecs["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    h_fluid = var["h_fluid"]

    T_sky = stepConditions["T_sky"]
    T_amb = stepConditions["T_amb"]
    T_back = stepConditions["T_back"]

    C_B = componentSpecs["C_B"]

    #T_PV = var["T_PV"]
    h_rad = var["h_rad"]
    S = var["S"]
    Fprime = var["Fp"]

    j = var["j"]
    b = var["b"]
    m = var["m"] 
    
    iota = componentSpecs["iota"]

    # K = b2 * ( -D_tube*Fprime*((l_c/D_tube)*(h_rad+1/R_t)+(iota/D_tube)/(R_b*Fprime))-2*k_abs*lambd_abs*m*math.tanh(m*L_af) )
    # T = b1 * ( D_tube*Fprime*((l_c/D_tube)*(h_rad+1/R_t)+(iota/D_tube)/(R_b*Fprime))+2*k_abs*lambd_abs*m*math.tanh(m*L_af))
    # E = D_tube*Fprime*((l_c/D_tube)*(S+h_rad*T_sky+T_amb/R_t)+((iota/D_tube)*(1-b3)*T_back)/(R_b*Fprime))+2*k_abs*lambd_abs*m*math.tanh(m*L_af)*((b/j) - b3*T_back) + (l_c/R_t)*(Fprime-1)*b3*T_back

    b1 = var["b1"]
    b2 = var["b2"]
    b3 = var["b3"]
    b4 = var["b4"]

    # T_PV - T_B 
    # if stepConditions["compt"] <=2:
    #     g4 = l_c*Fprime*(S+h_rad*T_sky+T_amb/R_t)
    #     g3 = (l_c/R_inter)*(Fprime-1)*b3
    #     g2 = (l_c/R_inter)*(Fprime-1)*b2
    #     g1 = (l_c/R_inter)*(Fprime-1)*b1
    # else:
    #     T_abs = var["T_abs_mean"]
    #     g4 = l_c*(S+h_rad*T_sky+T_amb/R_t + T_abs/R_inter)
    #     g3 = -(l_c/R_inter)*b3
    #     g2 = -(l_c/R_inter)*b2
    #     g1 = -(l_c/R_inter)*b1

    S_star = var["S_star"]


    g5 = l_c*Fprime*S_star
    g4 = (l_c/R_inter)*(Fprime-1)*b4
    g3 = (l_c/R_inter)*(Fprime-1)*b3
    g2 = (l_c/R_inter)*(Fprime-1)*b2
    g1 = (l_c/R_inter)*(Fprime-1)*b1

    # T_B - T_back
    h5 = 0.
    h4 = (-iota/R_b)*b4
    h3 = (iota/R_b)*(1-b3)
    h2 = (-iota/R_b)*b2
    h1 = (-iota/R_b)*b1

    # 2q'_absfin
    i5 = 2*k_abs*lambd_abs*m*math.tanh(m*L_af)*(b/j)
    i4 = -2*k_abs*lambd_abs*m*math.tanh(m*L_af)*b4
    i3 = -2*k_abs*lambd_abs*m*math.tanh(m*L_af)*b3
    i2 = -2*k_abs*lambd_abs*m*math.tanh(m*L_af)*b2
    i1 = -2*k_abs*lambd_abs*m*math.tanh(m*L_af)*b1

    K =  g2+h2+i2
    T = -(g1+h1+i1)
    E = (g3+h3+i3)*T_back + (g4+h4+i4)*T_sky + (g5+h5+i5)

    var["Ka_Bt"] = K
    var["Th_Bt"] = T
    var["Ep_Bt"] = E 

def KTE_tf(componentSpecs,stepConditions,var):
    """Calculates Ka_tf, Th_tf, and Ep_tf factors and stores them in var["Ka_tf"], var["Th_tf"], and var["Ep_tf"]

    $$
    \kappa_{tf} = \frac{d_2}{\chi R_B} + \frac{d_0}{\chi^2} - \frac{1}{\chi}
    $$

    $$
    \theta_{tf} = 1 - \frac{d_1}{\chi R_B}
    $$

    $$
    \epsilon_{tf} = \frac{d_3}{\chi R_B}T_{back} + \frac{d_4}{\chi R_B}T_{sky} + \frac{d_1}{\chi R_B}T_{back} + \frac{d_0\gamma}{\chi}T_{back}
    $$

    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables

    Returns:
        None
    """

    Ka_Bt = var["Ka_Bt"]
    Th_Bt = var["Th_Bt"]
    Ep_Bt = var["Ep_Bt"]

    b1 = var["b1"]
    b2 = var["b2"]
    b3 = var["b3"]
    b4 = var["b4"]

    d0 = 1/var["e0"]
    d1 = var["d1"]
    d2 = var["d2"]
    d3 = var["d3"]
    d4 = var["d4"]

    C_B = componentSpecs["C_B"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    h_rad_tube_abs = var["h_rad_tube_abs"]

    R_B = 1/(C_B+p_ext_tube_rad*h_rad_tube_abs)
    p_int_tube = componentSpecs["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    gamma = calc_gamma(componentSpecs, var)

    T_sky = stepConditions["T_sky"]
    T_back = stepConditions["T_back"]

    # var["Ka_tf"] = Ka_Bt - (gamma*b2)/c3
    # var["Th_tf"] = 1/c3 + gamma*b1 + Th_Bt
    # var["Ep_tf"] = Ep_Bt + gamma*(1-b3)*T_back

    var["Ka_tf"] = d1*b2 + d2 + Ka_Bt
    var["Th_tf"] = 1 - d1*b1 + Th_Bt
    # calculé : puis test avec changement des signes ci-dessus
    #     var["Th_tf"] = 1 - d1*b1 + Th_Bt
    var["Ep_tf"] = Ep_Bt + (d1*b3 + d3) * T_back + (d1*b4 + d4) * T_sky 

    # var["Ka_tf"] = (d0*b2)/(chi*R_B) + (d0/chi**2) - 1/chi
    # var["Th_tf"] = 1 - (d0*b1)/(chi*R_B)
    # var["Ep_tf"] = ( (d0*b3)/(chi*R_B) + (d0*gamma)*chi )*T_back

    # tentative de correction
    var["Ka_tf"] = var["Ka_Bt"] - (var["e1"] * var["b2"] + var["e2"]) * (gamma + var["h_rad_tube_sky"] * componentSpecs["p_tube_sky"])
    var["Th_tf"] = 1 + var["Th_Bt"] + var["e1"] * var["b1"] * (gamma + var["h_rad_tube_sky"] * componentSpecs["p_tube_sky"])
    var["Ep_tf"] = var["Ep_Bt"] + (gamma - (var["e1"] * var["b3"] + var["e3"]) * (gamma + var["h_rad_tube_sky"] * componentSpecs["p_tube_sky"])) * stepConditions["T_back"] + (var["h_rad_tube_sky"] * componentSpecs["p_tube_sky"] - (var["e1"] * var["b4"] + var["e4"]) * (gamma + var["h_rad_tube_sky"] * componentSpecs["p_tube_sky"])) * stepConditions["T_sky"]

def ab_f(componentSpecs,stepConditions,var):
    """Calculates the a_f and b_f factors and stores them in var["a_f"] and var["b_f"]

    $$
    a_f = \frac{N_{harp}}{\dot{m}C_p}\frac{\kappa_{tf}}{\theta_{tf}}
    $$

    $$
    b_f = \frac{N_{harp}}{\dot{m}C_p}\frac{\epsilon_{tf}}{\theta_{tf}}
    $$

    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
    
    Returns:
        None
    """
    N_harp = componentSpecs["N_harp"]
    mdot = stepConditions["mdot"]
    Cp = var["Cp"]
    
    Ka_tf = var["Ka_tf"]
    Th_tf = var["Th_tf"]
    Ep_tf = var["Ep_tf"]

    a = (N_harp/(mdot*Cp))*(Ka_tf/Th_tf)
    b = (N_harp/(mdot*Cp))*(Ep_tf/Th_tf)

    var["a_f"] = a
    var["b_f"] = b

def Cp(componentSpecs,stepConditions,var,hyp):
    """Calculates the specific heat capacity of the fluid and stores it in var["Cp"]
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None
    """
    T_m = (var["T_fluid_in"]+var["T_fluid_out"])/2 # K

    p_fluid = hyp["p_fluid"] # bar
    fluid = hyp["fluid"]
    glycol_rate = hyp["glycol_rate"] # %

    var["Cp"] = PropsSI('C','P', p_fluid*100000, 'T', T_m, f'INCOMP::{fluid}[{glycol_rate}]')

def Biot(lambd,k,h,delta):
    """Calculates the Biot number
    
    Args:
        lambd (float): thickness of the material [m]
        k (float): thermal conductivity of the material [W/m/K]
        h (float): heat transfer coefficient [W/m2/K]
        delta (float): width [m]
        
    Returns:
        float: Biot number"""
    return ((lambd*h)/k)*(1+lambd/delta)

def Bi_f3(componentSpecs,var):
    h_back = var["h_back_fins"]
    var["Bi_f3"] = Biot(componentSpecs["lambd_ail"],componentSpecs["k_ail"],h_back,componentSpecs["delta_f3"])

def gamma0int(N,L_fin,lambd,k,delta,delta_int,L_tube,h):
    """Calculates the gamma_0_int factor and returns it
    
    $$
    \gamma_{0,int} = \frac{\alpha}{\lambda}\frac{\sinh(\alpha L_{fin}/\lambda)+\beta\alpha/\lambda\cosh(\alpha L_{fin}/\lambda)}{\cosh(\alpha L_{fin}/\lambda)+\beta\sinh(\alpha L_{fin}/\lambda)}
    $$
    
    Args:
        N (int): number of fins
        L_fin (float): length of the fin
        lambd (float): thickness of the fin
        k (float): thermal conductivity of the fin
        delta (float): width of the fin
        delta_int (float): contact length between the fin and the tube
        L_tube (float): length of the tube
        h (float): heat transfer coefficient
        
    Returns:
        Bi (float): Biot number
        gamma_0_int (float): gamma_0_int factor"""
    
    Bi = Biot(lambd,k,h,delta)

    alpha = math.sqrt(2*Bi)
    beta = math.sqrt(Bi/2)*(1/(1+lambd/delta))
    arg = (alpha*L_fin)/lambd

    numerateur = (alpha/lambd)*math.sinh(arg) + ((beta*alpha)/lambd)*math.cosh(arg)
    denominateur = math.cosh(arg) + beta*math.sinh(arg)

    return Bi,k*(numerateur/denominateur)*((lambd*N*delta_int)/L_tube)

def gamma_0_int(componentSpecs,var):
    """Calculates the gamma_0_int factor and stores it in var["gamma_0_int"]"""

    var["Bi_f0"],var["gamma_0_int"] = gamma0int(componentSpecs["N_f0"],componentSpecs["L_f0"],componentSpecs["lambd_ail"],componentSpecs["k_ail"],componentSpecs["delta_f0"],componentSpecs["delta_f0_int"],componentSpecs["L_tube"],var["h_back_fins"])

def gamma1int(N,L_fin,lambd,k,delta,delta_int,L_tube,h):
    """Calculates the gamma_1_int factor and returns it
    
    $$
    \gamma_{1,int} = \frac{\alpha}{\lambda}\frac{\sinh(\alpha L_{fin}/\lambda)}{\cosh(\alpha L_{fin}/\lambda)}
    $$
    
    Args:
        N (int): number of fins
        L_fin (float): length of the fin
        lambd (float): thickness of the fin
        k (float): thermal conductivity of the fin
        delta (float): width of the fin
        delta_int (float): contact length between the fin and the tube
        L_tube (float): length of the tube
        h (float): heat transfer coefficient
        
    Returns:
        Bi (float): Biot number
        gamma_1_int (float): gamma_1_int factor"""
    
    Bi = Biot(lambd,k,h,delta)

    return Bi,2*k*((lambd*N*delta_int)/L_tube)*math.tanh(math.sqrt(2*Bi)*(L_fin/lambd))*(math.sqrt(2*Bi)/lambd)

def gamma_1_int(componentSpecs,var):

    var["Bi_f1"],var["gamma_1_int"] = componentSpecs["coeff_f1"]*gamma1int(componentSpecs["N_f1"],componentSpecs["L_f1"],componentSpecs["lambd_ail"],componentSpecs["k_ail"],componentSpecs["delta_f1"],componentSpecs["delta_f1_int"],componentSpecs["L_tube"],var["h_back_fins"])

def gamma_2_int(componentSpecs,var):

    Bi = Biot(componentSpecs["lambd_ail"],componentSpecs["k_ail"],var["h_back_fins"],componentSpecs["delta_f2"])
    var["Bi_f2"] = Bi
    a = componentSpecs["lambd_ail"]
    delta = componentSpecs["delta_f2"]

    alpha = math.sqrt(2*Bi)
    beta = math.sqrt(Bi/2)*(1/(1+a/delta))

    L_a = componentSpecs["L_f2"]
    N_ail = componentSpecs["N_f2"]
    k = componentSpecs["k_ail"]

    L_tube = componentSpecs["L_tube"]

    arg = (alpha*L_a)/a
    numerateur = (alpha/a)*math.sinh(arg) + ((beta*alpha)/a)*math.cosh(arg)
    denominateur = math.cosh(arg) + beta*math.sinh(arg)

    delta_f2 = componentSpecs["delta_f2"]
    
    gamma_int = k*(numerateur/denominateur)*((a*N_ail*delta_f2)/(L_tube*delta_f2))

    var["gamma_2_int"] = gamma_int

# Main functions

def one_loop(componentSpecs,stepConditions,var,hyp):
    """Procedure which calculates the variables for one loop
    
    \enumerate{
        \item gamma_{0,int}, gamma_{1,int}, gamma_{2,int}, Bi_{f3}
        \item h_{rad,g}, h_{rad}, h_{top,mean}
        \item a_0, a_1, a_3, a_2
        \item X_{celltemp}, eta_{PV}, S, and S^{*}
        \item F', j, m, b
        \item c_0, c_2, e_1, e_2, e_3, e_4, f_0, b_1, b_2, b_3, b_4
        \item \kappa_{Bt}, \theta_{Bt}, \epsilon_{Bt}
        \item d_1, d_2, d_3, d_4
        \item \kappa_{tf}, \theta_{tf}, \epsilon_{tf}
        \item a_f and b_f
        \item T_{fluid,out}
        \item q_{tube-fluid}
        \item \overline{T_{fluid}}
        \item \overline{T_{Base}}
        \item \overline{T_{tube}}
        \item \overline{T_{absfin}}
        \item \overline{T_{abs}}
        \item Qdot_{tube-back}
        \item Qdot_{absfin-back}
        \item \overline{T_{ins-tube}}
        \item \overline{T_{ins-absfin}}
        \item \overline{T_{ins}}
        \item Qdot_{ins,conv}
        \item Qdot_{ins,rad}
        \item h_{back,mean}, h_{rad,back}, h_{back,tube}, h_{rad,tube-sky}, h_{rad,back-tube}, h_{back,fins}, h_{rad,f}
        \item \overline{T_{PV}}, \overline{T_{PV,Base}}, \overline{T_{PV,absfin}}
        \item mtemp.T_glass_mean()
        \item mtr.qp_PV_Base()
        \item mtr.qp_Base_back()
        \item mtr.qp_fin()
        \item Cp()
    }

    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        T_fluid_in (float): fluid inlet temperature
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypothesis
    
    Returns:
        None"""

    if componentSpecs["fin_0"] == 1:
        gamma_0_int(componentSpecs,var)
    else:
        var["gamma_0_int"] = 0
    if componentSpecs["fin_1"] == 1:
        gamma_1_int(componentSpecs,var)
    else:
        var["gamma_1_int"] = 0
    if componentSpecs["fin_2"] == 1:
        gamma_2_int(componentSpecs,var)
    else:
        var["gamma_2_int"] = 0
    if componentSpecs["fin_3"] == 1:
        Bi_f3(componentSpecs,var)
    else:
        pass
    
    mht.h_rad_g(componentSpecs,stepConditions,var,hyp)
    mht.h_rad(componentSpecs,stepConditions,var,hyp)

    if componentSpecs["fin_0"] == 1 or componentSpecs["fin_1"] == 1 or componentSpecs["fin_2"] == 1:
        mht.h_top_mean(componentSpecs,stepConditions,var,hyp)
    else:
        mht.h_top_g(componentSpecs,stepConditions,var,hyp)

    a0(componentSpecs,stepConditions,var)
    a1(componentSpecs,stepConditions,var)
    a3(componentSpecs,stepConditions,var)
    a2(componentSpecs,stepConditions,var)


    X_celltemp(componentSpecs,var)
    eta_PV(componentSpecs,stepConditions,var)
    S(componentSpecs,stepConditions,var)
    S_star(componentSpecs,stepConditions,var)
    Fp(componentSpecs,var)
    j(componentSpecs,var)
    m(componentSpecs,var)
    b(componentSpecs,stepConditions,var)

    e0(componentSpecs,var)
    e1(componentSpecs,var)
    e2(componentSpecs,var)
    e3(componentSpecs,var)
    e4(componentSpecs,var)

    f0(componentSpecs,var)
    b1(componentSpecs,var)
    b2(componentSpecs,var)
    b3(componentSpecs,var)
    b4(componentSpecs,var)

    KTE_Bt(componentSpecs,stepConditions,var)

    d1(componentSpecs,var)
    d2(componentSpecs,var)
    d3(componentSpecs,var)
    d4(componentSpecs,var)

    KTE_tf(componentSpecs,stepConditions,var)

    ab_f(componentSpecs,stepConditions,var)
    mtemp.T_fluid_out(componentSpecs,var)
    mtr.q_tube_fluid(componentSpecs,stepConditions,var)
    mtemp.T_fluid_mean(componentSpecs,var)
    mtemp.T_Base_mean(componentSpecs,stepConditions,var)
    mtemp.T_tube_mean(componentSpecs,stepConditions,var)
    mtemp.T_absfin_mean(componentSpecs,stepConditions,var)
    mtemp.T_abs_mean(componentSpecs,var)


    mtr.Qdot_tube_back(componentSpecs,stepConditions,var)
    mtr.Qdot_absfin_back(componentSpecs,stepConditions,var)   
    mtemp.T_ins_tube_mean(componentSpecs,var)
    mtemp.T_ins_absfin_mean(componentSpecs,var)
    mtemp.T_ins_mean(componentSpecs,var)

    mtr.Qdot_ins_conv(componentSpecs,stepConditions,var)
    mtr.Qdot_ins_rad(componentSpecs,stepConditions,var)

    if hyp["calc_h_back_mean"]==1:
        mht.h_back_mean(componentSpecs,stepConditions,var,hyp)
    else:
        mht.h_back_abs(componentSpecs,stepConditions,var,hyp)

    mht.h_rad_back(componentSpecs,stepConditions,var,hyp)
    mht.h_back_tube(componentSpecs,stepConditions,var,hyp)
    mht.h_rad_tube_sky(componentSpecs,stepConditions,var,hyp)
    mht.h_rad_back_tube(componentSpecs,stepConditions,var,hyp)
    mht.h_back_fins(componentSpecs,stepConditions,var,hyp)

    mht.h_rad_tube_abs(componentSpecs,stepConditions,var,hyp)

    mtemp.T_PV_mean(componentSpecs,stepConditions,var)
    mtemp.T_PV_Base_mean(componentSpecs,stepConditions,var)
    mtemp.T_PV_absfin_mean(componentSpecs,var)
    mtemp.T_glass_mean(componentSpecs,stepConditions,var)

    mtr.qp_PV_Base(componentSpecs,var)
    mtr.qp_Base_back(componentSpecs,stepConditions,var)
    mtr.qp_fin(componentSpecs,var)

    Cp(componentSpecs,stepConditions,var,hyp)

def compute_power(componentSpecs,stepConditions,var):
    mtr.Qdot_top_conv(componentSpecs,stepConditions,var)
    mtr.Qdot_top_rad(componentSpecs,stepConditions,var)
    mtr.Qdot_sun_glass(componentSpecs,stepConditions,var)
    mtr.Qdot_sun_PV(componentSpecs,stepConditions,var)
    mtr.Qdot_glass_PV(componentSpecs,stepConditions,var)
    mtr.Qdot_PV_sky(componentSpecs,stepConditions,var)
    mtr.Qdot_PV_plate(componentSpecs,var)
    # Qdot_abs_back1(componentSpecs,stepConditions,var)
    mtr.Qdot_PV_Base(componentSpecs,var)
    mtr.Qdot_PV_absfin(componentSpecs,var)
    mtr.Qdot_Base_back(componentSpecs,stepConditions,var)
    mtr.Qdot_tube_fluid(componentSpecs,stepConditions,var)
    mtr.Qdot_tube_back(componentSpecs,stepConditions,var)
    mtr.Qdot_tube_sky(componentSpecs,stepConditions,var)
    mtr.Qdot_Base_tube(componentSpecs,var)
    # qp_fluid_back(componentSpecs,var)
    mtr.qp_fin(componentSpecs,var)
    mtr.Qdot_absfin_Base(componentSpecs,var)
    mtr.Qdot_abs_back2(componentSpecs,var)
    mtr.Qdot_fluid_back(componentSpecs,stepConditions,var)
    mtr.Qdot_tube_back(componentSpecs,stepConditions,var)
    mtr.Qdot_absfin_back(componentSpecs,stepConditions,var)
    mtr.Qdot_absfin_back_conv(componentSpecs,stepConditions,var)
    mtr.Qdot_absfin_back_rad(componentSpecs,stepConditions,var)
    # mtr.Qdot_absfin_tube(componentSpecs,stepConditions,var)
    mtr.Qdot_tube_back_wo_ins_conv(componentSpecs,stepConditions,var)
    mtr.Qdot_tube_back_wo_ins_rad(componentSpecs,stepConditions,var)
    mtr.Qdot_ins_tube_back_conv(componentSpecs,stepConditions,var)
    mtr.Qdot_ins_tube_back_rad(componentSpecs,stepConditions,var)
    mtr.Qdot_ins_absfin_back_conv(componentSpecs,stepConditions,var)
    mtr.Qdot_ins_absfin_back_rad(componentSpecs,stepConditions,var)

    if componentSpecs["fin_0"]==1 or componentSpecs["fin_1"]==1:
        mtr.Qdot_f01(componentSpecs,stepConditions,var)
    else:
        var["Qdot_f01"] = 0.

    mtr.power_balance_1(componentSpecs,var)
    mtr.power_balance_3(componentSpecs,var)

    if componentSpecs["fin_0"] == 1:
        mtr.qp_f0(componentSpecs,stepConditions,var)
    if componentSpecs["fin_1"] == 1:
        mtr.qp_f1(componentSpecs,stepConditions,var)
    if componentSpecs["fin_2"] == 1:
        pass

    mtemp.T_B_check(componentSpecs,stepConditions,var)

def initialize_var(var,componentSpecs,stepConditions,hyp,i):
    # Initialize the var dictionary with all necessary keys and values
    
    var = {'Slice' : i,
           'T_PV0':0,
           'Cp': hyp['Cp0'],
            'h_rad_tube_abs':hyp['h_rad_tube_abs0'],
            }
    
    if i==0:
        var= {**var,**{
            "h_top_g": hyp["h_top0"],
            "h_back": hyp["h_back0"],
            "h_rad_back": hyp["h_rad_back0"],
            "h_back_tube": hyp["h_back_tube0"],
            "h_rad_tube_sky": hyp["h_rad_tube_sky0"],
            "h_rad_back_tube": hyp["h_rad_back_tube0"],
            "h_back_fins": hyp["h_back_fins0"],

            'T_PV':stepConditions["guess_T_PV"],
            'T_fluid_in':stepConditions["T_fluid_in0"],
            'T_fluid_out':stepConditions["T_fluid_in0"]
        }}
    else:
        var = {**var,**{
            'h_top_g':var['h_top_g'],
            'h_back':var['h_back'],
            'h_rad_back':var['h_rad_back'],
            'h_rad_tube_sky':var['h_rad_tube_sky'],
            'h_back_tube':var['h_back_tube'],
            'h_back_fins':var['h_back_fins'],
            'h_rad_back_tube':var['h_rad_back_tube'],

            'T_PV':var['T_PV'],
            'T_fluid_in':var['T_fluid_out'],
            'T_fluid_out':var['T_fluid_out']
            }}

    var["T_tube_mean"] = (var["T_PV"]+var["T_fluid_in"])/2
    var["T_glass"] = var["T_PV"]
    mht.h_fluid(componentSpecs,stepConditions,var,hyp)

    return var

def simu_one_steady_state(componentSpecs, stepConditions, hyp):
    """Procedure which calculates the variables for one steady state
    
    \enumerate{
        \item Initialize the var dictionary with all necessary keys and values
        \item Loop until convergence is reached
        \item Update heat transfer coefficients
        \item Calculate the power balance
        \item Update the fluid inlet temperature
        \item Append the current iteration to the list of iterations
        }
        
    Args:
        componentSpecs (dict): dictionary containing the parameters
        stepConditions (dict): dictionary containing the meteo inputs
        hyp (dict): dictionary containing the hypothesis
        
    Returns:
        slices_df (DataFrame): dataframe containing the variables for each slice of the panel
        df_one (DataFrame): dataframe containing the variables for the last iteration
        its_data_list (list): list of dataframes containing the variables for each iteration
    """

    N_meander = componentSpecs["N_meander"]
    N_harp = componentSpecs["N_harp"]
    slices_data = []
    its_data_list = []

    var = {}

    for i in range(N_meander):
        var = initialize_var(var, componentSpecs, stepConditions, hyp, i)
        compt = 0
        stepConditions["compt"] = compt
        its_data = []

        while compt <= 3 or abs(var["T_PV"] - var["T_PV0"]) >= 0.00001:
            compt += 1
            stepConditions["compt"] = compt

            one_loop(componentSpecs, stepConditions, var, hyp)
            compute_power(componentSpecs, stepConditions, var)

            row_data = {**{'G': stepConditions["G"], 'Gp': stepConditions["Gp"],'T_amb': stepConditions["T_amb"], 'T_sky': stepConditions['T_sky'], 'T_back':stepConditions['T_back'], 'T_back_rad':stepConditions['T_back_rad'], 'u': stepConditions['u'],
                           'mdot': stepConditions['mdot'],'T_fluid_in':stepConditions['T_fluid_in0']}, **var}
            its_data.append(row_data)

        # Append only once after the loop
        its_data_list.append(pd.DataFrame(its_data))

        slices_data.append(its_data[-1])

    slices_df = pd.DataFrame(slices_data)
    df_mean = slices_df.mean()
    df_sum = slices_df.sum()
    df_one = pd.DataFrame()

    for key in slices_df.keys():
        if key in ['mdot','G','Gp','T_amb','T_sky','T_back','T_back_rad','u']:
            df_one[key] = [stepConditions[key]]
        elif key == "T_fluid_in":
            df_one[key] = [stepConditions["T_fluid_in0"]]
        elif key == "T_fluid_out":
            df_one[key] = [slices_df["T_fluid_out"].iloc[-1]]
        elif key in mean_list:
            df_one[key] = [df_mean[key]]
        elif key in add_list:
            df_one[key] = [df_sum[key]*N_harp]

    return slices_df, df_one, its_data_list

def simu_one_steady_state_all_he(panelSpecs,stepConditions,hyp, method_anomaly = 0):
    
    res = {}
    save_stepConditions = stepConditions.copy()
    save_T_fluid_in0 = stepConditions["T_fluid_in0"]

    # Test the main part

    panelSpecs['main']['name'] = 'part1'
    slices_df, df_one, its_data_list = simu_one_steady_state(panelSpecs['main'],stepConditions,hyp)
    res['main'] = {'slices_df':slices_df.copy(),'df_one':df_one.copy(),'its_data_list':its_data_list.copy()}

    hyp['h_back_prev'] = df_one['h_back'].values[0] # h_back de l'absorbeur
    hyp['h_top_man'] = df_one["h_top_g"].values[0]

    if len(panelSpecs['decomp'])>1:
        decomp = 1

        for part,part_name in list(panelSpecs['decomp'].items())[1:]:

            panelSpecs[part]['name'] = part

            if method_anomaly== 1 and (panelSpecs[part]['is_anomaly'] == 1 or panelSpecs[part]['is_inlet_man'] == 1 or panelSpecs[part]['is_outlet_man'] == 1):
                hyp['method_h_back_abs'] = 'free'
                
            slices_df, df_one, its_data_list = simu_one_steady_state(panelSpecs[part],stepConditions,hyp)
            res[part] = {'slices_df':slices_df.copy(),'df_one':df_one.copy(),'its_data_list':its_data_list.copy()}

            stepConditions["T_fluid_in0"] = df_one["T_fluid_out"].values[0]

    else:
        decomp = 0

    df_one = pd.DataFrame()

    for measure in res["main"]['df_one'].keys():
        if measure in ['mdot','G','Gp','T_amb','T_sky','T_back','T_back_rad','u']:
            df_one[measure] = [stepConditions[measure]]
        elif measure == "T_fluid_in":
            df_one[measure] = [save_T_fluid_in0]
        elif measure == "T_fluid_out":
            if decomp == 1:
                last_past = list(panelSpecs['decomp'].keys())[-1]
                df_one[measure] = [res[last_past]['df_one']['T_fluid_out'].values[0]]
            else:
                df_one[measure] = [res['main']['df_one']['T_fluid_out'].values[0]]
        elif measure in mean_list:
            av = 0
            Aire_tot = 0
            if decomp == 1:
                for part in res.keys():
                    if part == 'main':
                        continue
                    Aire = panelSpecs[part]["N_harp"]*panelSpecs[part]["W"]*panelSpecs[part]["L_tube"]
                    Aire_tot += Aire
                    av += res[part]['df_one'][measure].values[0]*Aire
                df_one[measure] = [av/Aire_tot]
            else:
                df_one[measure] = [res['main']['df_one'][measure].values[0]]
        elif measure in add_list:
            sum = 0
            if decomp == 1:
                for part in res.keys():
                    if part == 'main':
                        continue
                    sum += res[part]['df_one'][measure].values[0]
                df_one[measure] = [sum]
            else:
                df_one[measure] = [res['main']['df_one'][measure].values[0]]
    
    stepConditions = save_stepConditions.copy()

    return df_one,res

def simu_steadyStateConditions(panelSpecs,hyp,steadyStateConditions_df):
    
    # Dataframe object pour la liste des résultats sur tous les points de fonctionnement
    df_res = pd.DataFrame()

    compt_test = 0

    list_df = []
    list_list_df_historic = []
    list_res = []

    for i in range(0,len(steadyStateConditions_df)):

        stepConditions = {'G':steadyStateConditions_df["G"][i],"T_amb":steadyStateConditions_df["T_amb"][i],"T_back":steadyStateConditions_df["T_amb"][i],"u":steadyStateConditions_df["u"][i], "u_back" : steadyStateConditions_df["u_back"][i], "T_fluid_in0":steadyStateConditions_df["T_fluid_in"][i]}
        change_T_sky(stepConditions,hyp,'TUV')  # calculate Gp and T_sky

        stepConditions["mdot"] = steadyStateConditions_df["mdot"][i]

        # stepConditions["guess_T_PV"] = stepConditions["T_amb"] - 25
        stepConditions["guess_T_PV"] = (stepConditions["T_amb"]+stepConditions["T_fluid_in0"])/2

        df_one,res = simu_one_steady_state_all_he(panelSpecs,stepConditions,hyp)

        df_res = pd.concat([df_res,df_one],ignore_index=True)
        list_res.append(res)

        compt_test+=1

    # Analysing df

    # Be careful here you have zeros for some columns

    tab = pd.DataFrame()

    df_res['DT'] = df_res['T_fluid_out'] - df_res['T_fluid_in']
    df_res['Tm'] = (df_res['T_fluid_out'] + df_res['T_fluid_in'])/2
    df_res['T_m en °C'] = df_res['Tm']-273.15

    tab['G'] = df_res['G'] # a0
    tab['-(T_m - T_a)'] = -(df_res['Tm'] - df_res['T_amb']) # a1
    # tab['-(T_m - T_a)^2'] = -(df_res['Tm'] - df_res['T_amb'])**2 # a2
    tab['-(T_m - T_a)^2'] = 0.*df_res['Tm'] # a2
    tab['-up x (T_m - T_a)'] = (df_res['u'] - 3) * tab['-(T_m - T_a)'] # a3
    # tab['Gp'] = df_res['Gp'] # a4
    tab['Gp'] = 0. * df_res['Gp'] # a4
    tab['dTm/dt'] = 0. * df_res['Gp'] # a5
    tab['up x G'] = -(df_res['u'] - 3) * df_res['G'] # a6
    tab['up x Gp'] = -(df_res['u'] - 3) * df_res["Gp"] # a7
    tab['-(T_m - T_a)^4'] = -tab['-(T_m - T_a)']**4 # a8

    # coeff_density = [999.85,0.05332,-0.007564,0.00004323,-1.673e-7,2.447e-10]
    # coeff_density = list(reversed(coeff_density))

    coeff_Cp = [4.2184,-0.0028218,0.000073478,-9.4712e-7,7.2869e-9,-2.8098e-11,4.4008e-14]
    coeff_Cp = list(reversed(coeff_Cp))

    # df_res['density(T)'] = np.polyval(coeff_density,df_res['T_m en °C'])
    df_res['Cp(T)'] = np.polyval(coeff_Cp,df_res['T_m en °C'])*1000

    # df_res['mdot'] = df_res['density(T)']*(stepConditions["mdot"]/1000)

    df_res['Qdot'] = df_res['mdot']*df_res['Cp(T)']*df_res['DT']
    df_res['Qdot / AG'] = df_res['Qdot']/(panelSpecs['AG'])

    matrice = tab.to_numpy()
    B = df_res['Qdot / AG'].to_numpy()

    X = np.linalg.lstsq(matrice, B, rcond = -1)

    #_ = plt.plot(df['T_m*'].to_numpy(), B, 'o', label='Original data', markersize=2)
    #_ = plt.plot(df['T_m*'].to_numpy(), np.dot(matrice,X[0]), 'o', label='Fitted line',markersize=2)
    #_ = plt.legend()
    #plt.show()

    df_res_to_concat = df_res.drop(columns=["G","Gp"])

    df_res = pd.concat([tab,df_res_to_concat],axis=1)

    return df_res,X,list_res

def update_heat_transfer_coefficients(var, hyp, index):
    # This function should update heat transfer coefficients for the current index
    # based on prior values or initial guesses
    # Replace '...' with actual logic to calculate or update the coefficients
    if index == 0:
        ...
    else:
        ...

def recap_energy_balances(df_one):

    balances = {}
    balances_perc_of_average = {}

    balances['glass'] = df_one['Qdot_sun_glass'] - df_one['Qdot_top_conv'] - df_one['Qdot_top_rad'] - df_one['Qdot_glass_PV']
    balances['PV'] = df_one['Qdot_sun_PV'] + df_one['Qdot_glass_PV'] - df_one['Qdot_PV_sky'] - df_one['Qdot_PV_Base'] - df_one['Qdot_PV_absfin']
    balances['Base'] = df_one['Qdot_PV_Base'] + df_one['Qdot_absfin_Base'] - df_one['Qdot_Base_tube'] - df_one['Qdot_Base_back']
    balances['absfin'] = df_one['Qdot_PV_absfin'] - df_one['Qdot_absfin_Base'] - df_one['Qdot_absfin_back_conv'] - df_one['Qdot_absfin_back_rad'] 
    balances['tube'] = df_one['Qdot_Base_tube'] - df_one['Qdot_tube_fluid'] - df_one['Qdot_tube_back_conv'] - df_one['Qdot_tube_back_rad'] - df_one['Qdot_tube_sky']

    # list all the df_one columns in the 5 rows above in a list
    ht_labels = ['Qdot_sun_glass','Qdot_top_conv','Qdot_top_rad','Qdot_glass_PV','Qdot_sun_PV','Qdot_PV_sky','Qdot_PV_Base','Qdot_PV_absfin','Qdot_absfin_Base','Qdot_Base_tube','Qdot_Base_back','Qdot_tube_fluid','Qdot_tube_back_conv','Qdot_tube_back_rad','Qdot_tube_sky','Qdot_absfin_back_conv','Qdot_absfin_back_rad']
    ht_average = sum(abs(df_one[key].values[0]) for key in ht_labels)/len([key for key in ht_labels if df_one[key].values[0] != 0])

    for key in balances:
        balances_perc_of_average[key] = balances[key]/ht_average

    return balances, balances_perc_of_average

def recap_residuals(panelSpecs, df_one, res):

    # Initialize a dictionary to store the residuals
    residuals = {'Part': [], 'glass': [], 'PV': [], 'Base': [], 'absfin': [], 'tube': []}

    # Loop through the parts and calculate balances
    for part in panelSpecs['decomp'].keys():
        bal, bal_perc = recap_energy_balances(res[part]['df_one'])
        residuals['Part'].append(part)
        for key, value in bal_perc.items():
            residuals[key].append(value[0])

    # Calculate total balances
    bal, bal_perc = recap_energy_balances(df_one)
    residuals['Part'].append('Total')
    for key, value in bal_perc.items():
        residuals[key].append(value[0])

    # Create DataFrame
    df_residuals = pd.DataFrame(residuals)

    pd.set_option('display.float_format', lambda x: f'{x:.2e}')
    print(df_residuals)
    pd.reset_option(pat='display.float_format')

    return df_residuals

# Other functions

def find_a_i(df,componentSpecs):
    tab = pd.DataFrame()

    df['DT'] = df['T_fluid_out'] - df['T_fluid_in']
    df['Tm'] = (df['T_fluid_out'] + df['T_fluid_in'])/2

    df['T_m en °C'] = df['Tm']-273.15

    tab['G'] = df['G'] # a0
    tab['-(T_m - T_a)'] = -(df['Tm'] - df['T_amb']) # a1
    # tab['-(T_m - T_a)^2'] = -(df['Tm'] - df['T_amb'])**2 # a2
    tab['-(T_m - T_a)^2'] = 0.*df['Tm'] # a2
    tab['-up x (T_m - T_a)'] = (df['u'] - 3) * tab['-(T_m - T_a)'] # a3
    # tab['Gp'] = df['Gp'] # a4
    tab['Gp'] = 0. * df['Gp'] # a4
    tab['dT/dt'] = 0. * df['Gp'] # a5 = 0
    tab['up x G'] = -(df['u'] - 3) * df['G'] # a6
    tab['up x Gp'] = -(df['u'] - 3) * df["Gp"] # a7
    tab['-(T_m - T_a)^4'] = -tab['-(T_m - T_a)']**4 # a8

    # coeff_density = [999.85,0.05332,-0.007564,0.00004323,-1.673e-7,2.447e-10]
    # coeff_density = list(reversed(coeff_density))

    coeff_Cp = [4.2184,-0.0028218,0.000073478,-9.4712e-7,7.2869e-9,-2.8098e-11,4.4008e-14]
    coeff_Cp = list(reversed(coeff_Cp))

    # df['density(T)'] = np.polyval(coeff_density,df['T_m en °C'])
    df['Cp(T)'] = np.polyval(coeff_Cp,df['T_m en °C'])*1000

    # df['mdot'] = df['density(T)']*(stepConditions["mdot"]/1000)

    df['Qdot'] = df['mdot']*df['Cp(T)']*df['DT']
    df['Qdot / AG'] = df['Qdot']/(componentSpecs['AG'])

    matrice = tab.to_numpy()
    B = df['Qdot / AG'].to_numpy()

    X = np.linalg.lstsq(matrice, B, rcond = -1)

    return X

def change_u(componentSpecs,stepConditions,wind_speed):
    stepConditions["u"] = wind_speed
    
    a_w = componentSpecs["a_htop"]
    b_w = componentSpecs["b_htop"]

    new_h_wind = a_w*wind_speed+b_w

    componentSpecs["h_top"]=new_h_wind
    componentSpecs["R_t"]= componentSpecs["R_top"] + 1/componentSpecs["h_top"]

def change_T_sky(stepConditions,hyp,type):
    if type == "TUV":
        stepConditions["Gp"] = 4
        stepConditions["T_sky"] = stepConditions["T_amb"]
        # stepConditions["T_sky"] = ((stepConditions["Gp"]/hyp["sigma"]) + stepConditions["T_amb"]**4)**(1/4)
    
    else :
        Tsk = 0.0552*stepConditions["T_amb"]**1.5

        stepConditions["T_sky"] = Tsk
        stepConditions["Gp"] = hyp["sigma"]*(stepConditions["T_sky"]**4 - stepConditions["T_amb"]**4)

def change_N_ail(componentSpecs,N):
    componentSpecs["N_ail"] = N

def change_air_layer(componentSpecs,lambd_air):
    old_air_layer = componentSpecs["lambd_air"]
    k_air = componentSpecs["k_air"]

    old_R_T = componentSpecs["R_inter"]

    old_r_air = old_air_layer/k_air
    new_r_air = lambd_air/k_air

    componentSpecs["R_inter"] = old_R_T - old_r_air + new_r_air
    componentSpecs["lambd_air"] = lambd_air
    #print(componentSpecs["R_inter"])

def change_b_htop(componentSpecs,stepConditions,b_htop):
    componentSpecs["b_htop"] = b_htop

    change_u(componentSpecs,stepConditions,componentSpecs["u"])

def change_ins(componentSpecs,e_new,k_new):

    componentSpecs["R_2"]=e_new/k_new

def change_N_fins_per_EP(componentSpecs,N):
    componentSpecs["N_fins_on_abs"] = (6*N)/componentSpecs["N_harp"]
    componentSpecs["D"] = (0.160/N)