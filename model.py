import math
import copy
import pandas as pd
import numpy as np
import heat_transfer as bht
import ht

from CoolProp.CoolProp import PropsSI

# %run C:\Users\BU05\Documents\Modele1D_Type560\Type560.py

# Function for Excel
# Search a cell by its nom_variable and return a string like "B3"
def find_cell_by_name(wb,nom_variable):
    my_range = wb.defined_names[nom_variable]
    ch = my_range.value
    ch2 = ch.split("!")[1]
    ch3 = ch2.replace("$","")
    return ch3


## Relationship between wind and R_top
# Modify the dictionary componentSpecs by updating the wind speed
# To complete from Excel file
def change_u(componentSpecs,stepConditions,wind_speed):
    stepConditions["u"] = wind_speed
    
    a_w = componentSpecs["a_htop"]
    b_w = componentSpecs["b_htop"]

    new_h_wind = a_w*wind_speed+b_w

    componentSpecs["h_top"]=new_h_wind
    componentSpecs["R_t"]= componentSpecs["R_top"] + 1/componentSpecs["h_top"]

# return tanh or 1/tanh
def tanh_or_inverse(arg):
    return math.tanh(arg)

def h_fluid(componentSpecs,stepConditions,var,hyp):
    """Calculates the convective heat transfer coefficient between the fluid and the tube wall and stores it in var["h_fluid"]
    
    Args:
        componentSpecs (dict): dictionary containing the parameters of the PVT panel
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None
        
    """

    D_tube = componentSpecs["D_tube"]
    L_tube = componentSpecs["L_tube"]

    if componentSpecs["is_inlet_man"] == 1 or componentSpecs["is_outlet_man"] == 1 and componentSpecs["geometry"] == "harp":
        mdot = stepConditions["mdot"]/2
    else:
        mdot = stepConditions["mdot"]
    
    N_harp = componentSpecs["N_harp"]

    T_fluid = stepConditions["T_fluid_in0"]

    p_fluid = hyp["p_fluid"]
    fluid = hyp["fluid"]
    glycol_rate = hyp["glycol_rate"]

    k_fluid = PropsSI('L', 'P', p_fluid, 'T', T_fluid, f'INCOMP::{fluid}[{glycol_rate}]')
    rho_fluid = PropsSI('D', 'P', p_fluid, 'T', T_fluid, f'INCOMP::{fluid}[{glycol_rate}]')
    mu_fluid = PropsSI('V', 'P', p_fluid, 'T', T_fluid, f'INCOMP::{fluid}[{glycol_rate}]')
    Pr = PropsSI('Prandtl', 'P', p_fluid, 'T', T_fluid, f'INCOMP::{fluid}[{glycol_rate}]')

    flow_rate_per_riser = (mdot/N_harp)/rho_fluid # en m3/s
    tube_section = math.pi*(D_tube/2)**2

    fluid_speed = flow_rate_per_riser/tube_section

    Re = (rho_fluid*fluid_speed*D_tube)/mu_fluid

    

    if Re < 2000:
        if componentSpecs["tube_geometry"] == "rectangular":
            Nu = ht.conv_internal.Nu_laminar_rectangular_Shan_London(min(componentSpecs["H_tube"],componentSpecs["w_tube"])/max(componentSpecs["H_tube"],componentSpecs["w_tube"]))
        else:
            Nu = ht.conv_internal.Nu_conv_internal(Re,Pr,Method='Laminar - constant Q')
    else:
        Nu = ht.conv_internal.turbulent_Colburn(Re,Pr)

    var["Re"] = Re
    var["Nu"] = Nu
    var["h_fluid"]  = (k_fluid/D_tube)*Nu

def h_top(componentSpecs,stepConditions,var,hyp):
    """Calculates the convective heat transfer coefficient between the top of the panel and the ambient air and stores it in var["h_top_g"]
    
    $$
    h_{top} = \left( h_{free}^3 + h_{forced}^3 \right)^{1/3}
    $$

    Args:
        componentSpecs (dict): dictionary containing the parameters of the PVT panel
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None"""

    # Manifold
    if componentSpecs["is_inlet_man"] == 1 or componentSpecs["is_outlet_man"] == 1:
        if hyp['method_h_top_g_manifold'] == 'like_exchanger':
            var['h_top_g'] = hyp['h_top_man']
        else:
            raise ValueError("Method for h_top_g is not well defined for manifolds")

    # Heat exchanger
    else:
        T_glass = var["T_glass"]
        T_amb = stepConditions["T_amb"]

        if componentSpecs["orientation"]=="portrait":
            L_c = componentSpecs["L_pan"]
        else:
            L_c = componentSpecs["w_pan"]

        h_free = bht.top_h_simple(T_glass,T_amb,hyp["theta"],L_c)
        h_forced_turbulent = bht.h_top_forced_turbulent(T_glass,T_amb,stepConditions["u"],L_c)
        h_forced = bht.h_top_forced(T_glass,T_amb,stepConditions["u"],L_c)
        h_custom = bht.h_top_custom(T_glass,T_amb,stepConditions["u"],L_c)

        if hyp['method_h_top_g_exchanger'] == 'free_with_coeff':
            var["h_top_g"] = hyp["coeff_h_top_free"]*h_free
        elif hyp['method_h_top_g_exchanger'] == 'free':
            var["h_top_g"] = h_free
        elif hyp['method_h_top_g_exchanger'] == 'forced_turbulent_with_coeff':
            h_forced = hyp["coeff_h_top_forced"]*h_forced_turbulent
        elif hyp['method_h_top_g_exchanger'] == 'forced_turbulent':
            var["h_top_g"] = h_forced_turbulent
        elif hyp['method_h_top_g_exchanger'] == 'forced_with_coeff':
            var["h_top_g"] = hyp["coeff_h_top_forced"]*h_forced
        elif hyp['method_h_top_g_exchanger'] == 'forced':
            var["h_top_g"] = h_forced
        elif hyp['method_h_top_g_exchanger'] == 'free_forced_with_coeff':
            if stepConditions['u'] < 0.1:
                var["h_top_g"] = hyp["coeff_h_top_free"]*h_free
            else:
                var["h_top_g"] = hyp["coeff_h_top_forced"]*h_forced
        elif hyp['method_h_top_g_exchanger'] == 'mixed_with_coeff':
            var["h_top_g"] = ( (hyp["coeff_h_top_free"]*h_free)**3 + (hyp["coeff_h_top_forced"]*h_forced)**3 )**(1/3)
        elif hyp['method_h_top_g_exchanger'] == 'mixed':
            var["h_top_g"] = (h_free**3 + h_forced**3)**(1/3)
        elif hyp['method_h_top_g_exchanger'] == 'custom':
            var["h_top_g"] = h_custom
        else:
            raise ValueError("Method for h_top_g is not well defined for exchanger")

def h_back(componentSpecs,stepConditions,var,hyp):
    """Calculates the convective heat transfer coefficient between the back of the panel and the ambient air and stores it var["h_back"]
    
    Args:
        componentSpecs (dict): dictionary containing the parameters of the PVT panel
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None"""

    # Anomaly
    if  componentSpecs["is_anomaly"] == 1:

        if hyp['method_h_back_anomaly'] == "like_exchanger":
            var["h_back"] = hyp['h_back_prev']
        else:
            raise ValueError("Method for h_back is not well defined for anomalies")

    # Manifold
    elif componentSpecs["is_inlet_man"] == 1 or componentSpecs["is_outlet_man"] == 1 :

        L_c = componentSpecs['H_tube']

        if componentSpecs["insulated"] == 1:
            T_ref = var["T_ins_mean"]    
        else:
            T_ref = var["T_abs_mean"]

        if hyp['method_h_back_manifold'] == "free_cylinder":

                # res = bht.back_h_mixed(T_ref,stepConditions["T_back"],stepConditions["u_back"],hyp["theta"],L_c)
            var["h_back"] = bht.back_h_cylinder(T_ref,stepConditions["T_back"],L_c)

        elif hyp['method_h_back_manifold'] == "like_exchanger":
            var["h_back"] = hyp['h_back_prev']

        else:
            raise ValueError("Method for h_back is not well defined for manifolds")

    # Heat exchanger 
    else:

        # Parameters

        if componentSpecs["orientation"] == 'portrait':
            L_c = componentSpecs["L_abs"]
        else:
            L_c = componentSpecs["w_abs"]

        # If error
        if var["T_abs_mean"]==None:
            print('T_abs_mean = None in h_back()')
            var["h_back"] = 0.5
        
        # Case with fins
        elif componentSpecs["fin_0"] >= 1 or componentSpecs["fin_1"] >= 1 or componentSpecs["fin_2"] >= 1:

            D = componentSpecs["D"]

            if componentSpecs["N_ail"]<= 24:
                var["h_back"] = hyp["coeff_h_back"]*bht.back_h_simple(var["T_abs_mean"],stepConditions["T_back"],hyp["theta"],L_c)
            else:
                print('here')
                var["h_back"] = 1/(1/(hyp["coeff_h_back"]*bht.back_h_fins(var["T_abs_mean"],stepConditions["T_back"],hyp["theta"],L_c,D,componentSpecs["Heta"]))+0.01)

        # Cas without fins
        else:
            # theta est l'inclinaison du panneau componentSpecs rapport à l'horizontale
            
            if componentSpecs["insulated"] == 1:
                T_ref = var["T_ins_mean"]
            else:
                T_ref = var["T_abs_mean"]

            if hyp['method_h_back_hx'] == "free_with_coeff": 
                var["h_back"] = hyp["coeff_h_back"]*bht.back_h_simple(T_ref,stepConditions["T_back"],hyp["theta"],L_c)
            elif hyp['method_h_back_hx'] == "free":
                var["h_back"] = bht.back_h_simple(T_ref,stepConditions["T_back"],hyp["theta"],L_c)
            elif hyp['method_h_back_hx'] == "mixed":
                var["h_back"] = bht.back_h_mixed(T_ref,stepConditions["T_back"],stepConditions["u_back"],hyp["theta"],L_c)

            # res = hyp["coeff_h_back"]*bht.back_h_simple(T_ref,stepConditions["T_back"],hyp["theta"],L_c)
            # print('res',res)

            # if h1 == None:
            #     print('res = None in h_back()')
            #     print('L_c',L_c)
            #     print('T_abs',var["T_abs_mean"])
            #     print('theta',hyp["theta"])
            #     print('T_back',stepConditions["T_back"])
            #     print('h_back_calculated',res)
            # else:

            #     T_back_changed = hyp["T_back_changed"]
            #     if T_back_changed > 0:
            #         T_back_changed += 273.15
            #         h2 = hyp["coeff_h_back"]*bht.back_h_mixed(T_ref,T_back_changed,stepConditions["u_back"],hyp["theta"],L_c)
            #         h1 = (h2*T_back_changed)/stepConditions["T_back"]
            #     else:
            #         pass                          

            #     var["h_back"] = h1

def h_top_mean(componentSpecs,stepConditions,var,hyp):
    
    old_h_top = var["h_top_g"]
    h_top(componentSpecs,stepConditions,var,hyp)

    new_h_top = var["h_top_g"]

    var["h_top_g"] = (old_h_top+new_h_top)/2

def h_back_mean(componentSpecs,stepConditions,var,hyp):

    old_h_back = var["h_back"]
    h_back(componentSpecs,stepConditions,var,hyp)

    new_h_back = var["h_back"]

    var["h_back"] = (old_h_back+new_h_back)/2



def h_back_tube(componentSpecs,stepConditions,var,hyp):
    """Calculates the back heat transfer coefficient for the tube and stores it in var["h_back_tube"]
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None"""

    # Cas où le tuyau est bien lié à l'absorbeur (conductif)
    if componentSpecs["l_c"] > 0. :
        var["h_back_tube"] = var["h_back"]

    # Cas où le tuyau est décollé de l'échangeur (radiatif seul)
    # Alors on calcule dans tous les cas comme si c'était unn tube cylindrique
    else:

        L_c = componentSpecs['H_tube']

        if componentSpecs["insulated"] == 1 and stepConditions["compt"] >= 1:

            T_ref = var["T_ins_tube_mean"]

        else:
            T_ref = var["T_tube_mean"]

        # if componentSpecs["tube_geometry"]=="rectangular" or componentSpecs["tube_geometry"]=="square":
        #     res = bht.back_h_mixed(T_ref,stepConditions["T_back"],stepConditions["u_back"],hyp["theta"],L_c)
        # else:
        res = bht.back_h_cylinder(T_ref,stepConditions["T_back"],L_c)

        if componentSpecs["is_inlet_man"] == 1 and hyp["inlet_manifold_in_wind"] == 1:
            T_film = (var["T_tube_mean"]+stepConditions["T_amb"])/2
            D = componentSpecs["D_tube"]+componentSpecs["lambd_riser_back"]
            nu = bht.air_nu(T_film)
            Re = (stepConditions["u"]*D)/nu
            Pr = bht.air_Pr()
            k = bht.air_k(T_film)
            Nu_forced = ht.conv_external.Nu_external_cylinder(Re,Pr)
            h_forced = Nu_forced*(k/D)
        else:
            h_forced = 0.
                
        var["h_back_tube"] = (res**3 + h_forced**3)**(1/3)

def h_back_fins(componentSpecs,stepConditions,var,hyp):
    """Calculates the back heat transfer coefficient for the fins and stores it in var["h_back_fins"]
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None"""



    if hyp["h_back_fins_calc"] == "tube":
        var["h_back_fins"] = var["h_back_tube"]
    elif hyp["h_back_fins_calc"] == "abs":
        var["h_back_fins"] = var["h_back"]
    elif hyp["h_back_fins_calc"] == "TS":
        L_c = componentSpecs["L_fin"]
        D = componentSpecs["D"]
        h_free = hyp["coeff_h_back_fins_free"]*bht.back_h_fins(var["T_tube_mean"],stepConditions["T_back"],hyp["theta"],L_c,D,componentSpecs["Heta"])
        # h_free = 0.
        h_forced = hyp["coeff_h_back_fins_forced"]*bht.ht_fins_forced_wiki(componentSpecs["L_fin"],componentSpecs["D"],stepConditions["u_back"]+0.9*stepConditions["u"])
        var["h_back_fins"] = (h_free**3 + h_forced**3)**(1/3) + hyp["offset_h_back_fins"]
    else:
        pass

def h_rad_back_tube(componentSpecs,stepConditions,var,hyp):
    """Calculates the radiation heat transfer coefficient between the tube and the ambient air and stores it in var["h_rad_back_tube"]
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None
    """

    sigma = hyp["sigma"]
    T_back_rad = stepConditions["T_back"] # hypothèse T_amb = T_back   
    if componentSpecs["insulated"] == 1 and stepConditions["compt"] >= 1:
        T_ref = var["T_ins_tube_mean"]
        eps = componentSpecs["eps_ins"]
    else: 
        T_ref = var["T_tube_mean"]
        eps = componentSpecs["eps_hx_back"]

    h = eps*sigma*(T_ref+T_back_rad)*(T_ref**2+T_back_rad**2)
    var["h_rad_back_tube"]=h


def h_rad_back(componentSpecs,stepConditions,var,hyp):
    """Calculates the radiation heat transfer coefficient between the absorber and the ambient air and stores it in var["h_rad_back"]
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None
    """

    sigma = hyp["sigma"]
    T_back_rad = stepConditions["T_back"] # hypothèse T_amb = T_back   
    if componentSpecs["insulated"] == 1:
        T_ref = var["T_ins_mean"]
        eps = componentSpecs["eps_ins"]
    else: 
        T_ref = var["T_abs_mean"]
        eps = componentSpecs["eps_hx_back"]
   
    T_back_rad_changed = hyp["T_back_rad_changed"]

    if T_back_rad_changed > 0:
        T_back_rad_changed += 273.15
        sigma = hyp["sigma"]

        h2 = eps*sigma*(T_ref+T_back_rad_changed)*(T_ref**2+T_back_rad_changed**2)
        h1 = (h2*T_back_rad_changed)/T_back_rad
    else:
        h1 = eps*sigma*(T_ref+T_back_rad)*(T_ref**2+T_back_rad**2)

    var["h_rad_back"]=h1

def a0(componentSpecs,stepConditions,var):
    """Calculates the a0 factor and stores it in var["a0"]
    
    $$a_0 = \frac{1}{h_{top}+h_{rad}+1/R_g}\left(\alpha_gG+h_{top}T_{amb}+h_{rad}T_{sky}\right)$$

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
    
    $$a_3 = \left(h_{top}+h_{rad}\right)a_1$$
    
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
    var["a2"] = - var["a0"]/componentSpecs["R_g"]

def S_star(componentSpecs,stepConditions,var):
    var["S_star"] = var["S"] - var["a2"] + var["h_rad"]*stepConditions["T_sky"]

def h_rad_f(componentSpecs,stepConditions,var,hyp):
    """Calculates the radiation heat transfer coefficient between the tube and the absorber and stores it in var["h_rad_f"]
    
    $$h_{rad,f} = \epsilon_{hx} \sigma (T_{tube}+T_{abs})(T_{tube}^2+T_{abs}^2)$$
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None
    """
    # var["h_rad_f"] = hyp["h_rad_f"]

    if componentSpecs["l_c"] > 0. :
        h=0.

    else:

        sigma = hyp["sigma"]

        T_ref = var["T_tube_mean"]
        T_B = var["T_Base_mean"]

        eps = componentSpecs["eps_hx_top"]

        h = eps*sigma*(T_ref+T_B)*(T_ref**2+T_B**2)

    var["h_rad_f"] = hyp["coeff_h_rad_f"]*h

# Pour les ailettes, on ne prend pas en compte la composante radiative h_rad_back

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

def h_rad(componentSpecs,stepConditions,var,hyp):
    """Calculates the radiation heat transfer coefficient and stores it in var["h_rad"]
    
    $$
    h_{rad} = \epsilon \sigma (T_{PV}+T_{sky})(T_{PV}^2+T_{sky}^2)
    $$

    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None"""
    
    eps = componentSpecs["eps_PV"]
    sigma = hyp["sigma"]
    T_sky = stepConditions["T_sky"]

    T_PV = var["T_PV"]

    tau_g_IR = componentSpecs["tau_g_IR"]

    h = tau_g_IR*eps*sigma*(T_PV+T_sky)*(T_PV**2+T_sky**2)
    var["h_rad"]=h
    #var["h_rad"]=0.00001

def h_rad_g(componentSpecs,stepConditions,var,hyp):
    """Calculates the radiative heat transfer coefficient between the glass and the sky and stores it in var["h_rad_g"]
    
    Args:
        componentSpecs (dict): dictionary containing the parameters of the PVT panel
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None"""
    
    var["h_rad_g"] = bht.h_rad(componentSpecs["eps_g"],var["T_glass"],stepConditions["T_sky"])

def h_rad_tube_sky(componentSpecs,stepConditions,var,hyp):
    """Calculates the radiative heat transfer coefficient between the tube and the sky and stores it in var["h_rad_tube_sky"]
    
    Args:
        componentSpecs (dict): dictionary containing the parameters of the PVT panel
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None
    """
    eps = componentSpecs["eps_hx_top"]
    sigma = hyp["sigma"]
    T_sky = stepConditions["T_sky"]

    T_tube = var["T_tube_mean"]

    tau_g_IR = componentSpecs["tau_g_IR"]

    h = tau_g_IR * eps*sigma*(T_tube+T_sky)*(T_tube**2+T_sky**2)

    var["h_rad_tube_sky"]=h

# Depends on T_PV
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

    #T_PV = var["T_PV"]
    X_celltemp = var["X_celltemp"]

    eta = eta_nom*X_celltemp*X_rad*X_corr
    var["eta_PV"] = eta
    #var["eta_PV"] = 0.15

# net absorbed solar radiation (total absorbed - PV power production)
def S(componentSpecs,stepConditions,var):
    """Calculates the net absorbed solar radiation and stores it in var["S"]
    
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

def Fp(componentSpecs, var):
    """Calculates the Fp factor and stores it in var["Fp"]
    
    $$
    F_p = \frac{1}{h_{rad}R_{inter}+R_{inter}/R_t+1}
    $$
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters    
        var (dict): dictionary containing the variables
        
    Returns:
        None"""


    R_inter = componentSpecs["R_inter"]

    #T_PV = var["T_PV"]
    h_rad = var["h_rad"]
    a3 = var["a3"]

    Fp = 1/(a3*R_inter+h_rad*R_inter+1)
    var["Fp"] = Fp

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
    h_rad_f = var["h_rad_f"]

    Fprime = var["Fp"]
    
    # j = 1/(R_inter*Fprime)+1/(R_b*Fprime)-1/R_inter+h_rad_f/Fprime

    # j = 1/(R_inter*Fprime)+1/(R_b*Fprime)-1/R_inter

    j = 1/(Fprime*R_b) + 1/(R_inter*Fprime) - 1/R_inter

    if componentSpecs["fin_2"] == 1:

        gamma_int = var["gamma_2_int"]

        j += (gamma_int)/Fprime

    var["j"] = j

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

    h_rad_f = var["h_rad_f"]

    # if h_rad_f != 0:
    #     T_tube_mean = var["T_tube_mean"]
    #     b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime) + (h_rad_f*T_tube_mean)/Fprime
    # else:
    #     b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime)

    # b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime)

    S_star = var["S_star"]

    b = S_star + T_back/(R_b*Fprime)


    if componentSpecs["fin_2"]==1:
        gamma_int = var["gamma_2_int"]

        b += (gamma_int*T_back)/Fprime

    var["b"] = b

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
    
    q = k_abs*lambd_abs*m*((b/j)-T_B)*tanh_or_inverse(m*L_af)
    var["qp_fin"] = q

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
 

    h_rad_f = var["h_rad_f"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_tube_sky = var["h_rad_tube_sky"]
    p_tube_sky = componentSpecs["p_tube_sky"]

    # vérifier l'homogénéité
    var["e0"] = 1/chi + h_rad_f*p_ext_tube_rad + C_B + gamma + h_rad_tube_sky*p_tube_sky

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
    h_rad_f = var["h_rad_f"]

    e0 = var["e0"]

    var["e1"] = (1/e0)*(C_B+p_ext_tube_rad*h_rad_f)

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

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

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

# non utilisée
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

    h_rad_f = var["h_rad_f"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]

    h_rad_tube_sky = var["h_rad_tube_sky"]
    p_tube_sky = componentSpecs["p_tube_sky"]

    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    var["f0"] = C_B + h_rad_f*p_ext_tube_rad + gamma + h_rad_tube_sky*p_tube_sky

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
    h_rad_f = var["h_rad_f"]
    p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]

    e1 = var["e1"]
    f0 = var["f0"]
    # print('b1')
    # print('e1',e1)
    # print('gamma',gamma)

    var["b1"] = 1/(C_B + h_rad_f*p_ext_tube_rad - f0*e1)
                   
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
    h_rad_f = var["h_rad_f"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    e3 = var["e3"]
    f0 = var["f0"]
    # c'est vérifié c'est égal
    # var["b3"] = (-1/(C_B + h_rad_f*p_ext_tube_rad))*gamma
    # var["b3"] =  (-1/(C_B + h_rad_f*p_ext_tube_rad))*gamma

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
    
    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int
    h_rad_f = var["h_rad_f"]
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

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_f = var["h_rad_f"]

    e1 = var["e1"]

    # var["d1"] = (-gamma*d0)/R_B + h_rad_f*(1-(d0/R_B))
    var["d1"] = -e1*(gamma+h_rad_f*p_ext_tube_rad) + h_rad_f*p_ext_tube_rad

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

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_f = var["h_rad_f"]

    # var["d2"] = (-gamma*d0)/chi - (h_rad_f*d0)/chi

    e2 = var["e2"]

    var["d2"] = -e2*(gamma+h_rad_f*p_ext_tube_rad)

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

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_f = var["h_rad_f"]

    # var["d3"] = -gamma**2*d0 + gamma - (h_rad_f*d0)/(1/gamma)

    e3 = var["e3"]
    var["d3"] = -e3*(gamma+h_rad_f*p_ext_tube_rad) + gamma

def d4(componentSpecs, var):
    """Calculates the d4 factor and stores it in var["d4"]
    
    $$d_4 = \frac{1}{c_0\chi}\frac{1}{1-c_2}$$
    
    Args:
        componentSpecs (dict): dictionary containing the parameters
        var (dict): dictionary containing the variables
        
    Returns:
        None"""

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_f = var["h_rad_f"]

    # var["d3"] = -gamma**2*d0 + gamma - (h_rad_f*d0)/(1/gamma)

    e4 = var["e4"]

    var["d4"] = -e4*(gamma + h_rad_f*p_ext_tube_rad)


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

    # K = b2 * ( -D_tube*Fprime*((l_c/D_tube)*(h_rad+1/R_t)+(iota/D_tube)/(R_b*Fprime))-2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af) )
    # T = b1 * ( D_tube*Fprime*((l_c/D_tube)*(h_rad+1/R_t)+(iota/D_tube)/(R_b*Fprime))+2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af))
    # E = D_tube*Fprime*((l_c/D_tube)*(S+h_rad*T_sky+T_amb/R_t)+((iota/D_tube)*(1-b3)*T_back)/(R_b*Fprime))+2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*((b/j) - b3*T_back) + (l_c/R_t)*(Fprime-1)*b3*T_back

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
    i5 = 2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*(b/j)
    i4 = -2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*b4
    i3 = -2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*b3
    i2 = -2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*b2
    i1 = -2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*b1

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
    h_rad_f = var["h_rad_f"]

    R_B = 1/(C_B+p_ext_tube_rad*h_rad_f)
    p_int_tube = componentSpecs["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)
    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"];p_ext_tube_rad = componentSpecs["p_ext_tube_rad"]
    R_2 = componentSpecs["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

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

    # if componentSpecs["fin_0"] == 1:
    #     gamma_0_int = var["gamma_0_int"]
    # else:
    #     gamma_0_int = 0

    # if componentSpecs["fin_1"] == 1:
    #     gamma_1_int = var["gamma_1_int"]
    #     # print(gamma_1_int)
    # else:
    #     gamma_1_int = 0
        

    # k = componentSpecs["k_ail"]
    # C_B_f = (p_ext_tube*componentSpecs["k_riser"])/componentSpecs["lambd_riser_back"]
    # h_fluid = var["h_fluid"]

    # chi = 1/(h_fluid*p_int_tube)+1/C_B_f

    # gamma_back = p_ext_tube/(R_2+1/h_back)

    # zeta = (gamma_back + gamma_1_int + gamma_0_int)/(1+chi*(gamma_back+gamma_1_int+gamma_0_int))
    # # print(zeta)

    # a += (-N_harp/(mdot*Cp))*zeta
    # b += (N_harp/(mdot*Cp))*zeta*T_back

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

    var["T_absfin_mean"] = b/j + (T_B_mean-(b/j))*tanh_or_inverse(m*L_af)/(m*L_af)

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
    h_rad_f = var["h_rad_f"]

    var["Qdot_absfin_tube"] = L*p_ext_tube_rad*h_rad_f*(T_absfin_m-T_tube_m)

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


def T_ins_tube_mean(componentSpecs,var):
    R_2 = componentSpecs["R_2"]
    T_tube_mean = var["T_tube_mean"]
    Qdot_tube_back = var["Qdot_tube_back"]

    L = componentSpecs["L_tube"]
    p_ext_tube = componentSpecs["p_ext_tube"]

    var["T_ins_tube_mean"] = T_tube_mean - R_2*Qdot_tube_back/(L*p_ext_tube)

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

def T_ins_absfin_mean(componentSpecs,var):
    R_2 = componentSpecs["R_2"]
    T_absfin_mean = var["T_absfin_mean"]
    Qdot_absfin_back = var["Qdot_absfin_back"]

    L = componentSpecs["L_tube"]
    L_af = componentSpecs["L_af"]

    var["T_ins_absfin_mean"] = T_absfin_mean - R_2*Qdot_absfin_back/(L*2*L_af)

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


# def T_ins_mean(componentSpecs,stepConditions,var):
#     T_ref = var["T_abs_mean"]
#     R_2 = componentSpecs["R_2"]
#     h_back = var["h_back"] + var["h_rad_back"]
#     T_back = stepConditions["T_back"]

#     var["T_ins_mean"] = T_ref + (R_2/(R_2+1/h_back)) * (T_back - T_ref)


def T_ins_mean(componentSpecs,var):

    T_ins_tube_mean = var["T_ins_tube_mean"]
    T_ins_absfin_mean = var["T_ins_absfin_mean"]
    l_B = componentSpecs["l_B"]
    lambd_riser_back = componentSpecs["lambd_riser_back"]
    L_af = componentSpecs["L_af"]
    W = componentSpecs["W"]

    var["T_ins_mean"] = (l_B*T_ins_tube_mean + 2*L_af*T_ins_absfin_mean)/W



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

def Qdot_fluid(componentSpecs,stepConditions,var):

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
    h_rad_f = var["h_rad_f"]


    var["Qdot_Base_tube"] = L*(T_Base_m-T_tube_m)*(C_B+h_rad_f*p_ext_tube_rad)

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
        \item T_glass_mean()
        \item qp_PV_Base()
        \item qp_Base_back()
        \item qp_fin()
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
    
    h_rad_g(componentSpecs,stepConditions,var,hyp)
    h_rad(componentSpecs,stepConditions,var,hyp)

    if componentSpecs["fin_0"] == 1 or componentSpecs["fin_1"] == 1 or componentSpecs["fin_2"] == 1:
        h_top_mean(componentSpecs,stepConditions,var,hyp)
    else:
        h_top(componentSpecs,stepConditions,var,hyp)

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
    T_fluid_out(componentSpecs,var)
    q_tube_fluid(componentSpecs,stepConditions,var)
    T_fluid_mean(componentSpecs,var)
    T_Base_mean(componentSpecs,stepConditions,var)
    T_tube_mean(componentSpecs,stepConditions,var)
    T_absfin_mean(componentSpecs,stepConditions,var)
    T_abs_mean(componentSpecs,var)


    Qdot_tube_back(componentSpecs,stepConditions,var)
    Qdot_absfin_back(componentSpecs,stepConditions,var)   
    T_ins_tube_mean(componentSpecs,var)
    T_ins_absfin_mean(componentSpecs,var)
    T_ins_mean(componentSpecs,var)

    Qdot_ins_conv(componentSpecs,stepConditions,var)
    Qdot_ins_rad(componentSpecs,stepConditions,var)

    if hyp["calc_h_back_mean"]==1:
        h_back_mean(componentSpecs,stepConditions,var,hyp)
    else:
        h_back(componentSpecs,stepConditions,var,hyp)

    h_rad_back(componentSpecs,stepConditions,var,hyp)
    h_back_tube(componentSpecs,stepConditions,var,hyp)
    h_rad_tube_sky(componentSpecs,stepConditions,var,hyp)
    h_rad_back_tube(componentSpecs,stepConditions,var,hyp)
    h_back_fins(componentSpecs,stepConditions,var,hyp)

    h_rad_f(componentSpecs,stepConditions,var,hyp)

    T_PV_mean(componentSpecs,stepConditions,var)
    T_PV_Base_mean(componentSpecs,stepConditions,var)
    T_PV_absfin_mean(componentSpecs,var)
    T_glass_mean(componentSpecs,stepConditions,var)

    qp_PV_Base(componentSpecs,var)
    qp_Base_back(componentSpecs,stepConditions,var)
    qp_fin(componentSpecs,var)

    Cp(componentSpecs,stepConditions,var,hyp)


def compute_power(componentSpecs,stepConditions,var):
    Qdot_top_conv(componentSpecs,stepConditions,var)
    Qdot_top_rad(componentSpecs,stepConditions,var)
    Qdot_sun_glass(componentSpecs,stepConditions,var)
    Qdot_sun_PV(componentSpecs,stepConditions,var)
    Qdot_glass_PV(componentSpecs,stepConditions,var)
    Qdot_PV_sky(componentSpecs,stepConditions,var)
    Qdot_PV_plate(componentSpecs,var)
    # Qdot_abs_back1(componentSpecs,stepConditions,var)
    Qdot_PV_Base(componentSpecs,var)
    Qdot_PV_absfin(componentSpecs,var)
    Qdot_Base_back(componentSpecs,stepConditions,var)
    Qdot_fluid(componentSpecs,stepConditions,var)
    Qdot_tube_back(componentSpecs,stepConditions,var)
    Qdot_tube_sky(componentSpecs,stepConditions,var)
    Qdot_Base_tube(componentSpecs,var)
    # qp_fluid_back(componentSpecs,var)
    qp_fin(componentSpecs,var)
    Qdot_absfin_Base(componentSpecs,var)
    Qdot_abs_back2(componentSpecs,var)
    Qdot_fluid_back(componentSpecs,stepConditions,var)
    Qdot_tube_back(componentSpecs,stepConditions,var)
    Qdot_absfin_back(componentSpecs,stepConditions,var)
    Qdot_absfin_back_conv(componentSpecs,stepConditions,var)
    Qdot_absfin_back_rad(componentSpecs,stepConditions,var)
    # Qdot_absfin_tube(componentSpecs,stepConditions,var)
    Qdot_tube_back_wo_ins_conv(componentSpecs,stepConditions,var)
    Qdot_tube_back_wo_ins_rad(componentSpecs,stepConditions,var)
    Qdot_ins_tube_back_conv(componentSpecs,stepConditions,var)
    Qdot_ins_tube_back_rad(componentSpecs,stepConditions,var)
    Qdot_ins_absfin_back_conv(componentSpecs,stepConditions,var)
    Qdot_ins_absfin_back_rad(componentSpecs,stepConditions,var)

    if componentSpecs["fin_0"]==1 or componentSpecs["fin_1"]==1:
        Qdot_f01(componentSpecs,stepConditions,var)
    else:
        var["Qdot_f01"] = 0.

    power_balance_1(componentSpecs,var)
    power_balance_3(componentSpecs,var)

    if componentSpecs["fin_0"] == 1:
        qp_f0(componentSpecs,stepConditions,var)
    if componentSpecs["fin_1"] == 1:
        qp_f1(componentSpecs,stepConditions,var)
    if componentSpecs["fin_2"] == 1:
        pass

    T_B_check(componentSpecs,stepConditions,var)

def simu_one_steady_state_all_he(componentSpecs,stepConditions,hyp):
    
    res = {}

    save_T_fluid_in0 = stepConditions["T_fluid_in0"]

    # Test the main part

    slices_df, df_one, its_data_list = simu_one_steady_state(componentSpecs['main'],stepConditions,hyp)
    res['main'] = {'slices_df':slices_df.copy(),'df_one':df_one.copy(),'its_data_list':its_data_list.copy()}

    hyp['h_back_prev'] = df_one['h_back'].values[0] # h_back de l'absorbeur
    hyp['h_top_man'] = df_one["h_top_g"].values[0]

    if len(componentSpecs['decomp'])>1:
        decomp = 1

        for part,part_name in list(componentSpecs['decomp'].items())[1:]:

            slices_df, df_one, its_data_list = simu_one_steady_state(componentSpecs[part],stepConditions,hyp)
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
                last_past = list(componentSpecs['decomp'].keys())[-1]
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
                    Aire = componentSpecs[part]["N_harp"]*componentSpecs[part]["W"]*componentSpecs[part]["L_tube"]
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
    return df_one,res

def initialize_var(var,componentSpecs,stepConditions,hyp,i):
    # Initialize the var dictionary with all necessary keys and values
    
    var = {'Slice' : i,
           'T_PV0':0,
           'Cp': hyp['Cp0'],
            'h_rad_f':hyp['h_rad_f0'],
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
    h_fluid(componentSpecs,stepConditions,var,hyp)

    return var

def update_heat_transfer_coefficients(var, hyp, index):
    # This function should update heat transfer coefficients for the current index
    # based on prior values or initial guesses
    # Replace '...' with actual logic to calculate or update the coefficients
    if index == 0:
        ...
    else:
        ...

mean_list = ["T_glass","T_PV","T_PV_Base_mean","T_PV_absfin_mean","T_abs_mean","T_Base_mean","T_absfin_mean","T_ins_mean","T_ins_tube_mean","T_ins_absfin_mean","T_tube_mean","T_fluid_mean","h_top_g","h_rad","h_back","h_rad_back","h_back_tube","h_rad_back_tube","h_back_fins","h_rad_f","h_fluid","X_celltemp","eta_PV","S"]
add_list = ["Qdot_sun_glass","Qdot_sun_PV","Qdot_top_conv","Qdot_top_rad","Qdot_glass_PV","Qdot_PV_sky","Qdot_PV_plate","Qdot_PV_Base","Qdot_PV_absfin","Qdot_absfin_Base","Qdot_absfin_back","Qdot_absfin_back_conv","Qdot_absfin_back_rad","Qdot_Base_tube","Qdot_Base_back","Qdot_tube_sky","Qdot_tube_fluid","Qdot_tube_back","Qdot_ins_tube_back_conv","Qdot_ins_tube_back_rad","Qdot_ins_absfin_back_conv","Qdot_ins_absfin_back_rad","Qdot_tube_back_conv","Qdot_tube_back_rad","Qdot_absfin_back","Qdot_f01"]


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

def pos(x):
    if x>=0:
        return x
    else:
        return 0

def neg(x):
    if x<=0:
        return x
    else:
        return 0

def pos_df(df,col_name):
    l = []
    for x in df[col_name]:
        l.append(pos(x))
    return l

def neg_df(df,col_name):
    l = []
    for x in df[col_name]:
        l.append(neg(x))
    return l   

def write_stepConditions_from_steadyStateConditions(steadyStateConditions_df,i,hyp):
    stepConditions = {'G':steadyStateConditions_df["G"][i],"T_amb":steadyStateConditions_df["T_amb"][i],"T_back":steadyStateConditions_df["T_amb"][i],"u":steadyStateConditions_df["u"][i], "u_back" : steadyStateConditions_df["u_back"][i], "T_fluid_in0":steadyStateConditions_df["T_fluid_in"][i]}
    change_T_sky(stepConditions,hyp,'TUV')  # calculate Gp and T_sky

    stepConditions['T_back'] = stepConditions['T_amb']
    stepConditions['T_back_rad'] = stepConditions['T_amb']

    stepConditions["mdot"] = steadyStateConditions_df["mdot"][i]

    stepConditions["guess_T_PV"] = (stepConditions["T_amb"]+stepConditions["T_fluid_in0"])/2

    return stepConditions

def simu_steadyStateConditions(componentSpecs,hyp,steadyStateConditions_df):
    
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

        df_one,res = simu_one_steady_state_all_he(componentSpecs,stepConditions,hyp)

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
    df_res['Qdot / AG'] = df_res['Qdot']/(componentSpecs['AG'])

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

def simu_condi_mpe(componentSpecs,steadyStateConditions_df,l,h_back,L,hyp):
    
    variables = ['N_test','T_guess','G', 'Gp', 'T_amb', 'u', 'T_abs','T_fluid_in', 'T_fluid_out']
    
    # Dataframe object
    df = pd.DataFrame(columns = variables)

    sigma = hyp["sigma"]

    compt_test = 0

    for i in range(1,len(steadyStateConditions_df)+1):

        stepConditions = {}

        stepConditions["G"]=steadyStateConditions_df["G"][i]

        # T_amb = T_back
        stepConditions["T_amb"]=steadyStateConditions_df["ta"][i]+273.15

        change_T_sky(componentSpecs,'TUV')

        # Back temperature = ambiant temperature
        stepConditions["T_back"]=stepConditions["T_amb"]

        # Change wind_speed in componentSpecs and adapt R_t
        change_u(componentSpecs,stepConditions,steadyStateConditions_df["U"][i])

        stepConditions["mdot"] = steadyStateConditions_df["mdot"][i]

        T_f_in_list = [steadyStateConditions_df["tin"][i]+273.15]                

        T_f_out = stepConditions["T_back"] + (T_f_in_list[0]-stepConditions["T_back"])*math.exp(-(l*h_back*L)/((stepConditions["mdot"]/componentSpecs["N_harp"])*componentSpecs["Cp"]))
        
        # len(T_f_out_list) = 1

        to_add = {'N_test' : compt_test, 'T_guess' : 293.15, 'G' : stepConditions["G"], 'Gp' : stepConditions["Gp"], 'T_amb' : stepConditions["T_amb"], 'h_back' : h_back, 'u' : stepConditions["u"], 'T_fluid_in' : T_f_in_list[0], 'T_abs' : 293.15,'T_fluid_out' : T_f_out}

        df_to_add = pd.DataFrame.from_dict({'row' : to_add.values()},orient='index',columns=to_add.keys())

        df = pd.concat([df,df_to_add])

        compt_test+=1

    # Analysing df

    # Be careful here you have zeros for some columns

    df['DT'] = df['T_fluid_out'] - df['T_fluid_in']
    df['Tm'] = (df['T_fluid_out'] + df['T_fluid_in'])/2
    df['T_m*'] = (df['Tm'] - df['T_amb'])/df['G']
    df['G x (T_m*)^2'] = df['G'] * df['T_m*']**2 * 0
    df['up x T_m*'] = (df['u'] - 3) * df['T_m*']
    df['Gp/G'] = df['Gp']/df['G']
    df['up'] = df['u'] - 3
    df['up x Gp/G'] = (df['up'] * df['Gp'])/df['G']
    df['G^3 x (T_m*)^4'] = df['G']**3 * df['T_m*']**4 * 0

    df['T_m en °C'] = df['Tm']-273.15

    coeff_density = [999.85,0.05332,-0.007564,0.00004323,-1.673e-7,2.447e-10]
    coeff_density = list(reversed(coeff_density))

    coeff_Cp = [4.2184,-0.0028218,0.000073478,-9.4712e-7,7.2869e-9,-2.8098e-11,4.4008e-14]
    coeff_Cp = list(reversed(coeff_Cp))

    df['density(T)'] = np.polyval(coeff_density,df['T_m en °C'])
    df['Cp(T)'] = np.polyval(coeff_Cp,df['T_m en °C'])*1000

    df['mdot'] = df['density(T)']*(stepConditions["mdot"]/1000)

    df['Qdot'] = df['mdot']*df['Cp(T)']*df['DT']
    df['Qdot / (AG x G)'] = df['Qdot']/(componentSpecs['AG']*df['G'])

    ones = pd.DataFrame(np.ones(len(df['T_m*'])),columns=['Ones'])
    ones_column = ones["Ones"]
    df_mat = df[['T_m*','G x (T_m*)^2','up x T_m*','Gp/G','up','up x Gp/G','G^3 x (T_m*)^4']].join(ones_column)

    matrice = df_mat.to_numpy()
    B = df['Qdot / (AG x G)'].to_numpy()

    X=np.linalg.lstsq(matrice, B, rcond = -1)

    #_ = plt.plot(df['T_m*'].to_numpy(), B, 'o', label='Original data', markersize=2)
    #_ = plt.plot(df['T_m*'].to_numpy(), np.dot(matrice,X[0]), 'o', label='Fitted line',markersize=2)
    #_ = plt.legend()
    #plt.show()

    return df,X

def simu_condi_mpe_big(componentSpecs,stepConditions,steadyStateConditions_df,l,L,h_back_top,h_back_bottom,N_harp,hyp):
    
    variables = ['N_test', 'mdot', 'T_guess','G', 'Gp', 'T_amb', 'u', 'T_abs','T_fluid_in', 'T_fluid_out']
    
    # Dataframe object
    df_res = pd.DataFrame(columns = variables)

    sigma = hyp["sigma"]

    compt_test = 0

    for i in range(1,len(steadyStateConditions_df)+1):

        stepConditions["G"]=steadyStateConditions_df["G"][i]

        # T_amb = T_back
        stepConditions["T_amb"]=steadyStateConditions_df["ta"][i]+273.15

        change_T_sky(componentSpecs,'TUV')

        # Back temperature = ambiant temperature
        stepConditions["T_back"]=stepConditions["T_amb"]

        # Change wind_speed in componentSpecs and adapt R_t
        change_u(componentSpecs,stepConditions,steadyStateConditions_df["U"][i])

        stepConditions["mdot"] = steadyStateConditions_df["mdot"][i]

        T_f_in_list = [steadyStateConditions_df["tin"][i]+273.15]                

        T_f_out = stepConditions["T_back"] + (T_f_in_list[0]-stepConditions["T_back"])*math.exp(-(l*L*h_back_top+l*L*h_back_bottom)/((stepConditions["mdot"]/N_harp)*componentSpecs["Cp"]))
        
        # len(T_f_out_list) = 1

        to_add = {'N_test' : compt_test, 'mdot' : stepConditions["mdot"], 'T_guess' : 293.15, 'G' : stepConditions["G"], 'Gp' : stepConditions["Gp"], 'T_amb' : stepConditions["T_amb"], 'h_back' : h_back_top, 'u' : stepConditions["u"], 'T_fluid_in' : T_f_in_list[0], 'T_abs' : 293.15,'T_fluid_out' : T_f_out}

        df_to_add = pd.DataFrame.from_dict({'row' : to_add.values()},orient='index',columns=to_add.keys())

        df_res = pd.concat([df_res,df_to_add])

        compt_test+=1

    # Analysing df

    # Be careful here you have zeros for some columns

    tab = pd.DataFrame()

    df_res['DT'] = df_res['T_fluid_out'] - df_res['T_fluid_in']
    df_res['Tm'] = (df_res['T_fluid_out'] + df_res['T_fluid_in'])/2
    df_res['T_m en °C'] = df_res['Tm']-273.15

    # tab['G'] = 0. * df_res['G'] # a0
    tab['-(T_m - T_a)'] = -(df_res['Tm'] - df_res['T_amb']) # a1
    # tab['-(T_m - T_a)^2'] = -(df_res['Tm'] - df_res['T_amb'])**2 # a2
    # tab['-(T_m - T_a)^2'] = 0.*df_res['Tm'] # a2
    # tab['-up x (T_m - T_a)'] = (df_res['u'] - 3) * tab['-(T_m - T_a)'] # a3
    # tab['Gp'] = df_res['Gp'] # a4
    # tab['Gp'] = 0. * df_res['Gp'] # a4
    # tab['up x G'] = -(df_res['u'] - 3) * df_res['G'] # a6
    # tab['up x Gp'] = -(df_res['u'] - 3) * df_res["Gp"] # a7
    # tab['-(T_m - T_a)^4'] = 0. * (-tab['-(T_m - T_a)']**4) # a8

    # coeff_density = [999.85,0.05332,-0.007564,0.00004323,-1.673e-7,2.447e-10]
    # coeff_density = list(reversed(coeff_density))

    coeff_Cp = [4.2184,-0.0028218,0.000073478,-9.4712e-7,7.2869e-9,-2.8098e-11,4.4008e-14]
    coeff_Cp = list(reversed(coeff_Cp))

    # df_res['density(T)'] = np.polyval(coeff_density,df_res['T_m en °C'])
    df_res['Cp(T)'] = np.polyval(coeff_Cp,df_res['T_m en °C'])*1000

    # df_res['mdot'] = df_res['density(T)']*(stepConditions["mdot"]/1000)

    df_res['Qdot'] = df_res['mdot']*df_res['Cp(T)']*df_res['DT']
    df_res['Qdot / AG'] = df_res['Qdot']/(componentSpecs['AG'])

    tab = tab.astype('float64')

    matrice = tab.to_numpy()
    B = df_res['Qdot / AG'].to_numpy()

    X = np.linalg.lstsq(matrice, B, rcond = -1)

    #_ = plt.plot(df['T_m*'].to_numpy(), B, 'o', label='Original data', markersize=2)
    #_ = plt.plot(df['T_m*'].to_numpy(), np.dot(matrice,X[0]), 'o', label='Fitted line',markersize=2)
    #_ = plt.legend()
    #plt.show()

    df_res_to_concat = df_res.drop(columns=["G","Gp"])

    df_res = pd.concat([tab,df_res_to_concat],axis=1)

    return df_res,X


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


