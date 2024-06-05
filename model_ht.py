import math
import copy
import pandas as pd
import numpy as np
import heat_transfer as bht
import ht

from CoolProp.CoolProp import PropsSI

# INTERNAL CONVECTIVE

def h_fluid(componentSpecs,stepConditions,var,hyp):
    """Calculates the convective heat transfer coefficient between the fluid and the tube wall and stores it in var["h_fluid"]

    - If the exchanger is a flow distribution like, the considered flow rate value is the total flow rate in the PVT divided by 2 (because the flow is distributed/collected to/form the channels all along the manifold)
    - If the Re < 2000, the Nusselt number is calculated with the Shan-London correlation for laminar flow in a rectangular duct or with the 'Laminar - constant Q' method for a circular duct
    - If the Re >= 2000, the Nusselt number is calculated with the Colburn correlation for turbulent flow
    
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

    if (componentSpecs["is_inlet_man"] == 1 or componentSpecs["is_outlet_man"] == 1) and componentSpecs["geometry"] == "harp":
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

def get_CFD_value(componentSpecs, stepConditions, var, hyp, h, phi, T_1, T_2):
    big_it = hyp['big_it']
    CFD_ht = pd.read_csv(hyp['CFD_ht_path']+f'_{big_it}.csv',sep=';')
    var[h] = abs( CFD_ht.loc[CFD_ht['component'] == componentSpecs['name']][phi].values[0] / (var[T_1] - stepConditions[T_2]) )

# EXTERNAL CONVECTIVE

def h_top_g(componentSpecs,stepConditions,var,hyp):
    """Calculates the convective heat transfer coefficient between the top of the panel (glass) and the ambient air and stores it in var["h_top_g"]
    
    - If the component is a manifold:
        - If the method is 'like_exchanger', the coefficient is the same as the one found when solving the 'main'
    
    - If the component is a heat exchanger part:

        - These coefficients are computed :

            - $h_{free}$ = top_h_simple(...)
            - $h_{forced_turbulent}$ = h_top_forced_turbulent(...)
            - $h_{forced}$ = h_top_forced(...)
            - $h_{custom}$ = h_top_custom(...)

        - If the method is 'free_with_coeff', the coefficient is the free convection coefficient multiplied by a coefficient
        - If the method is 'free', the coefficient is the free convection coefficient
        - If the method is 'forced_turbulent_with_coeff', the coefficient is the forced convection coefficient in the turbulent regime multiplied by a coefficient
        - If the method is 'forced_turbulent', the coefficient is the forced convection coefficient in the turbulent regime
        - If the method is 'forced_with_coeff', the coefficient is the forced convection coefficient multiplied by a coefficient
        - If the method is 'forced', the coefficient is the forced convection coefficient
        - If the method is 'free_forced_with_coeff', the coefficient is the free convection coefficient multiplied by a coefficient if the wind speed is lower than 0.1 m/s, otherwise it is the forced convection coefficient multiplied by a coefficient
        - If the method is 'mixed_with_coeff', the coefficient is the cubic mean of the free convection coefficient and the forced convection coefficient multiplied by a coefficient
        - If the method is 'mixed', the coefficient is the cubic mean of the free convection coefficient and the forced convection coefficient
        - If the method is 'custom', the coefficient is the custom convection coefficient

    Args:
        componentSpecs (dict): dictionary containing the parameters of the PVT panel
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None"""
    
    if hyp['method_h_top_g_exchanger'] == 'CFD':
        get_CFD_value(componentSpecs, stepConditions, var, hyp, 'h_top_g', 'phi_top', 'T_glass', 'T_amb')
        return

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

def h_top_mean(componentSpecs,stepConditions,var,hyp):
    """Calculates the mean h_top between the value at iteration n-1 and the value calculated at iteration n and stores it in var["h_top_g"]
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None"""
    
    old_h_top = var["h_top_g"]
    h_top_g(componentSpecs,stepConditions,var,hyp)

    new_h_top = var["h_top_g"]

    var["h_top_g"] = (old_h_top+new_h_top)/2

def h_back_abs(componentSpecs,stepConditions,var,hyp):
    """Calculates the convective heat transfer coefficient between the absorber and the ambient air and stores it var["h_back"]

    - If the component is an anomaly, the coefficient is the same as the one found when solving the 'main'

    - If the component is a manifold:

        - The characteristic length L_c is the height of the tube H_tube

        - If the method is 'free_cylinder', the coefficient is calculated with the free convection coefficient for a cylinder
        - If the method is 'like_exchanger', the coefficient is the same as the one found when solving the 'main'

    - If the component is a heat exchanger part:

        - The characteristic length L_c is the length of the absorber L_abs if the orientation is 'portrait', otherwise it is the width of the absorber w_abs

        - If the component has fins:

            - If the number of fins is less than 24, the coefficient is calculated with the free convection coefficient for a simple geometry
            - If the number of fins is greater than 24, the coefficient is calculated with the free convection coefficient for a finned geometry

        - If the component has no fins:

            - If the method is 'free_with_coeff', the coefficient is the free convection coefficient multiplied by a coefficient
            - If the method is 'free', the coefficient is the free convection coefficient
            - If the method is 'mixed', the coefficient is the mixed convection coefficient
    
    Args:
        componentSpecs (dict): dictionary containing the parameters of the PVT panel
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None"""
    
    if hyp['method_h_back_hx'] == 'CFD':
        get_CFD_value(componentSpecs, stepConditions, var, hyp, 'h_back', 'phi_abs', 'T_abs_mean', 'T_amb')
        return

    # Anomaly
    # If this part of the PVT is an anomaly in the sif this part of the PVT is an anomaly in the sense that the exchanger
    # is detached from the PV over a short distance, the transfer coefficient at the absorber (which in this case is the PV backsheet)
    # is taken as that calculated for the “main”.

    if  componentSpecs["is_anomaly"] == 1:

        if hyp['method_h_back_anomaly'] == "like_exchanger":
            var["h_back"] = hyp['h_back_prev']
        else:
            raise ValueError("Method for h_back is not well defined for anomalies")

    # Manifold
    # If this part of the PVT is a manifold, the transfer coefficient at the absorber is taken as:
    # - that calculated from the correlation for a free cylinder
    # - or that for the “main”.
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
            print('T_abs_mean = None in h_back_abs()')
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

def h_back_mean(componentSpecs,stepConditions,var,hyp):
    """Calculates the mean h_back between the value at iteration n-1 and the value calculated at iteration n and stores it in var["h_back"]
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    Returns:
        None
    """

    old_h_back = var["h_back"]
    h_back_abs(componentSpecs,stepConditions,var,hyp)

    new_h_back = var["h_back"]

    var["h_back"] = (old_h_back+new_h_back)/2

def h_back_tube(componentSpecs,stepConditions,var,hyp):
    """Calculates the back heat transfer coefficient for the tube and stores it in var["h_back_tube"]

    - If the tube is well connected to the absorber (conductive), the coefficient is the same as the back heat transfer coefficient for the absorber
    - If not,
    
        - the characteristic length is the height of the tube H_tube
        - if it is a manifold and hyp['inlet_manifold_in_wind'] == 1, the coefficient is calculated as the cubic mean of ht.conv_external and the free convection coefficient for a cylinder
        - if not, the coefficient is calculated with the free convection coefficient for a cylinder
    
    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
        
    Returns:
        None"""
    
    if hyp['method_h_back_tube'] == 'CFD':
        get_CFD_value(componentSpecs, stepConditions, var, hyp, 'h_back_tube', 'phi_tube', 'T_tube_mean', 'T_amb')
        return

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

# RADIATIVE

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

    if hyp['method_h_rad_back_tube'] == 'CFD':
        get_CFD_value(componentSpecs, stepConditions, var, hyp, 'h_rad_back_tube', 'phi_rad_back_tube', 'T_abs_mean', 'T_amb')
        return

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

def h_rad_tube_abs(componentSpecs,stepConditions,var,hyp):
    """Calculates the radiation heat transfer coefficient between the tube and the absorber and stores it in var["h_rad_tube_abs"]
    
    If the tube is in contact with the absorber, it is set to 0 (conduction only). Otherwise, it is calculated as follows:

    $$h_{rad,f} = \epsilon_{hx} \sigma (T_{tube}+T_{abs})(T_{tube}^2+T_{abs}^2)$$
    
    A multiplying factor is optional in hypotheses file as "coeff_h_rad_tube_abs"

    Args:
        componentSpecs (dict): dictionary containing the PVT parameters
        stepConditions (dict): dictionary containing the meteo inputs
        var (dict): dictionary containing the variables
        hyp (dict): dictionary containing the hypotheses
    
    Returns:
        None
    """

    if componentSpecs["l_c"] > 0. :
        h=0.

    else:

        sigma = hyp["sigma"]

        T_ref = var["T_tube_mean"]
        T_B = var["T_Base_mean"]

        eps = componentSpecs["eps_hx_top"]

        h = eps*sigma*(T_ref+T_B)*(T_ref**2+T_B**2)

    var["h_rad_tube_abs"] = hyp["coeff_h_rad_tube_abs"]*h

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
