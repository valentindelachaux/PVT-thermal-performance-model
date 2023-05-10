import math
import copy
import pandas as pd
import numpy as np
import heat_transfer as bht

from CoolProp.CoolProp import PropsSI

# %run C:\Users\BU05\Documents\Modele1D_Type560\Type560.py

# par

# "sigma", "eta_nom","Eff_T","T_ref","Eff_G","G_ref","X_rad","X_corr","W","D_tube","N_harp",
# "L_tube","tau_alpha","eps","k_abs","lambd_abs","h_fluid","R_top","R_inter","R_b","C_B","G",
# "Gp","T_sky","T_amb","T_back","T_fluid_in0","Cp","mdot","theta"

# Function for Excel
# Search a cell by its nom_variable and return a string like "B3"
def find_cell_by_name(wb,nom_variable):
    my_range = wb.defined_names[nom_variable]
    ch = my_range.value
    ch2 = ch.split("!")[1]
    ch3 = ch2.replace("$","")
    return ch3
    


## Relationship between wind and R_top
# Modify the dictionary par by updating the wind speed
# To complete from Excel file
def change_u(par,par_p,wind_speed):
    par_p["u"] = wind_speed
    
    a_w = par["a_htop"]
    b_w = par["b_htop"]

    new_h_wind = a_w*wind_speed+b_w

    par["h_top"]=new_h_wind
    par["R_t"]= par["R_top"] + 1/par["h_top"]

# return tanh or 1/tanh
def tanh_or_inverse(arg):
    return math.tanh(arg)

def h_fluid(par,par_p,var,hyp):
    D_tube = par["D_tube"]

    if par["is_manifold"] == 1 and par["geometry"] == "harp":
        mdot = par_p["mdot"]/2
    else:
        mdot = par_p["mdot"]
    
    N_harp = par["N_harp"]

    T_fluid = par_p["T_fluid_in0"]

    p_fluid = hyp["p_fluid"]
    fluid = hyp["fluid"]
    glycol_rate = hyp["glycol_rate"]

    k_fluid = PropsSI('L', 'P', p_fluid, 'T', T_fluid, f'INCOMP::{fluid}[{glycol_rate}]')
    rho_fluid = PropsSI('D', 'P', p_fluid, 'T', T_fluid, f'INCOMP::{fluid}[{glycol_rate}]')
    mu_fluid = PropsSI('V', 'P', p_fluid, 'T', T_fluid, f'INCOMP::{fluid}[{glycol_rate}]')

    flow_rate_per_riser = (mdot/N_harp)/rho_fluid # en m3/s
    tube_section = math.pi*(D_tube/2)**2

    fluid_speed = flow_rate_per_riser/tube_section

    Re = (rho_fluid*fluid_speed*D_tube)/mu_fluid

    Pr_fluid = 7.

    if Re<2000:
        Nu_fluid = 0.7*0.023*(Re**0.8)*(Pr_fluid**0.4)
    else:
        Nu_fluid = 0.023*(Re**0.8)*(Pr_fluid**0.4)

    var["h_fluid"]  = (k_fluid/D_tube)*Nu_fluid

def h_top(par,par_p,var,hyp):

    if par["is_manifold"] == 1:
        var['h_top'] = hyp['h_top_man'] 

    else:
        T_PV = var["T_PV"]
        T_amb = par_p["T_amb"]

        if par["orientation"]=="portrait":
            L_c = par["L_pan"]
        else:
            L_c = par["w_pan"]

        h_free = hyp["coeff_h_top_free"]*bht.top_h_simple(T_PV,T_amb,hyp["theta"],L_c)



        h_forced = hyp["coeff_h_top_forced"]*bht.h_top_forced(T_PV,T_amb,par_p["u"],L_c)

        var["h_top"] = (h_free**3 + h_forced**3)**(1/3)


def h_back(par,par_p,var,hyp):


    # Anomaly
    if  par["is_anomaly"] == 1:
        var["h_back"] = hyp['h_back_prev']

    # Manifold
    elif par["is_manifold"] == 1 :
        if hyp['calc_h_back_manifold'] == 1:
            L_c = par['H_tube']

            if par["insulated"] == 1:
                R_2 = par["R_2"]
                h_back = var["h_back"]
                T_back = par_p["T_back"]

                T_ref = var["T_ins_mean"]
                var["T_ins_mean"] = T_ref + (R_2/(R_2+1/h_back)) * (T_back - T_ref)
                
            else:
                T_ref = var["T_abs_mean"]

            if par["tube_geometry"]=="rectangular" or par["tube_geometry"]=="square":
                res = bht.back_h_mixed(T_ref,par_p["T_back"],par_p["u_back"],hyp["theta"],L_c)
            else:
                res = bht.back_h_cylinder(T_ref,par_p["T_back"],L_c)
                
            var["h_back"] = res
            
        else:
            var["h_back"] = hyp['h_back_prev']

    # Echangeur 
    else:

        if par["orientation"] == 'portrait':
            L_c = par["L_abs"]
        else:
            L_c = par["w_abs"]

        # Si y a une erreur
        if var["T_abs_mean"]==None:
            print('T_abs_mean = None in h_back()')
            var["h_back"] = 0.5
        
        # Cas avec ailettes

        elif par["fin_0"] >= 1 or par["fin_1"] >= 1 or par["fin_2"] >= 1:

            D = par["D"]

            if par["N_ail"]<= 24:
                var["h_back"] = hyp["coeff_h_back"]*bht.back_h_simple(var["T_abs_mean"],par_p["T_back"],hyp["theta"],L_c)
            else:
                var["h_back"] = hyp["coeff_h_back"]*bht.back_h_fins(var["T_abs_mean"],par_p["T_back"],hyp["theta"],L_c,D,par["Heta"])

        # Cas sans ailette

        else:
            # theta est l'inclinaison du panneau par rapport à l'horizontale
            
            # Cas isolé
            if par["insulated"] == 1:
                R_2 = par["R_2"]
                h_back = var["h_back"]
                T_back = par_p["T_back"]

                T_ref = var["T_ins_mean"]
                var["T_ins_mean"] = T_ref + (R_2/(R_2+1/h_back)) * (T_back - T_ref)
            
            # Cas non-isolé
            else:
                T_ref = var["T_abs_mean"]

            h1 = hyp["coeff_h_back"]*bht.back_h_mixed(T_ref,par_p["T_back"],par_p["u_back"],hyp["theta"],L_c)

            # res = hyp["coeff_h_back"]*bht.back_h_simple(T_ref,par_p["T_back"],hyp["theta"],L_c)
            # print('res',res)

            if h1 == None:
                print('res = None in h_back()')
                print('L_c',L_c)
                print('T_abs',var["T_abs_mean"])
                print('theta',hyp["theta"])
                print('T_back',par_p["T_back"])
                print('h_back_calculated',res)
            else:

                T_back_changed = hyp["T_back_changed"]
                if T_back_changed > 0:
                    T_back_changed += 273.15
                    h2 = hyp["coeff_h_back"]*bht.back_h_mixed(T_ref,T_back_changed,par_p["u_back"],hyp["theta"],L_c)
                    h1 = (h2*T_back_changed)/par_p["T_back"]
                else:
                    pass                          

                var["h_back"] = h1
                # print("var['h_back']",var["h_back"])

def h_back_mean(par,par_p,var,hyp):

    old_h_back = var["h_back"]
    h_back(par,par_p,var,hyp)

    new_h_back = var["h_back"]

    if abs(new_h_back - old_h_back) > 0.2:
        var["h_back"] = (old_h_back+new_h_back)/2
    else:
        pass


def h_back_tube(par,par_p,var,hyp):

    # Cas où le tuyau est bien lié à l'absorbeur (conductif)
    if par["l_c"] > 0. :
        var["h_back_tube"] = var["h_back"]

    # Cas où le tuyau est décollé de l'échangeur (radiatif seul)
    # Alors on calcule dans tous les cas comme si c'était unn tube cylindrique
    else:

        L_c = par['H_tube']

        if par["insulated"] == 1 and par_p["compt"] >= 1:

            T_ref = var["T_ins_tube_mean"]

            # R_2 = par["R_2"]
            # h_back = var["h_back"]
            # T_back = par_p["T_back"]

            # T_ref = var["T_ins_mean"]
            # var["T_ins_mean"] = T_ref + (R_2/(R_2+1/h_back)) * (T_back - T_ref)
            
        else:
            T_ref = var["T_tube_mean"]

        # if par["tube_geometry"]=="rectangular" or par["tube_geometry"]=="square":
        #     res = bht.back_h_mixed(T_ref,par_p["T_back"],par_p["u_back"],hyp["theta"],L_c)
        # else:
        res = bht.back_h_cylinder(T_ref,par_p["T_back"],L_c)
                
        var["h_back_tube"] = res

def h_back_fins(par,par_p,var,hyp):
    if hyp["h_back_fins_calc"] == "tube":
        var["h_back_fins"] = var["h_back_tube"]
    elif hyp["h_back_fins_calc"] == "abs":
        var["h_back_fins"] = var["h_back"]
    elif hyp["h_back_fins_calc"] == "TS":
        L_c = par["L_fin"]
        D = par["D"]
        var["h_back_fins"] = hyp["coeff_h_back"]*bht.back_h_fins(var["T_tube_mean"],par_p["T_back"],hyp["theta"],L_c,D,par["Heta"])
    else:
        pass

def h_rad_back_tube(par,par_p,var,hyp):
    sigma = hyp["sigma"]
    T_rad_back = par_p["T_back"] # hypothèse T_amb = T_back   
    if par["insulated"] == 1 and par_p["compt"] >= 1:
        T_ref = var["T_ins_tube_mean"]
        eps = par["eps_ins"]
    else: 
        T_ref = var["T_tube_mean"]
        eps = par["eps_he"]

    h = eps*sigma*(T_ref+T_rad_back)*(T_ref**2+T_rad_back**2)
    var["h_rad_back_tube"]=h

def h_rad_back(par,par_p,var,hyp):

    sigma = hyp["sigma"]
    T_rad_back = par_p["T_back"] # hypothèse T_amb = T_back   
    if par["insulated"] == 1:
        T_ref = var["T_ins_mean"]
        eps = par["eps_ins"]
    else: 
        T_ref = var["T_abs_mean"]
        eps = par["eps_he"]
   
    T_rad_back_changed = hyp["T_rad_back_changed"]

    if T_rad_back_changed > 0:
        T_rad_back_changed += 273.15
        sigma = hyp["sigma"]

        h2 = eps*sigma*(T_ref+T_rad_back_changed)*(T_ref**2+T_rad_back_changed**2)
        h1 = (h2*T_rad_back_changed)/T_rad_back
    else:
        h1 = eps*sigma*(T_ref+T_rad_back)*(T_ref**2+T_rad_back**2)

    var["h_rad_back"]=h1

    # var["h_rad_back"] = 1E-12

def h_rad_f(par,par_p,var,hyp):
    # var["h_rad_f"] = hyp["h_rad_f"]

    if par["l_c"] > 0. :
        h=0.

    else:

        sigma = hyp["sigma"]

        T_ref = var["T_tube_mean"]
        T_B = var["T_Base_mean"]

        eps = par["eps_he"]

        h = eps*sigma*(T_ref+T_B)*(T_ref**2+T_B**2)

    var["h_rad_f"] = h

# Pour les ailettes, on ne prend pas en compte la composante radiative h_rad_back

def Biot(lambd,k,h,delta):
    return ((lambd*h)/k)*(1+lambd/delta)

def Bi_f0(par,var):
    h_back = var["h_back_fins"]
    var["Bi_f0"] = Biot(par["lambd_ail"],par["k_ail"],h_back,par["delta_f0"])

def Bi_f1(par,var):
    h_back = var["h_back_fins"]
    var["Bi_f1"] = Biot(par["lambd_ail"],par["k_ail"],h_back,par["delta_f1"])

def Bi_f2(par,var):
    h_back = var["h_back_fins"]
    var["Bi_f1"] = Biot(par["lambd_ail"],par["k_ail"],h_back,par["delta_f2"])

def Bi_f3(par,var):
    h_back = var["h_back_fins"]
    var["Bi_f3"] = Biot(par["lambd_ail"],par["k_ail"],h_back,par["delta_f3"])

##### Variables

## PV production

# Radiation heat transfer coefficient using equation 560.3
def h_rad(par,par_p,var,hyp):
    
    eps = par["eps"]
    sigma = hyp["sigma"]
    T_sky = par_p["T_sky"]

    T_PV = var["T_PV"]

    h = eps*sigma*(T_PV+T_sky)*(T_PV**2+T_sky**2)
    var["h_rad"]=h
    #var["h_rad"]=0.00001

# Depends on T_PV
def X_celltemp(par,var):
    Eff_T = par["Eff_T"]
    T_ref = par["T_ref"]


    T_PV = var["T_PV"]

    X = 1+Eff_T*(T_PV-T_ref)

    var["X_celltemp"]=X

def eta_PV(par,par_p,var):
    
    eta_nom = par["eta_nom"]
    G = par_p["G"]
    X_rad = par["X_rad"]
    X_corr = par["X_corr"]

    #T_PV = var["T_PV"]
    X_celltemp = var["X_celltemp"]

    eta = eta_nom*X_celltemp*X_rad*X_corr
    var["eta_PV"] = eta
    #var["eta_PV"] = 0.15


# net absorbed solar radiation (total absorbed - PV power production)
def S(par,par_p,var):
    tau_alpha = par["tau_alpha"]
    G = par_p["G"]

    #T_PV = var["T_PV"]
    eta_PV = var["eta_PV"]

    S = tau_alpha*G*(1-eta_PV)

    var["S"] = S

def Fp(par, var):
    R_inter = par["R_inter"]
    R_t = par["R_top"]+1/var["h_top"]

    #T_PV = var["T_PV"]
    h_rad = var["h_rad"]

    Fp = 1/(h_rad*R_inter+R_inter/R_t+1)
    var["Fp"] = Fp

# Plus utilisé apparemment
def gamma(par):
    alpha = par["alpha_ail"]
    beta = par["beta_ail"]
    a = par["lambd_ail"]
    L_a = par["Heta"]

    arg = (alpha*L_a)/a
    numerateur = (alpha/a)*math.sinh(arg) + ((beta*alpha)/a)*math.cosh(arg)
    denominateur = math.cosh(arg) + beta*math.sinh(arg)

    gamma = (numerateur/denominateur)
    par["gamma"] = gamma

def gamma0int(N,L_fin,lambd,k,delta,delta_int,L_tube,h):
    Bi = Biot(lambd,k,h,delta)

    alpha = math.sqrt(2*Bi)
    beta = math.sqrt(Bi/2)*(1/(1+lambd/delta))
    arg = (alpha*L_fin)/lambd

    numerateur = (alpha/lambd)*math.sinh(arg) + ((beta*alpha)/lambd)*math.cosh(arg)
    denominateur = math.cosh(arg) + beta*math.sinh(arg)

    return k*(numerateur/denominateur)*((lambd*N*delta_int)/L_tube)

def gamma_0_int(par,var):

    var["gamma_0_int"] = gamma0int(par["N_f0"],par["L_f0"],par["lambd_ail"],par["k_ail"],par["delta_f0"],par["delta_f0_int"],par["L_tube"],var["h_back_fins"])

def gamma1int(N,L_fin,lambd,k,delta,delta_int,L_tube,h):
    
    Bi = Biot(lambd,k,h,delta)

    return 2*k*((lambd*N*delta_int)/L_tube)*math.tanh(math.sqrt(2*Bi)*(L_fin/lambd))*(math.sqrt(2*Bi)/lambd)

def gamma_1_int(par,var):

    var["gamma_1_int"] = par["coeff_f1"]*gamma1int(par["N_f1"],par["L_f1"],par["lambd_ail"],par["k_ail"],par["delta_f1"],par["delta_f1_int"],par["L_tube"],var["h_back_fins"])

def gamma_2_int(par,var):

    Bi = var["Bi_f2"]
    a = par["lambd_ail"]
    delta = par["delta_f2"]

    alpha = math.sqrt(2*Bi)
    beta = math.sqrt(Bi/2)*(1/(1+a/delta))

    L_a = par["L_f2"]
    N_ail = par["N_f2"]
    k = par["k_ail"]

    L_tube = par["L_tube"]

    arg = (alpha*L_a)/a
    numerateur = (alpha/a)*math.sinh(arg) + ((beta*alpha)/a)*math.cosh(arg)
    denominateur = math.cosh(arg) + beta*math.sinh(arg)

    delta_f2 = par["delta_f2"]
    
    gamma_int = k*(numerateur/denominateur)*((a*N_ail*delta_f2)/(L_tube*delta_f2))

    var["gamma_2_int"] = gamma_int

def j(par,var):
    R_inter = par["R_inter"]
    R_b = par["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    h_rad_f = var["h_rad_f"]

    Fprime = var["Fp"]
    
    # j = 1/(R_inter*Fprime)+1/(R_b*Fprime)-1/R_inter+h_rad_f/Fprime

    j = 1/(R_inter*Fprime)+1/(R_b*Fprime)-1/R_inter

    if par["fin_2"] == 1:

        gamma_int = var["gamma_2_int"]

        j += (gamma_int)/Fprime

    var["j"] = j

def b(par,par_p,var):
    T_sky = par_p["T_sky"]
    T_amb = par_p["T_amb"]
    R_t = par["R_top"]+1/var["h_top"]
    T_back = par_p["T_back"]
    R_b = par["R_2"] + 1/(var["h_back"]+var["h_rad_back"])

    h_rad = var["h_rad"]
    S = var["S"]
    Fprime = var["Fp"]

    h_rad_f = var["h_rad_f"]

    # if h_rad_f != 0:
    #     T_tube_mean = var["T_tube_mean"]
    #     b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime) + (h_rad_f*T_tube_mean)/Fprime
    # else:
    #     b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime)

    b = S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime)


    if par["fin_2"]==1:
        gamma_int = var["gamma_2_int"]

        b += (gamma_int*T_back)/Fprime

    var["b"] = b

def m(par, var):
    lambd_abs = par["lambd_abs"]
    k_abs = par["k_abs"]

    Fprime = var["Fp"]

    j = var["j"]

    m = math.sqrt((Fprime*j)/(k_abs*lambd_abs))

    var["m"] = m   

# Need the absorber's fin base temperature T_B - function not used
def qp_fin(par, var):
    lambd_abs = par["lambd_abs"]
    k_abs = par["k_abs"]

    L_af = par["L_af"]

    #T_PV = var["T_PV"]
    T_B = var["T_Base_mean"]

    j = var["j"]
    b = var["b"]

    m = var["m"]
    
    q = k_abs*lambd_abs*m*((b/j)-T_B)*tanh_or_inverse(m*L_af)
    var["qp_fin"] = q

def c0(par, var):
    p_int_tube = par["p_int_tube"]
    h_fluid = var["h_fluid"]
    C_B = par["C_B"]   
    chi = 1/(h_fluid*p_int_tube)
 

    h_rad_f = var["h_rad_f"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    # vérifier l'homogénéité
    var["c0"] = 1/chi + h_rad_f*p_ext_tube_rad + C_B + gamma

def c2(par, var):

    C_B = par["C_B"]  

    h_rad_f = var["h_rad_f"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    c0 = var["c0"]

    var["c2"] = (C_B + gamma + h_rad_f*p_ext_tube_rad)/c0

def e1(par,var):

    C_B = par["C_B"]
    p_ext_tube_rad = par["p_ext_tube_rad"]
    h_rad_f = var["h_rad_f"]

    c0 = var["c0"]

    var["e1"] = (1/c0)*(C_B+p_ext_tube_rad*h_rad_f)

def e2(par,var):

    p_int_tube = par["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    c0 = var["c0"]

    var["e2"] = (1/c0)*(1/chi)

def e3(par,var):

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    c0 = var["c0"]

    var["e3"] = (1/c0)*gamma

def f0(par,var):

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    var["f0"] = par["C_B"] + par["h_rad_f"]*par["p_ext_tube_rad"] + gamma

def b1(par, var):
    C_B = par["C_B"]  
    h_rad_f = var["h_rad_f"]

    c2 = var["c2"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]

    c0 = var["c0"]

    # A CHECK c0 ou pas dans la multiplication
    # var["b1"] = (1/(C_B + h_rad_f*p_ext_tube_rad)) * (1/(1-c2))
    # var["b1"] = (1/(C_B+h_rad_f*p_ext_tube_rad))*c0

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    e1 = var["e1"]

    # print('b1')
    # print('e1',e1)
    # print('gamma',gamma)

    var["b1"] = 1/(C_B + h_rad_f*p_ext_tube_rad - (C_B + gamma + h_rad_f*p_ext_tube_rad)*e1)
                   
def b2(par, var):
    C_B = par["C_B"]  
    h_rad_f = var["h_rad_f"]
    


    p_int_tube = par["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    c2 = var["c2"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    e2 = var["e2"]

    # Attention la Y AVAIT UNE ERREUR DE SIGNE
    # var["b2"] = (c2/(1-c2))*(1/(C_B + h_rad_f*p_ext_tube_rad))*(1/chi)
    # var["b2"] = (1/(C_B+h_rad_f*p_ext_tube_rad))*(1/chi)
    var["b2"] = (C_B + gamma + h_rad_f*p_ext_tube_rad)*e2*var["b1"]

def b3(par, var):

    C_B = par["C_B"]  
    h_rad_f = var["h_rad_f"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    e3 = var["e3"]

    # c'est vérifié c'est égal
    # var["b3"] = (-1/(C_B + h_rad_f*p_ext_tube_rad))*gamma
    # var["b3"] =  (-1/(C_B + h_rad_f*p_ext_tube_rad))*gamma

    var["b3"] = ( (C_B + gamma + h_rad_f*p_ext_tube_rad)*e3 - gamma)*var["b1"]

def d1(par, var):
    d0 = 1/var["c0"]

    p_int_tube = par["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_f = var["h_rad_f"]

    e1 = var["e1"]

    # var["d1"] = (-gamma*d0)/R_B + h_rad_f*(1-(d0/R_B))
    var["d1"] = -e1*(gamma+h_rad_f*p_ext_tube_rad) + h_rad_f*p_ext_tube_rad

def d2(par, var):

    d0 = 1/var["c0"]

    p_int_tube = par["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_f = var["h_rad_f"]

    # var["d2"] = (-gamma*d0)/chi - (h_rad_f*d0)/chi

    e2 = var["e2"]

    var["d2"] = -e2*(gamma+h_rad_f*p_ext_tube_rad)

def d3(par, var):

    d0 = 1/var["c0"]

    p_int_tube = par["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    h_rad_f = var["h_rad_f"]

    # var["d3"] = -gamma**2*d0 + gamma - (h_rad_f*d0)/(1/gamma)

    e3 = var["e3"]
    var["d3"] = -e3*(gamma+h_rad_f*p_ext_tube_rad) + gamma

def KTE_Bt(par,par_p,var):
    lambd_abs = par["lambd_abs"]
    k_abs = par["k_abs"]
    W = par["W"]
    L_af = par["L_af"]
    l_B = par["l_B"]
    l_c = par["l_c"]
    p_int_tube = par["p_int_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]

    R_inter = par["R_inter"]
    R_t = par["R_top"]+1/var["h_top"]
    R_b = par["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    h_fluid = var["h_fluid"]

    T_sky = par_p["T_sky"]
    T_amb = par_p["T_amb"]
    T_back = par_p["T_back"]

    C_B = par["C_B"]

    #T_PV = var["T_PV"]
    h_rad = var["h_rad"]
    S = var["S"]
    Fprime = var["Fp"]

    j = var["j"]
    b = var["b"]
    m = var["m"] 
    
    iota = par["iota"]

    # K = b2 * ( -D_tube*Fprime*((l_c/D_tube)*(h_rad+1/R_t)+(iota/D_tube)/(R_b*Fprime))-2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af) )
    # T = b1 * ( D_tube*Fprime*((l_c/D_tube)*(h_rad+1/R_t)+(iota/D_tube)/(R_b*Fprime))+2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af))
    # E = D_tube*Fprime*((l_c/D_tube)*(S+h_rad*T_sky+T_amb/R_t)+((iota/D_tube)*(1-b3)*T_back)/(R_b*Fprime))+2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*((b/j) - b3*T_back) + (l_c/R_t)*(Fprime-1)*b3*T_back

    b1 = var["b1"]
    b2 = var["b2"]
    b3 = var["b3"]

    # T_PV - T_B 
    # if par_p["compt"] <=2:
    #     g4 = l_c*Fprime*(S+h_rad*T_sky+T_amb/R_t)
    #     g3 = (l_c/R_inter)*(Fprime-1)*b3
    #     g2 = (l_c/R_inter)*(Fprime-1)*b2
    #     g1 = (l_c/R_inter)*(Fprime-1)*b1
    # else:
    #     T_abs = var["T_abs_mean"]
    #     g4 = l_c*Fprime*(S+h_rad*T_sky+T_amb/R_t + T_abs/R_inter)
    #     g3 = -(l_c/R_inter)*b3
    #     g2 = -(l_c/R_inter)*b2
    #     g1 = -(l_c/R_inter)*b1

    g4 = l_c*Fprime*(S+h_rad*T_sky+T_amb/R_t)
    g3 = (l_c/R_inter)*(Fprime-1)*b3
    g2 = (l_c/R_inter)*(Fprime-1)*b2
    g1 = (l_c/R_inter)*(Fprime-1)*b1

    # T_B - T_back
    h4 = 0.
    h3 = (iota/R_b)*(1-b3)
    h2 = (-iota/R_b)*b2
    h1 = (-iota/R_b)*b1

    # 2q'_absfin
    i4 = 2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*(b/j)
    i3 = -2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*b3
    i2 = -2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*b2
    i1 = -2*k_abs*lambd_abs*m*tanh_or_inverse(m*L_af)*b1

    K =  g2+h2+i2
    T = -(g1+h1+i1)
    E = (g3+h3+i3)*T_back + (g4+h4+i4)

    var["Ka_Bt"] = K
    var["Th_Bt"] = T
    var["Ep_Bt"] = E 

def KTE_tf(par,par_p,var):

    Ka_Bt = var["Ka_Bt"]
    Th_Bt = var["Th_Bt"]
    Ep_Bt = var["Ep_Bt"]

    b1 = var["b1"]
    b2 = var["b2"]
    b3 = var["b3"]

    d0 = 1/var["c0"]
    d1 = var["d1"]
    d2 = var["d2"]
    d3 = var["d3"]

    C_B = par["C_B"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    h_rad_f = var["h_rad_f"]

    R_B = 1/(C_B+p_ext_tube*h_rad_f)
    p_int_tube = par["p_int_tube"]
    h_fluid = var["h_fluid"]
    chi = 1/(h_fluid*p_int_tube)
    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    T_back = par_p["T_back"]

    # var["Ka_tf"] = Ka_Bt - (gamma*b2)/c3
    # var["Th_tf"] = 1/c3 + gamma*b1 + Th_Bt
    # var["Ep_tf"] = Ep_Bt + gamma*(1-b3)*T_back

    var["Ka_tf"] = d1*b2 + d2 + Ka_Bt
    var["Th_tf"] = 1 - d1*b1 + Th_Bt
    # calculé : puis test avec changement des signes ci-dessus
    #     var["Th_tf"] = 1 - d1*b1 + Th_Bt
    var["Ep_tf"] = Ep_Bt + (d1*b3 + d3) * T_back

    # var["Ka_tf"] = (d0*b2)/(chi*R_B) + (d0/chi**2) - 1/chi
    # var["Th_tf"] = 1 - (d0*b1)/(chi*R_B)
    # var["Ep_tf"] = ( (d0*b3)/(chi*R_B) + (d0*gamma)*chi )*T_back

def Cp(par,par_p,var,hyp):
    T_m = (var["T_fluid_in"]+var["T_fluid_out"])/2

    p_fluid = hyp["p_fluid"]
    fluid = hyp["fluid"]
    glycol_rate = hyp["glycol_rate"]

    var["Cp"] = PropsSI('C','P', p_fluid, 'T', T_m, f'INCOMP::{fluid}[{glycol_rate}]')


def ab_f(par,par_p,var):
    N_harp = par["N_harp"]
    mdot = par_p["mdot"]
    Cp = var["Cp"]
    
    Ka_tf = var["Ka_tf"]
    Th_tf = var["Th_tf"]
    Ep_tf = var["Ep_tf"]

    a = (N_harp/(mdot*Cp))*(Ka_tf/Th_tf)
    b = (N_harp/(mdot*Cp))*(Ep_tf/Th_tf)

    var["a_f"] = a
    var["b_f"] = b

    # if par["fin_0"] == 1:
    #     gamma_0_int = var["gamma_0_int"]
    # else:
    #     gamma_0_int = 0

    # if par["fin_1"] == 1:
    #     gamma_1_int = var["gamma_1_int"]
    #     # print(gamma_1_int)
    # else:
    #     gamma_1_int = 0
        

    # k = par["k_ail"]
    # C_B_f = (p_ext_tube*par["k_riser"])/par["lambd_riser_back"]
    # h_fluid = var["h_fluid"]

    # chi = 1/(h_fluid*p_int_tube)+1/C_B_f

    # gamma_back = p_ext_tube/(R_2+1/h_back)

    # zeta = (gamma_back + gamma_1_int + gamma_0_int)/(1+chi*(gamma_back+gamma_1_int+gamma_0_int))
    # # print(zeta)

    # a += (-N_harp/(mdot*Cp))*zeta
    # b += (N_harp/(mdot*Cp))*zeta*T_back

# Eq. 560.36
def T_fluid_out(par, T_fluid_in,var):

    a = var["a_f"]
    b = var["b_f"]

    L_tube = par["L_tube"]

    res = (T_fluid_in+(b/a))*math.exp(a*L_tube) - b/a
    var["T_fluid_out"] = res

# Eq. 560.38
def q_tube_fluid(par,par_p,T_fluid_in,var):
    
    N_harp = par["N_harp"]
    L = par["L_tube"]
    mdot = par_p["mdot"]
    Cp = var["Cp"]    
    
    T_f_out = var["T_fluid_out"]
    res = (mdot*Cp*(T_f_out-T_fluid_in))/(L*N_harp)

    var["q_tube_fluid"] = res

def q_Base_tube(par,var):
    Ka_Bt = var["Ka_Bt"]
    Th_Bt = var["Th_Bt"]
    Ep_Bt = var["Ep_Bt"]

    T_fluid = var["T_fluid_mean"]

    q_tube_fluid = var["q_tube_fluid"]

    var["q_Base_tube"] = -Th_Bt*q_tube_fluid + Ka_Bt*T_fluid + Ep_Bt

# Eq. 560.40
def T_fluid_mean(par,T_fluid_in,var):

    L_tube = par["L_tube"]

    h_back_tube = var["h_back_tube"]
    if h_back_tube == None:
        print(var["T_tube_mean"])
        h_back_tube = 3.

    a = var["a_f"]
    b = var["b_f"]

    res = ((T_fluid_in+(b/a))/(a*L_tube))*math.exp(a*L_tube) - (T_fluid_in+(b/a))/(a*L_tube) - b/a
    var["T_fluid_mean"] = res

# Eq. 560.28 -> calculate the mean base temperature
def T_Base_mean(par, par_p,var): #T_fluid has already been used for q_f_p and T_f_mean calculations

    # C_B = par["C_B"]

    # p_int_tube = par["p_int_tube"]
    # h_fluid = var["h_fluid"]
    # chi = 1/(h_fluid*p_int_tube)

    # h_back = var["h_back"]+var["h_rad_back"]
    # p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    # R_2 = par["R_2"]
    # gamma_back = p_ext_tube/(R_2+1/h_back)
    # gamma_0_int = var["gamma_0_int"]
    # gamma_1_int = var["gamma_1_int"]
    # gamma = gamma_back + gamma_0_int + gamma_1_int

    # h_rad_f = var["h_rad_f"]

    # c0 = var["c0"]
    T_fluid = var["T_fluid_mean"]
    q_tube_fluid = var["q_tube_fluid"]
    T_back = par_p["T_back"]

    # res = (1/(C_B+h_rad_f))*(c0*q_tube_fluid - (1/chi)*T_fluid - gamma*T_back)
    # var["T_Base_mean"] = res

    # res = (1/(h_fluid*p_int_tube)+1/C_B)*q_tube_fluid + T_f_mean

    b1 = var["b1"]
    b2 = var["b2"]
    b3 = var["b3"]

    res = b1*q_tube_fluid + b2*T_fluid + b3*T_back
    var["T_Base_mean"] = res

# Eq. 560.42 -> calculate the mean fin temperature
def T_absfin_mean(par,par_p,var):

    lambd_abs = par["lambd_abs"]
    k_abs = par["k_abs"]
    W = par["W"]

    L_af = par["L_af"]

    R_t = par["R_top"] + 1/var["h_top"]
    R_b = par["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    h_fluid = var["h_fluid"]

    T_sky = par_p["T_sky"]
    T_amb = par_p["T_amb"]
    T_back = par_p["T_back"]

    h_rad = var["h_rad"]
    S = var["S"]
    Fprime = var["Fp"]

    j = var["j"]
    m = var["m"]

    T_B_mean = var["T_Base_mean"]

    first_term = (S+h_rad*T_sky+T_amb/R_t+T_back/(R_b*Fprime))/j

    second_term = (T_B_mean-first_term)*tanh_or_inverse(m*L_af)/(m*L_af)

    res = first_term + second_term
    var["T_absfin_mean"] = res

# Eq. 560.43 -> calculate the mean absorber temperature
def T_abs_mean(par,var):

    W = par["W"]
    l_B = par["l_B"]
    L_af = par["L_af"]

    T_Base_mean = var["T_Base_mean"]
    T_absfin_mean = var["T_absfin_mean"]

    res = (l_B*T_Base_mean+(L_af*2)*T_absfin_mean)/W
    var["T_abs_mean"] = res

def T_tube_mean(par,par_p,var):

    e1 = var["e1"]
    e2 = var["e2"]
    e3 = var["e3"]

    T_B = var["T_Base_mean"]
    T_fluid = var["T_fluid_mean"]
    T_back = par_p["T_back"]

    var["T_tube_mean"] = e1*T_B + e2*T_fluid + e3*T_back

# Eq. 560.1 -> calculte the mean PV surface temperature
def T_PV_mean(par,par_p,var):

    R_inter = par["R_inter"]
    R_t = par["R_top"] + 1/var["h_top"]
    h_fluid = var["h_fluid"]

    T_sky = par_p["T_sky"]
    T_amb = par_p["T_amb"]
    
    h_rad = var["h_rad"]
    S = var["S"]

    T_abs_mean = var["T_abs_mean"]

    res = (1/(1/R_t+h_rad+1/R_inter))*(S+T_amb/R_t+h_rad*T_sky+(T_abs_mean/R_inter))

    var["T_PV0"] = var["T_PV"]
    var["T_PV"] = res

def T_PV_Base_mean(par,par_p,var):

    R_inter = par["R_inter"]
    R_t = par["R_top"] + 1/var["h_top"]
    h_fluid = var["h_fluid"]

    T_sky = par_p["T_sky"]
    T_amb = par_p["T_amb"]
    
    h_rad = var["h_rad"]
    S = var["S"]

    T_Base_mean = var["T_Base_mean"]

    res = (1/(1/R_t+h_rad+1/R_inter))*(S+T_amb/R_t+h_rad*T_sky+(T_Base_mean/R_inter))

    var["T_PV_Base_mean"] = res

def T_PV_absfin_mean(par,var):
    L_af = par["L_af"]
    W = par["W"]

    T_PV_Base_mean = var["T_PV_Base_mean"]
    T_PV_mean = var["T_PV"]

    var["T_PV_absfin_mean"] = (W*T_PV_mean - (W-2*L_af)*T_PV_Base_mean)/(2*L_af)


# Eq. 560.47
def Q_top_conv(par,par_p,var):

    T_PV_m = var["T_PV"]
    T_amb = par_p["T_amb"]

    R_t = par["R_top"] + 1/var["h_top"]
    W = par["W"]
    L = par["L_tube"]

    var["Q_top_conv"] = (W*L)*(T_PV_m-T_amb)/R_t

def Q_top_rad(par,par_p,var):

    h_r = var["h_rad"]
    T_PV_m = var["T_PV"]
    T_sky = par_p["T_sky"]
    W = par["W"]

    L = par["L_tube"]

    var["Q_top_rad"] = W*L*h_r*(T_PV_m-T_sky)

def Q_PV_plate(par,var):

    R_inter = par["R_inter"]
    W = par["W"]

    T_PV_m = var["T_PV"]
    T_abs_m = var["T_abs_mean"]
    L = par["L_tube"]

    var["Q_PV_plate"] = (W*L)*(T_PV_m-T_abs_m)/R_inter

def power_balance_1(par,var):
    S = var["S"]
    Q1 = var["Q_top_conv"]
    Q2 = var["Q_top_rad"]
    Q3 = var["Q_PV_plate"]
    W = par["W"]
    L = par["L_tube"]

    var["power_balance_1"] = (W*L)*S-Q1-Q2-Q3

def Q_PV_Base(par,var):

    R_inter = par["R_inter"]
    l_B = par["l_B"]

    T_PV_mB = var["T_PV_Base_mean"]
    T_Base_m = var["T_Base_mean"]
    L = par["L_tube"]

    var["Q_PV_Base"] = L*l_B*((T_PV_mB-T_Base_m)/R_inter)

def Q_PV_absfin(par,var):
    R_inter = par["R_inter"]
    L_af = par["L_af"]

    T_PV_absfin_mean = var["T_PV_absfin_mean"]
    T_absfin_mean = var["T_absfin_mean"]
    L = par["L_tube"]

    var["Q_PV_absfin"] = L*2*L_af*((T_PV_absfin_mean-T_absfin_mean)/R_inter)

def qp_PV_Base(par,var):

    R_inter = par["R_inter"]
    l_c = par["l_c"]

    T_PV_m = var["T_PV"]
    T_Base_m = var["T_Base_mean"]
    L = par["L_tube"]

    var["qp_PV_Base"] = l_c*((T_PV_m-T_Base_m)/R_inter)

def Q_absfin_back(par,par_p,var):

    R_b = par["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    L_af = par["L_af"]

    T_absfin_m = var["T_absfin_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]

    var["Q_absfin_back"] = L*2*L_af*(T_absfin_m-T_back)/R_b

def Q_absfin_tube(par,par_p,var):

    T_absfin_m = var["T_absfin_mean"]
    T_tube_m = var["T_tube_mean"]
    L = par["L_tube"]

    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    h_rad_f = var["h_rad_f"]

    var["Q_absfin_tube"] = L*p_ext_tube*h_rad_f*(T_absfin_m-T_tube_m)

def Q_tube_back(par,par_p,var):

    T_tube_m = var["T_tube_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]

    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    var["Q_tube_back"] = L*gamma*(T_tube_m - T_back)

def Q_tube_back_wo_ins_conv(par,par_p,var):
    T_tube_m = var["T_tube_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]

    h_back_tube = var["h_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    gamma_back = p_ext_tube*h_back_tube
    gamma_0_int = var["gamma_0_int"]
    gamma_1_int = var["gamma_1_int"]
    gamma = gamma_back + gamma_0_int + gamma_1_int

    var["Q_tube_back_conv"] = L*gamma*(T_tube_m - T_back)

def Q_tube_back_wo_ins_rad(par,par_p,var):
    T_tube_m = var["T_tube_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]

    h_rad_back_tube = var["h_rad_back_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    gamma_back = p_ext_tube*h_rad_back_tube

    var["Q_tube_back_rad"] = L*gamma_back*(T_tube_m - T_back)


def T_ins_tube_mean(par,var):
    R_2 = par["R_2"]
    T_tube_mean = var["T_tube_mean"]
    Q_tube_back = var["Q_tube_back"]

    var["T_ins_tube_mean"] = T_tube_mean - R_2*Q_tube_back

def Q_ins_tube_back_conv(par,par_p,var):

    T_ins_tube_m = var["T_ins_tube_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    h_back_tube = var["h_back_tube"]

    var["Q_ins_tube_back_conv"] = L*p_ext_tube*h_back_tube*(T_ins_tube_m - T_back)

def Q_ins_tube_back_rad(par,par_p,var):

    T_ins_tube_m = var["T_ins_tube_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]
    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    h_rad_back_tube = var["h_rad_back_tube"]

    var["Q_ins_tube_back_rad"] = L*p_ext_tube*h_rad_back_tube*(T_ins_tube_m - T_back)

def T_ins_absfin_mean(par,var):
    R_2 = par["R_2"]
    T_absfin_mean = var["T_absfin_mean"]
    Q_absfin_back = var["Q_absfin_back"]

    var["T_ins_absfin_mean"] = T_absfin_mean - R_2*Q_absfin_back

def Q_ins_absfin_back_conv(par,par_p,var):

    T_ins_absfin_m = var["T_ins_absfin_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]
    L_af = par["L_af"]
    h_back = var["h_back"]

    var["Q_ins_absfin_back_conv"] = L*2*L_af*h_back*(T_ins_absfin_m - T_back)

def Q_ins_absfin_back_rad(par,par_p,var):

    T_ins_absfin_m = var["T_ins_absfin_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]
    L_af = par["L_af"]
    h_back = var["h_rad_back"]

    var["Q_ins_absfin_back_rad"] = L*2*L_af*h_back*(T_ins_absfin_m - T_back)


# def T_ins_mean(par,par_p,var):
#     T_ref = var["T_abs_mean"]
#     R_2 = par["R_2"]
#     h_back = var["h_back"] + var["h_rad_back"]
#     T_back = par_p["T_back"]

#     var["T_ins_mean"] = T_ref + (R_2/(R_2+1/h_back)) * (T_back - T_ref)


def T_ins_mean(par,var):

    T_ins_tube_mean = var["T_ins_tube_mean"]
    T_ins_absfin_mean = var["T_ins_absfin_mean"]
    w_tube = par["w_tube"]
    lambd_riser_back = par["lambd_riser_back"]
    L_af = par["L_af"]
    W = par["W"]

    var["T_ins_mean"] = ((w_tube+lambd_riser_back)*T_ins_tube_mean + 2*L_af*T_ins_absfin_mean)/W

def Q_ins_conv(par,par_p,var):

    T_ins_m = var["T_ins_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]

    h_back = var["h_back"]

    var["Q_ins_conv"] = L*h_back*(T_ins_m - T_back)

def Q_ins_rad(par,par_p,var):

    T_ins_m = var["T_ins_mean"]
    T_back = par_p["T_back"]
    L = par["L_tube"]

    h_back = var["h_rad_back"]

    var["Q_ins_rad"] = L*h_back*(T_ins_m - T_back)

def Q_fluid(par,par_p,var):

    h_fluid = var["h_fluid"]
    p_int_tube = par["p_int_tube"]

    chi = 1/(h_fluid*p_int_tube)

    T_tube_m = var["T_tube_mean"]
    T_fluid_m = var["T_fluid_mean"]
    L = par["L_tube"]

    var["Q_tube_fluid"]=(L/chi)*(T_tube_m-T_fluid_m)

    N_harp = par["N_harp"]
    L = par["L_tube"]
    mdot = par_p["mdot"]
    Cp = var["Cp"]    
    
    T_f_out = var["T_fluid_out"]
    T_f_in = par_p["T_fluid_in0"]

    var["Q_fluid"] = (mdot*Cp*(T_f_out-T_f_in))/(N_harp)

def Q_Base_tube(par,var):
    
    C_B = par["C_B"]
    
    T_Base_m = var["T_Base_mean"]
    T_tube_m = var["T_tube_mean"]

    l_B = par["l_B"]
    L = par["L_tube"]

    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    h_rad_f = var["h_rad_f"]


    var["Q_Base_tube"] = L*(T_Base_m-T_tube_m)*(C_B+h_rad_f*p_ext_tube_rad)

def Q_Base_back(par,par_p,var):

    R_b = par["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    iota = par["iota"]

    T_Base_m = var["T_Base_mean"]
    T_back = par_p["T_back"]

    L = par["L_tube"]

    var["Q_Base_back"] = L*iota*(T_Base_m-T_back)/R_b

def qp_Base_back(par,par_p,var):

    R_b = par["R_2"] + 1/(var["h_back"]+var["h_rad_back"])
    iota = par["iota"]

    T_Base_m = var["T_Base_mean"]
    T_back = par_p["T_back"]

    L = par["L_tube"]

    var["qp_Base_back"] = iota*(T_Base_m-T_back)/R_b

def Q_absfins_Base(par,var):
    q = var["qp_fin"]
    L = par["L_tube"]

    var["Q_absfins_Base"] = 2*L*q

def Q_abs_back2(par,var):
    var["Q_abs_back2"] = var["Q_absfins_Base"] - var["Q_PV_plate"] + var["Q_PV_Base"]

def power_balance_3(par,var):
    Q_PV_Base = var["Q_PV_Base"]
    Q_absfins_Base = var["Q_absfins_Base"]
    Q_fluid = var["Q_tube_fluid"]
    Q_Base_back = var["Q_Base_back"]
    Q_fluid_back = var["Q_fluid_back"]

    var["power_balance_3"] = Q_PV_Base + Q_absfins_Base - (Q_fluid + Q_fluid_back) - Q_Base_back

def PB_3(par,var):
    PB3 = var["qp_PV_Base"] - var["qp_Base_back"] + 2*var["qp_fin"]-var["qp_fluid"]
    var["PB_3"] = PB3
    # print(PB3)

def Q_fluid_back(par,par_p,var):

    p_ext_tube = par["p_ext_tube"];p_ext_tube_rad = par["p_ext_tube_rad"]
    R_2 = par["R_2"]

    k = par["k_ail"]
    h_fluid = var["h_fluid"]
    p_int_tube = par["p_int_tube"]

    chi = 1/(h_fluid*p_int_tube)

    L = par["L_tube"]

    if par["fin_0"]==1:
        gamma_0_int = var["gamma_0_int"]
    else:
        gamma_0_int = 0
    if par["fin_1"]==1:
        gamma_1_int = var["gamma_1_int"]
    else:
        gamma_1_int = 0

    R_2 = par["R_2"]
    h_back_tube = var["h_back_tube"]+var["h_rad_back_tube"]
    gamma_back = p_ext_tube/(R_2+1/h_back_tube)

    T_fluid_m = var["T_fluid_mean"]
    T_back = par_p["T_back"]

    zeta = (gamma_back)/(1+chi*(gamma_back+gamma_1_int+gamma_0_int))

    var["Q_fluid_back"] = L*zeta*(T_fluid_m-T_back)


def qp_f0(par,par_p,var):

    T_fluid_m = var["T_fluid_mean"]
    T_back = par_p["T_back"]

    gamma_0_int = var["gamma_0_int"]

    var["qp_f0"] = gamma_0_int*(T_fluid_m-T_back)

def Q_f01(par,par_p,var):

    L = par["L_tube"]

    if par["fin_0"]==1:
        gamma_0_int = var["gamma_0_int"]
    else:
        gamma_0_int = 0
    if par["fin_1"]==1:
        gamma_1_int = var["gamma_1_int"]
    else:
        gamma_1_int = 0

    gamma = gamma_0_int + gamma_1_int
    
    T_tube_m = var["T_tube_mean"]
    T_back = par_p["T_back"]

    Q = L*gamma*(T_tube_m-T_back)

    var["Q_f01"] = Q


def qp_f1(par,par_p,var):

    T_fluid_m = var["T_fluid_mean"]
    T_back = par_p["T_back"]

    gamma_1_int = var["gamma_1_int"]

    var["qp_f1"] = gamma_1_int*(T_fluid_m-T_back)

def qp_f2(par,par_p,var):

    T_abs_m = var["T_abs_mean"]
    T_back = par_p["T_back"]

    gamma_2_int = var["gamma_2_int"]

    var["qp_f2"] = gamma_2_int*(T_abs_m-T_back)

def Q_f2(par,var):

    var["Q_f2"] = var["qp_f2"] * par["L_tube"]

def one_loop(par,par_p,T_fluid_in,var,hyp):

    # print(par_p["T_sky"])

    var["h_back_fins"] = hyp["h_back_fins0"]

    if par["fin_0"] == 1:
        Bi_f0(par,var)
        gamma_0_int(par,var)
    else:
        var["gamma_0_int"] = 0
    if par["fin_1"] == 1:
        Bi_f1(par,var)
        gamma_1_int(par,var)
    else:
        var["gamma_1_int"] = 0
    if par["fin_2"] == 1:
        Bi_f2(par,var)
        gamma_2_int(par,var)
    else:
        var["gamma_2_int"] = 0
    if par["fin_3"] == 1:
        Bi_f3(par,var)
        # directement dans le calcul de KTE()
    else:
        pass

    h_rad(par,par_p,var,hyp) # T_PV
    h_top(par,par_p,var,hyp)

    X_celltemp(par,var) # T_PV
    eta_PV(par,par_p,var) # X_celltemp so only T_PV
    S(par,par_p,var) # eta_PV so only T_PV


    Fp(par,var) # h_rad
    # print('Fp',var["Fp"])
    # print('h_back',var["h_back"])
    j(par,var) # Fp
    # print('j',var["j"])
    m(par,var) # Fp and j
    # print('m',var["m"])
    b(par,par_p,var) # h_rad, S and Fp
    # print('b',var["b"])
    c0(par,var)
    # print('c0',var["c0"])
    c2(par,var)

    e1(par,var)
    # print('e1',var["e1"])
    e2(par,var)
    # print('e2',var["e2"])
    e3(par,var)
    # print('e3',var["e3"])

    b1(par,var)
    # print('b1',var['b1'])
    b2(par,var)
    # print('b2',var['b2'])
    b3(par,var)
    # print('b3',var['b3'])

    KTE_Bt(par,par_p,var)
    # print("Ka_Bt",var["Ka_Bt"])
    # print("Th_Bt",var["Th_Bt"])
    # print("Ep_Bt",var["Ep_Bt"])


    d1(par,var)
    d2(par,var)
    d3(par,var)

    KTE_tf(par,par_p,var)
    # print("Ka_tf",var["Ka_tf"])
    # print("Th_tf",var["Th_tf"])
    # print("Ep_tf",var["Ep_tf"])

    ab_f(par,par_p,var) # utilise les coeffs Ka, Th et Ep calculés dans KTE_tf

    T_fluid_out(par,T_fluid_in,var)
    # print('T_fluid_in0',par_p['T_fluid_in0'],'K')
    # print('T_fluid_out',var['T_fluid_out'],'K')
    q_tube_fluid(par,par_p,T_fluid_in,var)
    # print(par_p["mdot"])
    # print(var["Cp"])
    # print('q_tube_fluid',var['q_tube_fluid'],'W/m')
    T_fluid_mean(par,T_fluid_in,var)
    # print('T_fluid_mean',var['T_fluid_mean'],'K')
    T_Base_mean(par,par_p,var)
    # print('T_Base_mean',var['T_Base_mean'],'K')
    T_tube_mean(par,par_p,var)
    # print('T_tube_mean',var['T_tube_mean'],'K')
    T_absfin_mean(par,par_p,var)
    T_abs_mean(par,var)

    Q_tube_back(par,par_p,var)
    Q_absfin_back(par,par_p,var)
    
    T_ins_tube_mean(par,var)
    T_ins_absfin_mean(par,var)
    T_ins_mean(par,var)

    Q_ins_conv(par,par_p,var)
    Q_ins_rad(par,par_p,var)

    if par["fin_0"] == 1 or par["fin_1"] == 1 or par["fin_2"] == 1:
        h_back_mean(par,par_p,var,hyp)
    else:
        h_back(par,par_p,var,hyp)
        # print('in one_loop',var["h_back"])

    h_rad_back(par,par_p,var,hyp)
    h_back_tube(par,par_p,var,hyp)
    h_rad_back_tube(par,par_p,var,hyp)
    h_back_fins(par,par_p,var,hyp)

    h_rad_f(par,par_p,var,hyp)

    T_PV_mean(par,par_p,var)
    T_PV_Base_mean(par,par_p,var)
    T_PV_absfin_mean(par,var)

    qp_PV_Base(par,var)
    qp_Base_back(par,par_p,var)
    qp_fin(par,var)

    Cp(par,par_p,var,hyp)

    # print('end of iteration',par_p["compt"])
    # print('T_tube_mean',var['T_tube_mean'])
    # print('T_Base_mean',var['T_Base_mean'])
    # print('T_PV',var['T_PV'])
    # print("h_back",var["h_back"])


def compute_power(par,par_p,var):
    Q_top_conv(par,par_p,var)
    Q_top_rad(par,par_p,var)
    Q_PV_plate(par,var)
    # Q_abs_back1(par,par_p,var)
    Q_PV_Base(par,var)
    Q_PV_absfin(par,var)
    Q_Base_back(par,par_p,var)
    Q_fluid(par,par_p,var)
    Q_tube_back(par,par_p,var)
    Q_Base_tube(par,var)
    # qp_fluid_back(par,var)
    qp_fin(par,var)
    Q_absfins_Base(par,var)
    Q_abs_back2(par,var)
    Q_fluid_back(par,par_p,var)
    Q_tube_back(par,par_p,var)
    Q_absfin_back(par,par_p,var)
    Q_absfin_tube(par,par_p,var)
    Q_tube_back_wo_ins_conv(par,par_p,var)
    Q_tube_back_wo_ins_rad(par,par_p,var)
    Q_ins_tube_back_conv(par,par_p,var)
    Q_ins_tube_back_rad(par,par_p,var)
    Q_ins_absfin_back_conv(par,par_p,var)
    Q_ins_absfin_back_rad(par,par_p,var)

    if par["fin_0"]==1 or par["fin_1"]==1:
        Q_f01(par,par_p,var)

    power_balance_1(par,var)
    power_balance_3(par,var)

    if par["fin_0"] == 1:
        qp_f0(par,par_p,var)
    if par["fin_1"] == 1:
        qp_f1(par,par_p,var)
    if par["fin_2"] == 1:
        pass

def simu_one_steady_state_all_he(par,par_p,hyp):
    
    res = {}

    if par['manifold']['input_man'] == 1:

        save_T_fluid_in0 = par_p["T_fluid_in0"]
        # Test without manifolds
        df,df_one,list_df_historic = simu_one_steady_state(par['exchanger'],par_p,hyp)

        hyp['h_back_prev'] = df_one['h_back'].values[0] # h_back de l'absorbeur
        hyp['h_top_man'] = df_one['h_top'].values[0]   

        # Inlet manifold

        df,df_one,list_df_historic = simu_one_steady_state(par['manifold'],par_p,hyp)
        res['inlet_man'] = [df.copy(),df_one.copy(),list_df_historic.copy()]

        par_p["T_fluid_in0"] = df_one["T_fluid_out"].values[0]

        # Première anomalie

        if par['anomaly1']['input_an'] == 1:

            

            df,df_one,list_df_historic = simu_one_steady_state(par['anomaly1'],par_p,hyp)
            res['anomaly1'] = [df.copy(),df_one.copy(),list_df_historic.copy()]

            par_p["T_fluid_in0"] = df_one["T_fluid_out"].values[0]
        else:
            pass

        # Echangeur 

        df,df_one,list_df_historic = simu_one_steady_state(par['exchanger'],par_p,hyp)
        res['exchanger'] = [df.copy(),df_one.copy(),list_df_historic.copy()]

        par_p['T_fluid_in0'] = df_one['T_fluid_out'].values[0]

        # Outlet manifold

        df,df_one,list_df_historic = simu_one_steady_state(par['manifold'],par_p,hyp)
        res['outlet_man'] = [df.copy(),df_one.copy(),list_df_historic.copy()]

        par_p["T_fluid_in0"] = save_T_fluid_in0

    else:
        df,df_one,list_df_historic = simu_one_steady_state(par['exchanger'],par_p,hyp)
        res['exchanger'] = [df.copy(),df_one.copy(),list_df_historic.copy()]

    df_c = pd.DataFrame()

    for str in res.keys():
        df_c = pd.concat([df_c,res[str][1]],axis=0)

    df_mean = df_c.mean()
    df_sum = df_c.sum()

    df_one = pd.DataFrame()

    for str in res["exchanger"][1].keys():
        if str in ['mdot','G','Gp','T_amb','u']:
            df_one[str] = [par_p[str]]
        elif str == "T_fluid_in":
            df_one[str] = [par_p["T_fluid_in0"]]
        elif str == "T_fluid_out":
            if par['manifold']['input_man'] == 1:
                df_one[str] = [res['outlet_man'][1]['T_fluid_out'].values[0]]
            else:
                df_one[str] = [res['exchanger'][1]['T_fluid_out'].values[0]]
        elif str in ["T_PV","T_PV_Base_mean","T_PV_absfin_mean","T_abs_mean","T_Base_mean","T_absfin_mean","T_ins_mean","T_ins_tube_mean","T_ins_absfin_mean","T_tube_mean","T_fluid_mean","h_top","h_rad","h_back","h_rad_back","h_back_tube","h_rad_back_tube","h_back_fins","h_rad_f","h_fluid","X_celltemp","eta_PV","S"]:
            av = 0
            Aire_tot = 0
            for typ in res.keys():
                if typ == 'inlet_man' or typ == 'outlet_man':
                    typ_par = "manifold"
                else:
                    typ_par = typ
                Aire = par[typ_par]["N_harp"]*par[typ_par]["W"]*par[typ_par]["L_tube"]
                Aire_tot += Aire
                av += res[typ][1][str].values[0]*Aire
            df_one[str] = [av/Aire_tot]
        elif str in ["Q_S","Q_top_conv","Q_top_rad","Q_PV_plate","Q_PV_Base","Q_PV_absfin","Q_absfins_Base","Q_Base_tube","Q_tube_fluid","Q_ins_tube_back_conv","Q_ins_tube_back_rad","Q_ins_absfin_back_conv","Q_ins_absfin_back_rad","Q_tube_back_conv","Q_tube_back_rad","Q_absfin_back","Q_absfin_tube","Q_fluid","Q_f01"]:
            sum = 0
            for typ in res.keys():
                if typ == 'inlet_man' or typ == 'outlet_man':
                    typ_par = "manifold"
                else:
                    typ_par = typ
                sum += res[typ][1][str].values[0]
            df_one[str] = [sum]

    return df_one,res


# par and var are dictionnaries
# Division of the panel in N rectangles (N=16)
def simu_one_steady_state(par,par_p,hyp):

    list_T_PV = [par_p["guess_T_PV"]]
    list_T_f_out = [par_p["T_fluid_in0"]]

    list_var_conv = []

    N_meander=par["N_meander"]

    list_df_historic = [] # liste de N df correspondant aux historiques de convergence pour chaque tranche de panneau
    df = pd.DataFrame()

    for i in range(N_meander):
        # on crée un dictionnaire var pour chaque tranche de panneau
        var = {}
        # Initialisatoin des h dans le dictionnaire var

        h_fluid(par,par_p,var,hyp)
        if i==0:
            var["h_back"] = hyp["h_back0"]
            var["h_rad_back"] = hyp["h_rad_back0"]
            var["h_back_tube"] = hyp["h_back_tube0"]
            var["h_rad_back_tube"] = hyp["h_rad_back_tube0"]
            var["h_rad_f"] = hyp["h_rad_f0"]

        else:
            var["h_back"] = h_back_prior
            var["h_rad_back"] = h_rad_back_prior
            var["h_back_tube"] = h_back_tube_prior
            var["h_rad_back_tube"] = h_rad_back_tube_prior
            var["h_rad_f"] = hyp["h_rad_f0"]

        var["Cp"] = hyp["Cp0"]

        new_guess_T_PV = list_T_PV[i]
        var["T_PV0"] = 0
        var["T_PV"] = new_guess_T_PV
        

        var['Slice'] = i

        T_f_in = list_T_f_out[i]
        var['T_fluid_in'] = T_f_in
        var["T_tube_mean"] = (var["T_PV"]+var["T_fluid_in"])/2
        
        # print('boucle ',i)
        compt = 0
        par_p["compt"] = compt
        while compt<= 8 or abs(var["T_PV"]-var["T_PV0"])>=0.001:
        # while compt<= 2 or abs(var["PB_3"])>=0.01:
            compt+=1
            par_p["compt"] = compt
            one_loop(par,par_p,T_f_in,var,hyp)

            compute_power(par,par_p,var)

            par_var = {'mdot' : par_p['mdot'],'G':par_p["G"],'Gp':par_p["Gp"],'T_amb':par_p["T_amb"],'u':par_p['u'],'h_top' : var['h_top'], 'h_back' : var['h_back'], 'h_rad_back' : var["h_rad_back"],'h_back_tube' : var['h_back_tube'], 'h_rad_back_tube' : var["h_rad_back_tube"],'h_back_fins' : var["h_back_fins"],"h_rad_f":var["h_rad_f"],'h_fluid' : var['h_fluid']}
            var_copy = copy.deepcopy(var)
            to_add_conv = {**par_var, **var_copy}
            list_var_conv.append(to_add_conv)

        one_loop(par,par_p,T_f_in,var,hyp)
        #break
        # T_PV_27(par,var)
        
        compute_power(par,par_p,var)

        par_var = {'mdot' : par_p['mdot'],'G':par_p["G"],'Gp':par_p["Gp"],'T_amb':par_p["T_amb"],'u':par_p['u'],'h_top' : var['h_top'], 'h_back' : var['h_back'], 'h_rad_back' : var["h_rad_back"], 'h_back_tube' : var['h_back_tube'], 'h_rad_back_tube' : var["h_rad_back_tube"],'h_back_fins' : var["h_back_fins"],"h_rad_f":var["h_rad_f"],'h_fluid' : var['h_fluid']}
        var_copy = copy.deepcopy(var)
        to_add = {**par_var, **var_copy}

        df_to_add = pd.DataFrame.from_dict({'row' : to_add.values()},orient='index',columns=to_add.keys())
        df = pd.concat([df,df_to_add])

        list_T_PV.append(var["T_PV"])
        list_T_f_out.append(var["T_fluid_out"])

        list_df_historic.append(pd.DataFrame(list_var_conv))

        h_back_prior = var["h_back"]
        h_rad_back_prior = var["h_rad_back"]
        h_back_tube_prior = var["h_back_tube"]
        h_rad_back_tube_prior = var["h_rad_back_tube"]

    df_mean = df.mean()
    df_sum = df.sum()

    df_one = pd.DataFrame()

    for str in df.keys():
        if str in ['mdot','G','Gp','T_amb','u']:
            df_one[str] = [par_p[str]]
        elif str == "T_fluid_in":
            df_one[str] = [par_p["T_fluid_in0"]]
        elif str == "T_fluid_out":
            df_one[str] = [list_T_f_out[N_meander]]
        elif str in ["T_PV","T_PV_Base_mean","T_PV_absfin_mean","T_abs_mean","T_Base_mean","T_absfin_mean","T_ins_mean","T_ins_tube_mean","T_ins_absfin_mean","T_tube_mean","T_fluid_mean","h_top","h_rad","h_back","h_rad_back","h_back_tube","h_rad_back_tube","h_back_fins","h_rad_f","h_fluid","X_celltemp","eta_PV","S"]:
            df_one[str] = [df_mean[str]]
        elif str in ["Q_top_conv","Q_top_rad","Q_PV_plate","Q_PV_Base","Q_PV_absfin","Q_absfins_Base","Q_Base_tube","Q_tube_fluid","Q_ins_tube_back_conv","Q_ins_tube_back_rad","Q_ins_absfin_back_conv","Q_ins_absfin_back_rad","Q_tube_back_conv","Q_tube_back_rad","Q_absfin_back","Q_absfin_tube","Q_fluid","Q_f01"]:
            df_one[str] = [par["N_harp"]*df_sum[str]]

    df_one["Q_S"] = [par["W"]*par["L_tube"]*df_mean["S"]]

    return df,df_one,list_df_historic

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

def simu_condi(par,hyp,condi_df):
    
    # Dataframe object pour la liste des résultats sur tous les points de fonctionnement
    df_res = pd.DataFrame()

    compt_test = 0

    list_df = []
    list_list_df_historic = []
    list_res = []

    for i in range(0,len(condi_df)):

        par_p = {'G':condi_df["G"][i],"T_amb":condi_df["T_amb"][i],"T_back":condi_df["T_amb"][i],"u":condi_df["u"][i], "u_back" : condi_df["u_back"][i], "T_fluid_in0":condi_df["T_fluid_in"][i]}
        change_T_sky(par_p,hyp,'TUV')  # calculate Gp and T_sky

        par_p["mdot"] = condi_df["mdot"][i]

        par_p["guess_T_PV"] = par_p["T_amb"] - 25
        # par_p["guess_T_PV"] = (par_p["T_amb"]+par_p["T_fluid_in0"])/2

        df_one,res = simu_one_steady_state_all_he(par,par_p,hyp)

        df_res = pd.concat([df_res,df_one],ignore_index=True)
        list_df.append(res['exchanger'][0])
        list_list_df_historic.append(res['exchanger'][2])
        list_res.append(res)

        compt_test+=1

    # Analysing df

    # Be careful here you have zeros for some columns

    tab = pd.DataFrame()

    df_res['DT'] = df_res['T_fluid_out'] - df_res['T_fluid_in']
    df_res['T_m'] = (df_res['T_fluid_out'] + df_res['T_fluid_in'])/2
    df_res['T_m en °C'] = df_res['T_m']-273.15

    tab['G'] = df_res['G'] # a0
    tab['-(T_m - T_a)'] = -(df_res['T_m'] - df_res['T_amb']) # a1
    # tab['-(T_m - T_a)^2'] = -(df_res['T_m'] - df_res['T_amb'])**2 # a2
    tab['-(T_m - T_a)^2'] = 0.*df_res['T_m'] # a2
    tab['-up x (T_m - T_a)'] = (df_res['u'] - 3) * tab['-(T_m - T_a)'] # a3
    # tab['Gp'] = df_res['Gp'] # a4
    tab['Gp'] = 0. * df_res['Gp'] # a4
    tab['dT_m/dt'] = 0. * df_res['Gp'] # a5
    tab['up x G'] = -(df_res['u'] - 3) * df_res['G'] # a6
    tab['up x Gp'] = -(df_res['u'] - 3) * df_res["Gp"] # a7
    tab['-(T_m - T_a)^4'] = -tab['-(T_m - T_a)']**4 # a8

    # coeff_density = [999.85,0.05332,-0.007564,0.00004323,-1.673e-7,2.447e-10]
    # coeff_density = list(reversed(coeff_density))

    coeff_Cp = [4.2184,-0.0028218,0.000073478,-9.4712e-7,7.2869e-9,-2.8098e-11,4.4008e-14]
    coeff_Cp = list(reversed(coeff_Cp))

    # df_res['density(T)'] = np.polyval(coeff_density,df_res['T_m en °C'])
    df_res['Cp(T)'] = np.polyval(coeff_Cp,df_res['T_m en °C'])*1000

    # df_res['mdot'] = df_res['density(T)']*(par_p["mdot"]/1000)

    df_res['Qdot'] = df_res['mdot']*df_res['Cp(T)']*df_res['DT']
    df_res['Qdot / AG'] = df_res['Qdot']/(par['AG'])

    matrice = tab.to_numpy()
    B = df_res['Qdot / AG'].to_numpy()

    X = np.linalg.lstsq(matrice, B, rcond = -1)

    #_ = plt.plot(df['T_m*'].to_numpy(), B, 'o', label='Original data', markersize=2)
    #_ = plt.plot(df['T_m*'].to_numpy(), np.dot(matrice,X[0]), 'o', label='Fitted line',markersize=2)
    #_ = plt.legend()
    #plt.show()

    df_res_to_concat = df_res.drop(columns=["G","Gp"])

    df_res = pd.concat([tab,df_res_to_concat],axis=1)

    return df_res,X,list_res,list_df,list_list_df_historic

def find_a_i(df,par):
    tab = pd.DataFrame()

    df['DT'] = df['T_fluid_out'] - df['T_fluid_in']
    df['T_m'] = (df['T_fluid_out'] + df['T_fluid_in'])/2

    df['T_m en °C'] = df['T_m']-273.15

    tab['G'] = df['G'] # a0
    tab['-(T_m - T_a)'] = -(df['T_m'] - df['T_amb']) # a1
    # tab['-(T_m - T_a)^2'] = -(df['T_m'] - df['T_amb'])**2 # a2
    tab['-(T_m - T_a)^2'] = 0.*df['T_m'] # a2
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

    # df['mdot'] = df['density(T)']*(par_p["mdot"]/1000)

    df['Qdot'] = df['mdot']*df['Cp(T)']*df['DT']
    df['Qdot / AG'] = df['Qdot']/(par['AG'])

    matrice = tab.to_numpy()
    B = df['Qdot / AG'].to_numpy()

    X = np.linalg.lstsq(matrice, B, rcond = -1)

    return X

def simu_condi_mpe(par,condi_df,l,h_back,L,hyp):
    
    variables = ['N_test','T_guess','G', 'Gp', 'T_amb', 'u', 'T_abs','T_fluid_in', 'T_fluid_out']
    
    # Dataframe object
    df = pd.DataFrame(columns = variables)

    sigma = hyp["sigma"]

    compt_test = 0

    for i in range(1,len(condi_df)+1):

        par_p = {}

        par_p["G"]=condi_df["G"][i]

        # T_amb = T_back
        par_p["T_amb"]=condi_df["ta"][i]+273.15

        change_T_sky(par,'TUV')

        # Back temperature = ambiant temperature
        par_p["T_back"]=par_p["T_amb"]

        # Change wind_speed in par and adapt R_t
        change_u(par,par_p,condi_df["U"][i])

        par_p["mdot"] = condi_df["mdot"][i]

        T_f_in_list = [condi_df["tin"][i]+273.15]                

        T_f_out = par_p["T_back"] + (T_f_in_list[0]-par_p["T_back"])*math.exp(-(l*h_back*L)/((par_p["mdot"]/par["N_harp"])*par["Cp"]))
        
        # len(T_f_out_list) = 1

        to_add = {'N_test' : compt_test, 'T_guess' : 293.15, 'G' : par_p["G"], 'Gp' : par_p["Gp"], 'T_amb' : par_p["T_amb"], 'h_back' : h_back, 'u' : par_p["u"], 'T_fluid_in' : T_f_in_list[0], 'T_abs' : 293.15,'T_fluid_out' : T_f_out}

        df_to_add = pd.DataFrame.from_dict({'row' : to_add.values()},orient='index',columns=to_add.keys())

        df = pd.concat([df,df_to_add])

        compt_test+=1

    # Analysing df

    # Be careful here you have zeros for some columns

    df['DT'] = df['T_fluid_out'] - df['T_fluid_in']
    df['T_m'] = (df['T_fluid_out'] + df['T_fluid_in'])/2
    df['T_m*'] = (df['T_m'] - df['T_amb'])/df['G']
    df['G x (T_m*)^2'] = df['G'] * df['T_m*']**2 * 0
    df['up x T_m*'] = (df['u'] - 3) * df['T_m*']
    df['Gp/G'] = df['Gp']/df['G']
    df['up'] = df['u'] - 3
    df['up x Gp/G'] = (df['up'] * df['Gp'])/df['G']
    df['G^3 x (T_m*)^4'] = df['G']**3 * df['T_m*']**4 * 0

    df['T_m en °C'] = df['T_m']-273.15

    coeff_density = [999.85,0.05332,-0.007564,0.00004323,-1.673e-7,2.447e-10]
    coeff_density = list(reversed(coeff_density))

    coeff_Cp = [4.2184,-0.0028218,0.000073478,-9.4712e-7,7.2869e-9,-2.8098e-11,4.4008e-14]
    coeff_Cp = list(reversed(coeff_Cp))

    df['density(T)'] = np.polyval(coeff_density,df['T_m en °C'])
    df['Cp(T)'] = np.polyval(coeff_Cp,df['T_m en °C'])*1000

    df['mdot'] = df['density(T)']*(par_p["mdot"]/1000)

    df['Qdot'] = df['mdot']*df['Cp(T)']*df['DT']
    df['Qdot / (AG x G)'] = df['Qdot']/(par['AG']*df['G'])

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

def simu_condi_mpe_big(par,par_p,condi_df,l,L,h_back_top,h_back_bottom,N_harp,hyp):
    
    variables = ['N_test', 'mdot', 'T_guess','G', 'Gp', 'T_amb', 'u', 'T_abs','T_fluid_in', 'T_fluid_out']
    
    # Dataframe object
    df_res = pd.DataFrame(columns = variables)

    sigma = hyp["sigma"]

    compt_test = 0

    for i in range(1,len(condi_df)+1):

        par_p["G"]=condi_df["G"][i]

        # T_amb = T_back
        par_p["T_amb"]=condi_df["ta"][i]+273.15

        change_T_sky(par,'TUV')

        # Back temperature = ambiant temperature
        par_p["T_back"]=par_p["T_amb"]

        # Change wind_speed in par and adapt R_t
        change_u(par,par_p,condi_df["U"][i])

        par_p["mdot"] = condi_df["mdot"][i]

        T_f_in_list = [condi_df["tin"][i]+273.15]                

        T_f_out = par_p["T_back"] + (T_f_in_list[0]-par_p["T_back"])*math.exp(-(l*L*h_back_top+l*L*h_back_bottom)/((par_p["mdot"]/N_harp)*par["Cp"]))
        
        # len(T_f_out_list) = 1

        to_add = {'N_test' : compt_test, 'mdot' : par_p["mdot"], 'T_guess' : 293.15, 'G' : par_p["G"], 'Gp' : par_p["Gp"], 'T_amb' : par_p["T_amb"], 'h_back' : h_back_top, 'u' : par_p["u"], 'T_fluid_in' : T_f_in_list[0], 'T_abs' : 293.15,'T_fluid_out' : T_f_out}

        df_to_add = pd.DataFrame.from_dict({'row' : to_add.values()},orient='index',columns=to_add.keys())

        df_res = pd.concat([df_res,df_to_add])

        compt_test+=1

    # Analysing df

    # Be careful here you have zeros for some columns

    tab = pd.DataFrame()

    df_res['DT'] = df_res['T_fluid_out'] - df_res['T_fluid_in']
    df_res['T_m'] = (df_res['T_fluid_out'] + df_res['T_fluid_in'])/2
    df_res['T_m en °C'] = df_res['T_m']-273.15

    # tab['G'] = 0. * df_res['G'] # a0
    tab['-(T_m - T_a)'] = -(df_res['T_m'] - df_res['T_amb']) # a1
    # tab['-(T_m - T_a)^2'] = -(df_res['T_m'] - df_res['T_amb'])**2 # a2
    # tab['-(T_m - T_a)^2'] = 0.*df_res['T_m'] # a2
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

    # df_res['mdot'] = df_res['density(T)']*(par_p["mdot"]/1000)

    df_res['Qdot'] = df_res['mdot']*df_res['Cp(T)']*df_res['DT']
    df_res['Qdot / AG'] = df_res['Qdot']/(par['AG'])

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


def change_T_sky(par_p,hyp,type):
    if type == "TUV":
        par_p["Gp"] = 4
        par_p["T_sky"] = par_p["T_amb"]
        # par_p["T_sky"] = ((par_p["Gp"]/hyp["sigma"]) + par_p["T_amb"]**4)**(1/4)
    
    else :
        Tsk = 0.0552*par_p["T_amb"]**1.5

        par_p["T_sky"] = Tsk
        par_p["Gp"] = hyp["sigma"]*(par_p["T_sky"]**4 - par_p["T_amb"]**4)

def change_N_ail(par,N):
    par["N_ail"] = N

def change_a(par,var,a):
    par["lambd_ail"]=a
    Bi_f0(par,var)
    Bi_f1(par,var)
    Bi_f2(par,var)
    Bi_f3(par,var)

def change_air_layer(par,lambd_air):
    old_air_layer = par["lambd_air"]
    k_air = par["k_air"]

    old_R_T = par["R_inter"]

    old_r_air = old_air_layer/k_air
    new_r_air = lambd_air/k_air

    par["R_inter"] = old_R_T - old_r_air + new_r_air
    par["lambd_air"] = lambd_air
    #print(par["R_inter"])

def change_b_htop(par,b_htop):
    par["b_htop"] = b_htop

    change_u(par,par_p,par["u"])

def change_ins(par,e_new,k_new):

    par["R_2"]=e_new/k_new

def change_N_fins_per_EP(par,N):
    par["N_fins_on_abs"] = (6*N)/par["N_harp"]
    par["D"] = (0.160/N)

