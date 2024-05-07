import math
from datetime import datetime
import openpyxl as opxl
from openpyxl.utils.dataframe import dataframe_to_rows
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import model as ty

import scipy.optimize as sco

import os

from IPython.core.display import HTML

def lin(x,a,b):
    return a*x+b

def linear_interpolation_condi(u_list,condi):
    def lin(x,a,b):
        return a*x+b

    popt_list = []
    pcov_list = []

    for i in range(len(u_list)):
        popt,pcov = sco.curve_fit(lin,condi.loc[condi["u"]==u_list[i]]['Tm - T_a'],condi.loc[condi["u"]==u_list[i]]['Qdot / AG'])
        popt_list.append(popt)
        pcov_list.append(pcov)

    return popt_list,pcov_list

def linear_interpolation_df(u_list,df_res):

    popt_list = []
    pcov_list = []

    for i in range(len(u_list)):
        popt,pcov = sco.curve_fit(lin,-df_res.loc[df_res["u"]==u_list[i]]['-(Tm - T_a)'],df_res.loc[df_res["u"]==u_list[i]]['Qdot / AG'])
        popt_list.append(popt)
        pcov_list.append(pcov)

    return popt_list,pcov_list

def AG(par):
    par['AG'] = par['L_pan']*par['w_pan']

def l_B(par):
    par['l_B'] = par['l_c']

def W(par):
    par["W"] = par["L_abs"]/par["N_meander"]

def tube(par):
    H_tube = par['H_tube']
    w_tube = par['w_tube']

    if par['tube_geometry'] == 'circular':
        assert(H_tube == w_tube)
        par['D_tube'] = H_tube
        par['p_int_tube'] = 2*math.pi*(H_tube/2)
    else:
        par['D_tube'] = (2*H_tube*w_tube)/(H_tube+w_tube)
        par['p_int_tube'] = 2*(H_tube+w_tube)

    par['Dext_tube'] = par['D_tube'] + par['lambd_riser_back']

    par['delta'] = (par['W']-par['Dext_tube'])/2

def L_af(par):
    par['L_af'] = (par['W'] - par['l_B'])/2 + 1E-10

def insulated(par):
    if par['lambd_ins']>0:
        par['insulated'] = 1
    else:
        par['insulated'] = 0

def iota(par):
    par['iota'] = 0.

# Does not depend on T_PV
def X_rad(par):
    # Eff_G = par["Eff_G"]
    # G_ref = par["G_ref"]
    # G_T = par["G"]

    # X = 1+Eff_G*(G_T-G_ref)
    par["X_rad"]=1.

# sheet = Comparaisons panneaux

def C_B(par):
    lambd_riser_plate = par["lambd_riser_plate"]
    l_c = par["l_c"]
    k_riser_plate = par["k_riser_plate"]

    par["C_B"] = (l_c*k_riser_plate)/(lambd_riser_plate+1*1E-3)

def R_g(par):
    par["R_g"] = par["lambd_upper_glass"]/par["k_glass"]

def R_top(par):
    par["R_top"] = par['lambd_upper_glass']/par['k_glass'] + par['lambd_upper_EVA']/par['k_EVA']

def R_inter(par):
    par["R_inter"] = par['lambd_si']/par['k_si'] + par['lambd_lower_EVA']/par['k_EVA'] + par['lambd_PVDF']/par['k_PVDF'] + par['lambd_PET']/par['k_PET'] + par['lambd_adh']/par['k_adh'] + par['lambd_lower_glass']/par['k_glass'] + par['lambd_conductive_plate']/par['k_conductive_plate'] + par['lambd_air']/par['k_air'] + par['lambd_abs']/par['k_abs']

def R_2(par):
    par["R_2"] = par['lambd_ins']/par['k_ins'] + 1*1E-12


def create_dict_from_excel(path,sheet):

    df = pd.read_excel(path,sheet_name=sheet)
    df = df[['Label','Value','Unit']]

    for i in range(len(df)):
        if type(df.loc[i,'Label']) != str:
            if math.isnan(df.loc[i,'Label']):
                df.drop(index=i,inplace=True)
                
    df.reset_index(drop=True,inplace=True)
    df.set_index('Label',inplace=True)

    res = {}

    for i in range(len(df)):
        if df.iloc[i]['Unit'] == 'mm':
            res[df.iloc[i].name] = df.iloc[i]['Value']/1000
        else:
            res[df.iloc[i].name] = df.iloc[i]['Value']

    return res

# Retourne le dictionnaire de paramètres correpondant au PVT
def import_geometry(path):

    par = {}
    
    par_decomp = create_dict_from_excel(path,'decomp')
    par_decomp = {k: v for k, v in par_decomp.items() if pd.notnull(v)}

    par_pv = create_dict_from_excel(path,'PV')

    par['decomp'] = par_decomp
    par['pv'] = par_pv

    count_manifold = 0
    par['consider_manifolds'] = 0

    for part in par_decomp.keys():
        par[part] = create_dict_from_excel(path,par_decomp[part])
        par[part].update(par_pv)

        if par_decomp[part] == 'manifold' and count_manifold == 0:
            count_manifold += 1
            par['consider_manifolds'] = 1
            par[part]['is_inlet_man'] = 1
            par[part]['is_outlet_man'] = 0
            par[part]['is_exchanger'] = 0

        elif par_decomp[part] == 'manifold' and count_manifold == 1:
            count_manifold += 1
            par[part]['is_inlet_man'] = 0
            par[part]['is_outlet_man'] = 1
            par[part]['is_exchanger'] = 0

        else:
            par[part]['is_inlet_man'] = 0
            par[part]['is_outlet_man'] = 0
            par[part]['is_exchanger'] = 1
        

    for dic in par.values():
        if type(dic) == dict:
            for i in range(len(dic)):
                if list(dic.keys())[i] in [key for key,value in par_decomp.items()] or type(dic[list(dic.keys())[i]]) != float:
                    pass
                else:
                    if math.isnan(dic[list(dic.keys())[i]]):
                        dic[list(dic.keys())[i]] = 0
                    else:
                        pass
        else:
            pass
    
    for el in [par[key] for key,value in par_decomp.items()]:
        AG(el)
        tube(el)
        l_B(el)
        L_af(el)
        insulated(el)

        iota(el)
        X_rad(el)
        C_B(el)
        R_g(el)
        R_top(el)
        R_inter(el)
        R_2(el)

    par["AG"] = par["main"]["AG"]

    return par

def find_row_from_db(sheet_i,pan):
    for i in range(1,100):
        if sheet_i['A'+str(i)].value == pan:
            return i
        else:
            pass

def create_list_coeff(sheet,ligne):
    parc = ["S","U","V","W","X","Y","Z","AA","AB"]
    res = []
    for i in range(len(parc)):
        s = sheet[parc[i]+str(ligne)].value
        if s == None:  
            s = 0
        else:
            pass    
        res.append(s)
    return res

def disp_html(df):
    display(HTML(df.to_html()))

# coeff = liste eta0hem,a1,a2,a3,a',a5,a6,a7,a8

# Validé 
def comp_power_coeff_Tout(coeff,AG,T_fluid_out,G,Gp,T_fluid_in,T_amb,u):
    return AG*((coeff[0]-coeff[6]*(u-3))*G - (coeff[1]+coeff[3]*(u-3))*((T_fluid_in+T_fluid_out)/2 - T_amb) + coeff[4]*Gp - coeff[7]*(u-3)*Gp - coeff[2]*((T_fluid_in+T_fluid_out)/2 - T_amb)**2 -coeff[8]*((T_fluid_in+T_fluid_out)/2 - T_amb)**4)

###
def comp_power_coeff_Tm(coeff,AG,Tm,G,Gp,T_amb,u):
    return AG*((coeff[0]-coeff[6]*(u-3))*G-(coeff[1]+coeff[3]*(u-3))*((Tm)/2 - T_amb)-coeff[4]*Gp-coeff[7]*(u-3)*Gp)

# Create par

def create_inputs(cond): #cond c'est le nom du fichier
    ## Meteo inputs for tests

    # Simulations with different meteos

    G_list = [966]
    #G_list = [0]
    #G_list = [1000,1000]
    coeff_G_p_list = [0]
    #G_p_list = [0]
    u_list = [0,0.7,1.4,2.7]
    #u_list = [0]
    T_amb_list = [265,270,275,280,285,290,295,300,305]
    #T_amb_list = [280,300]

    #compt = 2

    T_f_in_list = [263,268,273,278,283,288,293,298,303,308,313,318,323]

    T_guess_list = [293]

    return G_list,coeff_G_p_list,u_list,T_amb_list,T_f_in_list,T_guess_list

def create_inputs_from_excel(cond,par,hyp): # cond est le nom du fichier de données Excel
    # TUV
    condi = pd.read_excel(cond,"Data")

    # Enlève la première ligne avec les unités
    condi.drop(index=condi.index[0], 
        axis=0, 
        inplace=True)
    
    dict = {'ta':'T_amb','U':'u','tin':'T_fluid_in','te':'T_fluid_out','tm':'Tm','tm-ta':'Tm - T_a'}

    column_headers = list(condi.columns.values)

    for i in range(len(column_headers)):
        head = column_headers[i]
        if head in list(dict.keys()):
            condi.rename(columns = {head:dict[head]}, inplace = True)

    condi['T_fluid_in'] = condi['T_fluid_in']+273.15
    condi['T_fluid_out'] = condi['T_fluid_out']+273.15
    condi['Tm'] = condi['Tm']+273.15
    condi['T_amb'] = condi['T_amb']+273.15

    condi["u_back"] = hyp["u_back"]

    condi = condi.reset_index()

    condi = condi.drop(columns=["index","Date","UTC"])

    condi = condi.astype('float64')

    condi['-(Tm - T_a)'] = -(condi['Tm'] - condi['T_amb']) # a1

    condi['Qdot / AG'] = condi['Qdot'] / par['AG']   

    condi['-(Tm - T_a)'] = -( condi['Tm'] - condi['T_amb'] ) # a1
    condi['-(Tm - T_a)^2'] = -condi['-(Tm - T_a)']**2 # a2
    condi['-up x (Tm - T_a)'] = ( (condi['u'] - 3) * condi['-(Tm - T_a)'] ) # a3
    # condi['Gp/G'] = df['Gp']/df['G']
    condi['Gp'] = 0. # * df['Gp'] à la place de 0 a4
    condi['-dTm/dt'] = 0. # a5
    condi['-up x G'] = -(condi['u'] - 3) * condi['G'] # a6
    condi['-up x Gp'] = -(condi['u'] - 3) * condi['Gp'] # df['Gp'] à la place de 0 a7
    condi['-(Tm - T_a)^4'] = -condi['-(Tm - T_a)']**4 # a8 

    return condi

def create_out():
    path = os.getcwd()
    pathout = path+'\\Outputs'
    exist = os.path.exists(pathout)

    if exist == False:
        os.makedirs(pathout)

    wbo = opxl.Workbook()
    return pathout,wbo

def find_cell_by_name(wb,nom_variable):
    my_range = wb.defined_names[nom_variable]
    ch = my_range.value
    ch2 = ch.split("!")[1]
    ch3 = ch2.replace("$","")
    return ch3

def create_par():
    par = {}

    par["version"] = 1

    ## Physics constants

    par["sigma"] = 1 # Stefan-Boltzmann constant

    ## PV

    par["eta_nom"] = 1 # nominal efficiency of PV panel

    par["Eff_T"] = 1
    par["T_ref"] = 25 # reference temperature, often 25°C

    par["Eff_G"] = 1
    par["G_ref"] = 1000 # reference irradiance, often 1000 W/m2

    #par["X_rad"] = 1
    par["X_corr"] = 1 

    ## Heat exchanger specs


    par["L_pan"] = 1
    par["w_pan"] = 1
    par["L_abs"] = 1
    par["w_abs"] = 1

    par["W"] = 1 # the width (x-direction) between adjacent fluid tubes
    par["D_tube"] = 1 #
    par["p_int_tube"] = 1 #
    par["p_ext_tube"] = 1 #
    par["w_tube"] = 1
    par["l_c"] = 1
    par["l_B"] = 1
    par["L_af"] = 1  
    par["iota"] = 1
    par["N_meander"] = 1
    par["N_harp"] = 1 # number of identical tubes carrying fluid through the collector
    par["N_harp_actual"] = 1 # actual numer of identical tubes in parallel (harp geometry)
    par["L_tube"] = 1 # the length of the collector along the flow direction = L_tube in Excel

    par["D"] = 1

    par["insulated"] = 1

    ## Additionnal fins = ailettes

    par["ailette"] = 1

    par["geometry"] = 1

    par["fin_0"] = 1
    par["N_f0"] = 1
    par["L_f0"] = 1
    par["delta_f0"] = 1

    par["fin_1"] = 1
    par["N_f1"] = 1
    par["L_f1"] = 1
    par["delta_f1"] = 1
    par["delta_f1_int"] = 1
    par["coeff_f1"] = 1

    par["fin_2"] = 1
    par["N_f2"] = 1
    par["L_f2"] = 1
    par["delta_f2"] = 1

    par["fin_3"] = 1
    par["N_f3"] = 1
    par["L_f3"] = 1
    par["delta_f3"] = 1

    par["Heta"] = 1

    par["N_ail"] = 1

    ## Thermal / radiance

    par["tau_alpha"] = 1 # transmittance-absorptance product for the solar collector
    par["eps"] = 1 # emissivity of the top surface of the collector (PV surface)
    par["eps_he"] = 1

    # Geometry and thermal conductivities

    par["lambd_glass"] = 1
    par["k_glass"] = 1

    par["k_air"] = 1
    par["lambd_air"] = 1

    par["k_abs"] = 1 # thermal conductivity of the plate material
    par["lambd_abs"] = 1 # thickness of the absorber plate

    par["lambd_riser_plate"] = 1
    par["lambd_riser_back"] = 1
    par["k_riser_plate"] = 1
    par["k_riser_back"] = 1

    par["k_ail"] = 1 # thermal conductivity of the fin
    par["lambd_ail"] = 1 # thickness of the fin

    par["k_ins"] = 1
    par["lambd_ins"] = 1

    par["coeff_h_top_free"] = 1
    par["coeff_h_top_forced"] = 1
    par["coeff_h_back"] = 1

    par["h_rad_back"] = 1

    par["h_rad_f"] = 1

    # Ci-dessous les résistances du panneau sont calculées directement dans le fichier Inputs.xlsx

    par["R_inter"] = 1 # instead of 1/h_outer in the document
    par["R_inter"] = 1 # = R_1 = R_inter = resistance to heat transfer from the PV cells to the absorber plate
    par["R_2"] = 1

    par["C_B"] = 1 # the conductance between the absorber plate and the bonded tube

    ## Initialisation d'une météo

    par["G"] = 1 # total solar radiation (beam + diffuse) incident upon the collector surface = POA irradiance
    par["Gp"] = 1 # infra-red 
    par["coeff_G_p"] = 1
    par["T_sky"] = 1 # sky temperature for long-wave radiation calculations
    par["T_amb"] = 1 
    par["T_back"] = 1
    par["u"] = 1 # wind speed

    ## Fluid

    par["T_fluid_in0"] = 1
    par["Cp"] = 1 # specific heat of the fluid flowing through the PV/T collector
    par["mdot"] = 1 # flow rate of fluid through the solar collector

    par["k_fluid"] = 1
    par["rho_fluid"] = 1
    par["mu_fluid"] = 1

    ## Installation

    par["theta"] = 1 # angle of incidence

    par["orientation"] = 1

    ## Type de test

    par["test"] = 1

    # Excel parameters

    list_parameters = []
    *list_parameters, = par

    path = os.getcwd()
    print(path)

    inp = r'\Inputs.xlsm'
    fichier_i = path+inp
    wbi = opxl.load_workbook(fichier_i,data_only=True)
    sheet_i = wbi["Main"]

    # Find parameters in Excel file Inputs.xlsx

    for i in range(len(list_parameters)):
        nom_var = list_parameters[i]
        cell = find_cell_by_name(wbi,nom_var)
        valeur = sheet_i[cell].value
        
        par[nom_var]=valeur

    wbi.close()

    ### Computation of some parameters from inputs

    # Calculate AG
    par["AG"] = par["L_pan"]*par["w_pan"]

    # Calculate delta : demi-intervalle entre deux risers (extérieur à extérieur)
    # utilisé dans gamma_2_int et Q_abs_back

    par["Dext_tube"] = par["D_tube"] + par["lambd_riser_back"]

    par["delta"] = (par["W"]-par["Dext_tube"])/2
    # Calculate X_rad which depends on G
    ty.X_rad(par)

    # Calculate the conductance between the absorber and the fluid through any riser
    ty.C_B(par)

    # Calculate h_fluid
    ty.h_fluid(par)

    return par

def write_stepConditions_from_steadyStateConditions(steadyStateConditions_df,i,hyp):

    stepConditions = {'G':steadyStateConditions_df["G"][i],"T_amb":steadyStateConditions_df["T_amb"][i],"T_back":steadyStateConditions_df["T_amb"][i],"u":steadyStateConditions_df["u"][i], "u_back" : steadyStateConditions_df["u_back"][i], "T_fluid_in0":steadyStateConditions_df["T_fluid_in"][i]}
    ty.change_T_sky(stepConditions,hyp,'TUV')  # calculate Gp and T_sky

    stepConditions['T_back'] = stepConditions['T_amb']
    stepConditions['T_back_rad'] = stepConditions['T_amb']

    stepConditions["mdot"] = steadyStateConditions_df["mdot"][i]

    stepConditions["guess_T_PV"] = (stepConditions["T_amb"]+stepConditions["T_fluid_in0"])/2

    return stepConditions

# Pre-processing and processing functions for parametric studies

def pre_proc(test):
    if test == "lambd_air" or test == "air_layer_TUV":
        return np.linspace(0,0.005,20)
    elif test == "T_guess":
        return np.linspace(280,300,11)
    elif test == "L_f2" or test == "L_f2_TUV":
        return np.linspace(0.001,0.2,30)
    elif test == "coeff_h_top_TUV":
        return [0.9,1,1.05,1.1,1.15,1.2]
    elif test == "coeff_h_back_TUV":
        return np.linspace(1,2,21)
    elif test == "N_riser":
        return [3,6,9,12,15,18,21,24,30,40,50,60]
    elif test == "N_fins_per_EP":
        return [6,7,8,10,11,12,13,14,15,16,18,20,22]
    elif test == "coeff_h_top":
        return np.linspace(0.001,1,20)
    elif test == "b_htop":
        return np.linspace(1,5,6)
    elif test == "parametric_insulation":
        return np.linspace(0,0.1,31)
    elif test == "iota":
        return [0.,0.02,0.05,0.08,0.1,0.2,0.5,1]
    elif test == "D_tube":
        return [0.001,0.002,0.004,0.008,0.01,0.012,0.02,0.05]
    elif test == "k_riser":
        return [1,50,100,150,200,250,300,350,400]
    elif test == "L_f2":
        return np.linspace(0,0.6,10)
    elif test == "absorber":
        return [4,50,100,150,200,250,300,400]
    elif test == "e_abs":
        return np.linspace(0.0001,0.003,50)
    elif test == "a_htop_TUV":
        return [0.5,1,1.5,2,3,4,5,6,7,8,9,10]
    elif test == "N_f1_TUV":
        return [15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,120,150]
    else:
        return []


def proc(par,test,i,test_list):
    if test == "lambd_air" or test == "air_layer_TUV":
        ty.change_air_layer(par,test_list[i])
    elif test == "T_guess":
        T_guess_list = [test_list[i]]
    elif test == "L_f2" or test =="L_f2_TUV":
        par["L_f2"] = test_list[i]
        par["Heta"] = test_list[i]
    elif test == "coeff_h_top_TUV":
        par["coeff_h_top"] = test_list[i]
    elif test == "coeff_h_back_TUV":
        par["coeff_h_back"] = test_list[i]
    elif test == "N_riser":
        if par["geometry"]=="meander":
            par["N_meander"]=test_list[i]
            par["W"]=par["L_abs"]/test_list[i]
        elif par["geometry"] == "harp":
            par["N_harp"]=test_list[i]
            par["N_harp_actual"]=test_list[i]
            par["W"]=par["w_abs"]/test_list[i]

        par["l_i"]=par["W"]
        par["delta"] = (par["W"]-par["Dext_tube"])/2
        par["L_af"]=(par["W"]-par["l_B"])/2
    elif test == "N_fins_per_EP":
        ty.change_N_fins_per_EP(par,test_list[i])
    elif test == "coeff_h_top":
        par["coeff_h_top"] = test_list[i]
    elif test == "b_htop":
        ty.change_b_htop(par,test_list[i])
    elif test == "parametric_insulation":
        ty.change_ins(par,test_list[i])
    elif test == "iota":
        par["iota"] = test_list[i]
    elif test == "D_tube":
        par["D_tube"] = test_list[i]
        par["iota"] = test_list[i]
        ty.h_fluid(par)
    elif test == "k_riser":
        par["k_riser"] = test_list[i]
        ty.C_B(par)
    elif test == "L_f2":
        par["L_f2"] = test_list[i]
    elif test == "absorber":
        par["k_abs"]=test_list[i]
    elif test == "e_abs":
        R_inter = par["R_inter"]
        old_R_abs = par["lambd_abs"]/par["k_abs"]
        par["lambd_abs"]=test_list[i]
        new_R_abs = par["lambd_abs"]/par["k_abs"]
        par["R_inter"] = R_inter - old_R_abs + new_R_abs
    elif test == "a_htop_TUV":
        par["a_htop"] = test_list[i]
    elif test == "N_f1_TUV":
        par["N_f1"] = test_list[i]
        par["N_ail"] = test_list[i]
        par["D"] = (par["w_abs"]-par["N_ail"]*par["lambd_ail"])/(par["N_ail"]-1)
    elif test == "N_f0_TUV":
        par["N_f0"] = test_list[i]
        par["N_ail"] = test_list[i]
        par["D"] = (par["w_abs"]-par["N_ail"]*par["lambd_ail"])/(par["N_ail"]-1)
    else:
        pass

def display_a_i(X):

    huit = len(X[0])-1

    index_coeff = ['eta0,hem','a1','a2','a3','a4','a5','a6','a7','a8']

    for l in range(len(index_coeff)):
        print(index_coeff[l],' : ',round(X[0][l],3))

    print(round(X[0][0] - (X[0][6]*(1.3-3)),3)*100,'%')

    print(round(X[0][1] + X[0][3]*(1.3-3),1))

def A0_A1(X):

    return (round(X[0][0] - (X[0][5]*(1.3-3)),3),round(X[0][1] + X[0][3]*(1.3-3),1))