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

def simu_condi_mpe(componentSpecs,steadyStateConditions_df,l,h_back,L,hyp):
    
    variables = ['N_test','T_guess','G', 'Gp', 'T_amb', 'u', 'T_abs','T_fluid_in', 'T_fluid_out']
    
    # Dataframe object
    df = pd.DataFrame(columns = variables)

    sigma = scc.sigma

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

    sigma = scc.sigma

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
