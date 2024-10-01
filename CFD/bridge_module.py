## IMPORTS

import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_dir, '..', '..', 'PVT-thermal-performance-model')))
sys.path.append(os.path.abspath(os.path.join(current_dir, '..', '..', 'RD-systems-and-test-benches')))
sys.path.append(os.path.abspath(os.path.join(current_dir, '..', '..', 'parallel-flow-distribution-pressure-loss')))
sys.path.append(os.path.abspath(os.path.join(current_dir, '..', '..', 'parallel-flow-distribution-pressure-loss', 'ansys')))

import time
import math
from datetime import datetime
from io import StringIO

import copy
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from numpy.polynomial.polynomial import Polynomial
import scipy.integrate as integrate
import scipy.optimize as sco
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error
import openpyxl
import networkx as nx
import sklearn.metrics
from IPython.core.display import HTML
import plotly.io as pio
from openpyxl.utils.dataframe import dataframe_to_rows
import openpyxl as opxl
import matplotlib.ticker as mtick
from matplotlib.cm import get_cmap
import time

import model_ht as modht

import utils.data_processing as dp

import shutil

from CoolProp.CoolProp import PropsSI

import jou_gen as jg
import ansys_py_bridge as apb
import ansys.fluent.core as pyfluent
import model as ty
import proc as pr
import plot_functions_here as pfun
import heat_transfer as bht
import fluids as fds
import ht
import general as gen

## Fonctions

# Flatten the dictionary and convert to a DataFrame
def flatten_dict(d, parent_key='', sep='_'):
    items = []
    for k, v in d.items():
        new_key = f'{parent_key}{sep}{k}' if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def flatten_dict_in_df(d):
    # Apply the flatten function to each key in the main dictionary
    flattened_data = {k: flatten_dict(v) for k, v in d.items()}

    # Convert to a DataFrame
    df = pd.DataFrame.from_dict(flattened_data, orient='index')

    return df

def check_folder(fp):
    if not os.path.exists(fp):
        os.makedirs(fp)
    return fp

def create_Inputs_fp(SR_caoMesh_fp, testConditionsCode, method):

    i = 0
    Inputs_fp = os.path.join(SR_caoMesh_fp, f'{testConditionsCode}_{method}_try{i}_Inputs')

    while os.path.exists(Inputs_fp):
        i += 1
        Inputs_fp = os.path.join(SR_caoMesh_fp, f"{testConditionsCode}_{method}_try{i}_Inputs")
    os.makedirs(Inputs_fp)

    return Inputs_fp, i

def import_inputs(geometry_path, hypotheses_path, testConditions_path, Inputs_PyFluent_path) :
     
    panelSpecs = pr.import_geometry(geometry_path)
    hyp = pr.create_dict_from_excel(hypotheses_path,'Main')
    steadyStateConditions_df = pd.read_excel(testConditions_path,header=2)
    Inputs_PyFluent = read_Inputs_PyFluent(Inputs_PyFluent_path)

    return panelSpecs, hyp, steadyStateConditions_df, Inputs_PyFluent

def save_BC(solver, folder_path_case) :

    bc_types = ['velocity_inlet', 'pressure_inlet', 'pressure_outlet', 'wall']
    bc_dict = {}

    for bc_type in bc_types:
        if len(getattr(solver.setup.boundary_conditions, bc_type).keys()) > 0:
            for key in getattr(solver.setup.boundary_conditions, bc_type).keys():
                bc_dict[bc_type] = {key : getattr(solver.setup.boundary_conditions, bc_type)[key].get_state(), **bc_dict.get(bc_type, {})}
        else:
            bc_dict[bc_type] = {}

    bc = {
        'velocity-inlet' : flatten_dict_in_df(bc_dict['velocity_inlet']),
        'pressure-inlet' : flatten_dict_in_df(bc_dict['pressure_inlet']),
        'pressure-outlet' : flatten_dict_in_df(bc_dict['pressure_outlet']),
        'wall' : flatten_dict_in_df(bc_dict['wall']),
        }

    for key in bc.keys():
        dp.write_df_in_sheet(bc[key], os.path.join(folder_path_case, 'boundary_conditions.xlsx'), sheet_name = str(key), with_index = True)

def read_Inputs_PyFluent(Inputs_PyFluent_path):
    Inputs_PyFluent = pd.read_excel(Inputs_PyFluent_path)
    Inputs_PyFluent.drop(columns=['comment'], inplace=True, errors='ignore')
    Inputs_PyFluent.set_index('named_expression', inplace=True)

    return Inputs_PyFluent

def create_folder_paths(folder_path, case, case_name = 'cas_') : 
    folder_path_geom = os.path.join(folder_path, f'{case_name}{case}')
    if not os.path.exists(folder_path_geom):
        os.makedirs(folder_path_geom)
    return folder_path_geom

def create_save_paths(hyp, folder_path_case, big_it):
    hyp['CFD_ht_path'] = folder_path_case

    file_path_result_CFD = os.path.join(hyp['CFD_ht_path'], f'phis_{big_it}.csv') ## On ne peut pas l'appeler autrement pour le moment, il y a un appel dans la fonction simu_on_steadyState
    df_one_fp = os.path.join(hyp['CFD_ht_path'], f'df_one_{big_it}.csv')
    df_one_per_part_fp = os.path.join(hyp['CFD_ht_path'], f'df_one_per_part_{big_it}.csv')
    res_fp = os.path.join(hyp['CFD_ht_path'], f'res_{big_it}.xlsx')
    file_path_result_PyFluent = os.path.join(hyp['CFD_ht_path'], f'PyFluent_{big_it}.csv')

    return file_path_result_CFD, file_path_result_PyFluent, df_one_per_part_fp, df_one_fp, res_fp

def init_mesh(tui, Mesh_fp, S2S_fp, caoMeshName) :
    caoMeshCode = caoMeshName.split('__')[0]

    jg.change_mesh(tui, Mesh_fp, caoMeshName)
    read_or_create_radiation(tui, S2S_fp, caoMeshCode)

def update_operating_conditions(tui, p_op, T_amb, boussinesq = True):
    
    tui.define.operating_conditions.operating_pressure(p_op)
    tui.define.operating_conditions.operating_temperature(T_amb)

    if boussinesq:
        tui.define.operating_conditions.operating_density("no")

    else:
        tui.define.operating_conditions.operating_density("yes", PropsSI('D', 'T', 273.15, 'P', 101325, 'air'))

def update_air_properties(solver, p_op, T_op):

    air_dict = solver.setup.materials.fluid["air"].get_state()

    air_dict['density']['option'] = 'boussinesq'
    air_dict['density']['boussinesq'] = PropsSI('D', 'T', T_op, 'P', p_op, 'air')
    air_dict['specific_heat']['option'] = 'constant'
    air_dict['specific_heat']['constant'] = PropsSI('C', 'T', T_op, 'P', p_op, 'air')
    air_dict['thermal_conductivity']['option'] = 'constant'
    air_dict['thermal_conductivity']['constant'] = PropsSI('conductivity', 'T', T_op, 'P', p_op, 'air')
    air_dict['therm_exp_coeff']['option'] = 'constant'
    air_dict['therm_exp_coeff']['constant'] = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', T_op, 'P', p_op, 'air')
    air_dict['viscosity']['option'] = 'constant'
    air_dict['viscosity']['constant'] = PropsSI('viscosity', 'T', T_op, 'P', p_op, 'air')

    solver.setup.materials.fluid["air"].set_state(air_dict)

def change_air_to_incompressible_ideal_gas(solver):

    air_dict = solver.setup.materials.fluid["air"].get_state()

    air_dict['density']['option'] = 'incompressible-ideal-gas'

    solver.setup.materials.fluid["air"].set_state(air_dict)

# # Change conditions for new hx_flat.yd/yu.x surfaces
def adapt_for_fins(tui):

    field_index_list = [2,2,2,4,4,4]

    for i, surface_name in enumerate(["hx_flat_yd_air.1", "hx_flat_yd_air.2", "hx_flat_yd_air.3", "hx_flat_yu_air.1", "hx_flat_yu_air.2", "hx_flat_yu_air.3"]):
        jg.change_bc_wall(tui, surface_name, "conductive_temperature_field", "\"e_hx+220[kg m s^-3 K^-1]/hint_hx\"", f"\"T_field_{field_index_list[i]}\"")

def init_solver(fp, server_code):
    """Initialize the tui and solver objectfs for PyFluent
    
    Args:
        folder_path_ansys (str): Path to the folder containing the server_info file
        server_code (str): Code of the server generated by Fluent
        
    Returns:
        tui (object): TUI object
        solver (object): Solver object
    """

    solver_path = os.path.join(fp, 'server', f'server_info-{server_code}.txt')
    solver = pyfluent.connect_to_fluent(server_info_file_name=solver_path)
    tui = solver.tui
    return tui, solver

def modify_invariable(tui, Inputs_PyFluent):
    theta = Inputs_PyFluent.loc['theta', 'value']
    jg.change_gravity(tui, theta)
    T_amb = Inputs_PyFluent.loc['T_amb', 'value']
    jg.change_named_expression(tui, 'T_amb', T_amb, 'K')

def import_1Dresults_into_Pyfluent(Inputs_PyFluent, df_one, res, T_fluid_in0, method='bridge') :

    nb_hx = int(Inputs_PyFluent.loc['nb_hx','value'])

    Inputs_PyFluent.loc['T_fluid_in_man', 'value'] = T_fluid_in0

    e_PV = Inputs_PyFluent.loc['e_PV','value'] # m
    L_PV = Inputs_PyFluent.loc['L_PV','value'] # m
    w_PV = Inputs_PyFluent.loc['w_PV','value'] # m

    volume_PV_slice = L_PV * w_PV * e_PV

    Inputs_PyFluent.loc['Heat_Generation_Rate', 'value'] = - df_one['Qdot_PV_sky'].values[0] / 4.75 / volume_PV_slice

    if method == 'bridge':

        for i in range(1, nb_hx + 1) :      
            Inputs_PyFluent.loc[f'T_fluid_in_{i}', 'value'] = res[f'part{i+1}']['df_one']['T_fluid_in'].loc[0]
            Inputs_PyFluent.loc[f'a_f_{i}', 'value'] = res[f'part{i+1}']['slices_df']['a_f'].values[0]
            Inputs_PyFluent.loc[f'b_f_{i}', 'value'] = res[f'part{i+1}']['slices_df']['b_f'].values[0]

        Inputs_PyFluent.loc['T_fluid_in_man', 'value'] = res['part1']['df_one']['T_fluid_mean'].loc[0]
        Inputs_PyFluent.loc['T_fluid_out_man', 'value'] = res['part7']['df_one']['T_fluid_mean'].loc[0]

    elif method == 'uniform':
            
        for i in range(1, nb_hx + 1) :        

            Inputs_PyFluent.loc[f'T_fluid_in_{i}', 'value'] = T_fluid_in0
            Inputs_PyFluent.loc[f'a_f_{i}', 'value'] = 1e-15
            Inputs_PyFluent.loc[f'b_f_{i}', 'value'] = 0.

        Inputs_PyFluent.loc['T_fluid_in_man', 'value'] = T_fluid_in0
        Inputs_PyFluent.loc['T_fluid_out_man', 'value'] = T_fluid_in0
    
def import_Inputs_PyFluent_to_FluentCase(tui, Inputs_PyFluent, method = 'init') :

    if method == 'init':

        for named_expression in Inputs_PyFluent.index :
            value = Inputs_PyFluent.loc[named_expression, 'value']
            unit = Inputs_PyFluent.loc[named_expression, 'unit']
            jg.change_named_expression(tui, named_expression, value, unit)

    elif method == 'update':

        nb_hx = int(Inputs_PyFluent.loc['nb_hx','value'])
        for part in range(1, nb_hx + 1):
            jg.change_named_expression(tui, f'T_fluid_in_{part}', Inputs_PyFluent.loc[f'T_fluid_in_{part}', 'value'], 'K')
            jg.change_named_expression(tui, f'a_f_{part}', Inputs_PyFluent.loc[f'a_f_{part}', 'value'], 'm^-1')
            jg.change_named_expression(tui, f'b_f_{part}', Inputs_PyFluent.loc[f'b_f_{part}', 'value'], 'K m^-1')

        jg.change_named_expression(tui, 'Heat_Generation_Rate', Inputs_PyFluent.loc['Heat_Generation_Rate', 'value'], 'W/m^3')
        jg.change_named_expression(tui, 'T_fluid_in_man', Inputs_PyFluent.loc['T_fluid_in_man', 'value'], 'K')
        jg.change_named_expression(tui, 'T_fluid_out_man', Inputs_PyFluent.loc['T_fluid_out_man', 'value'], 'K')

    else:
        raise ValueError('method must be either init or update')

def replace_keys_in_dict(original_dict, key_mapping):

    return {key_mapping.get(old_key, old_key): value for old_key, value in original_dict.items()}

def write_excel_from_dict(dict_of_dataframes, file_name):
    """
    Writes an Excel file with multiple sheets from a dictionary.
    
    Args:
    - file_name: Name of the Excel file to be created (including .xlsx extension).
    - dict_of_dataframes: A dictionary where keys are sheet names and values are Pandas DataFrames.
    """

    with pd.ExcelWriter(file_name, engine='openpyxl') as writer:
        for sheet_name, df in dict_of_dataframes.items():
            df.to_excel(writer, sheet_name=sheet_name, index=True)

def write_res_to_excel(res, res_fp, df_one_per_part_fp):

    res_for_excel = {}

    for part in res.keys():
        res_for_excel[f'{part}_slices_df'] = res[part]['slices_df']

    write_excel_from_dict(res_for_excel, res_fp)

    df_one_per_part = pd.concat([res[part]['df_one'] for part in res.keys() if part!='main'], axis=0)
    df_one_per_part.index = [f'part{i}' for i in range(1, len(df_one_per_part) + 1)]

    df_one_per_part.to_csv(df_one_per_part_fp, sep=';')

def save_simu_1D(df_one, res, df_one_fp, res_fp, df_one_per_part_fp) :

    df_one.to_csv(df_one_fp, sep = ";")
    write_res_to_excel(res, res_fp, df_one_per_part_fp)

        # keys = list(res['part1']['df_one'].keys())
        # values_dict = {key: [] for key in keys}

        # for i in range(1, nb_hx+3):
        #     part_key = f'part{i}'
        #     if part_key in res:
        #         for key in keys:
        #             values_dict[key].append(res[part_key]['df_one'][key].values[0])

        # result_1D_df_one = pd.DataFrame(values_dict)

        # result_1D_df_one.index = [f'part{i}' for i in range(1, nb_hx + 3)]
        # result_1D_df_one.to_csv(df_one_fp, sep=';')

        # keys = list(res['part1']['slices_df'].keys())
        # values_dict = {key: [] for key in keys}

        # for part in range(1, nb_hx+3):
        #     part_key = f'part{part}'
        #     if part_key in res:
        #         for key in keys:
        #             values_dict[key].append(res[part_key]['slices_df'][key].values[0])

        # result_1D_slices_df = pd.DataFrame(values_dict)
        # result_1D_slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
        # result_1D_slices_df

def convert_CFD_results_into_phis(tui, folder_path_case, file_path_result_CFD, file_path_result_PyFluent, Inputs_PyFluent, panelSpecs, hyp, big_it) :
        
        nb_hx = int(Inputs_PyFluent.loc['nb_hx', 'value'])

        jg.write_report(tui,'ht',folder_path_case,f'ht_tot_report_{big_it}')
        jg.write_report(tui,'rad_ht',folder_path_case,f'ht_rad_report_{big_it}')

        Qdot_tot = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_{big_it}.csv'))
        Qdot_rad = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_{big_it}.csv'))

        Qdot = pd.merge(Qdot_tot,Qdot_rad, on='Component')
        Qdot['conv_ht'] = Qdot['ht'] - Qdot['rad_ht']
        Qdot.set_index('Component', inplace=True)

        Qdot.to_csv(os.path.join(folder_path_case,f'all_ht_report_{big_it}.csv'), sep=';')

        if (hyp['fins_CFD'] == 1.) or (hyp['fins_CFD'] == "1"):

            parts_tube_back = [
                ['manifold_yu'],
                ['hx_bend_yu_air', 'hx_bend_yu_pv'],
                ['hx_flat_yu_air', 'hx_flat_yu_air.1', 'hx_flat_yu_air.2', 'hx_flat_yu_air.3'],
                ['hx_bend_mid_air', 'hx_bend_mid_pv'],
                ['hx_flat_yd_air', 'hx_flat_yd_air.1', 'hx_flat_yd_air.2', 'hx_flat_yd_air.3'],
                ['hx_bend_yd_air', 'hx_bend_yd_pv'],
                ['manifold_yd']
            ]

        else:

            parts_tube_back = [
                ['manifold_yu'],
                ['hx_bend_yu_air', 'hx_bend_yu_pv'],
                ['hx_flat_yu_air'],
                ['hx_bend_mid_air', 'hx_bend_mid_pv'],
                ['hx_flat_yd_air'],
                ['hx_bend_yd_air', 'hx_bend_yd_pv'],
                ['manifold_yd']
            ]

        parts_top = [
            [],
            [],
            ['hx_flat_yu_pv-pv_backsheet-cd-cd1-pv-corps'],
            [],
            ['hx_flat_yd_pv-pv_backsheet-cd-cd1-pv-corps'],
            [],
            []
        ]

        parts_tube_fluid = [
            ['manifold_yu'],
            ['hx_bend_yu_air', 'hx_bend_yu_pv'],
            ['hx_flat_yu_air', 'hx_flat_yu_pv-pv_backsheet-cd-cd1-pv-corps'],
            ['hx_bend_mid_air', 'hx_bend_mid_pv'],
            ['hx_flat_yd_air', 'hx_flat_yd_pv-pv_backsheet-cd-cd1-pv-corps'],
            ['hx_bend_yd_air', 'hx_bend_yd_pv'],
            ['manifold_yd']
        ]

        PV = ['pv_front', 'pv_backsheet']

        Qdot_tube_back = []
        Qdot_top = []
        Qdot_top_rad = []
        Qdot_tube_fluid = []
        Qdot_PV_sky = []

        phi = pd.read_csv(r"D:\ANSYS Fluent Projects\pvt_slice_outdoor\Fluent_GMI\tests\phis_template.csv", sep=";")
        phi.set_index('Component', inplace=True)

        Areas_top = [ pr.top_area_tube_contact_PV(panelSpecs[part])/4.75 for part in panelSpecs['decomp'].keys() if part != 'main']
        Areas_back = [ pr.back_area_tube_conv_and_rad(panelSpecs[part])/4.75 for part in panelSpecs['decomp'].keys() if part != 'main']

        for i in range(1, nb_hx + 3):
            if i == 3 or i == 5:

                phi.loc[f'part{i}', 'phi_tube'] = Qdot[Qdot.index.isin(parts_tube_back[i - 1])]['ht'].sum() / Areas_back[i - 1]

                phi.loc[f'part{i}', 'phi_top'] = Qdot[Qdot.index.isin(parts_top[i - 1])]['ht'].sum() / Areas_top[i - 1] if i in [3,5] else 1e-10

                phi.loc[f'part{i}', 'phi_abs'] = 1e-10
        
        phi.to_csv(file_path_result_CFD, sep=';', index=True)

def save_residuals_data_casdat_from_CFD(tui, folder_path_case, case, big_it) :

        jg.write_residuals_file(tui, folder_path_case, f'residuals_cas{case}_it{big_it}')
        jg.compute_mass_flow_rate(tui, 'face-inlet-under-panel', folder_path_case, f'mass_flow_rate_cas{case}_it{big_it}')
        jg.compute_temp_avg(tui, 'face-outlet-under-panel', folder_path_case, f'temp_outlet_cas{case}_it{big_it}')
        jg.compute_surface_temperatures(tui, 'pvt_slice_outdoor_Fluent_GMI_fins', os.path.join(folder_path_case, f'surface_temperatures_it{big_it}'))

        tui.file.write_case_data(os.path.join(r"D:\ANSYS_save", f"casdat_cas{case}_it{big_it}.dat.h5"))
        shutil.move(os.path.join(r"D:\ANSYS_save", f"casdat_cas{case}_it{big_it}.cas.h5"), os.path.join(folder_path_case, f"casdat_cas{case}_it{big_it}.cas.h5"))
        shutil.move(os.path.join(r"D:\ANSYS_save", f"casdat_cas{case}_it{big_it}.dat.h5"), os.path.join(folder_path_case, f"casdat_cas{case}_it{big_it}.dat.h5"))

def copy_inputs(Inputs, fp):

    info_names = ['panelSpecs', 'Model_hypotheses', 'steadyStateConditions_dict', 'Inputs_PyFluent']

    for i, info in enumerate(Inputs):
        if i >= 2:
            with open(os.path.join(fp, f'{info_names[i-2]}.pkl'), 'wb') as f:
                pickle.dump(info, f)

def read_or_create_radiation(tui, S2S_fp, caoMeshCode):
    target_fp = os.path.join(S2S_fp, f'{caoMeshCode}.s2s.h5')
    if not os.path.exists(target_fp):
        jg.create_radiation(tui, S2S_fp, caoMeshCode)
    else:
        jg.read_radiation(tui, S2S_fp, caoMeshCode)

p_op = 101325

def simu_bridge_cases(tui, solver, folder_path, Inputs, nb_it = 50, nb_big_it = 5, T_out_diff_convergence_control = 0.05, method = 'bridge'):

    caoMeshCode, testConditionsCode, panelSpecs, hyp, steadyStateConditions_dict, Inputs_PyFluent = Inputs

    SR_caoMesh_fp = check_folder(os.path.join(folder_path, 'SimulationResults', caoMeshCode))

    nb_case = len(steadyStateConditions_dict)
    nb_hx = int(Inputs_PyFluent.loc['nb_hx','value'])

    Inputs_fp, i_try = create_Inputs_fp(SR_caoMesh_fp, testConditionsCode,method)

    for case in range(nb_case):

        panelSpecs_case = copy.deepcopy(panelSpecs)
        hyp_case = copy.deepcopy(hyp)
        steadyStateConditions_dict_case = copy.deepcopy(steadyStateConditions_dict)
        Inputs_PyFluent_case = copy.deepcopy(Inputs_PyFluent)

        folder_path_case = check_folder( os.path.join(SR_caoMesh_fp, f"{testConditionsCode}_{method}_try{i_try}_case{case}") )
        jg.change_report_file_path(tui, "report-ht-hx-rfile", os.path.join(folder_path_case, f'report_case{case}'))

        if case == 0:
            copy_inputs(Inputs, Inputs_fp)

        log_df = pd.DataFrame(columns=['iteration', 'time', 'log', 'errors'])

        now = datetime.now().strftime("%Y-%m-%d %H-%M")
        with open(os.path.join(folder_path_case, f'{now}.txt'), 'w') as f:
            pass

        stepConditions = steadyStateConditions_dict_case[case].copy()

        big_it = hyp_case['big_it']
        T_amb = stepConditions['T_amb']
        T_fluid_in0 = stepConditions['T_fluid_in0']
        
        file_path_result_CFD, file_path_result_PyFluent, df_one_per_part_fp, df_one_fp, res_fp = create_save_paths(hyp_case, folder_path_case, big_it)

        df_one, res = ty.simu_one_steady_state_all_he(panelSpecs_case, stepConditions, hyp_case)
        T_out_prev = 0.
        T_out_now = df_one['T_fluid_out'].loc[0]

        T_m = df_one['T_fluid_mean'].loc[0]

        update_operating_conditions(tui, p_op, T_amb)
        update_air_properties(solver, p_op, (T_amb + T_m)/2)

        hyp_case['method_h_top_g_exchanger'] = 'CFD'
        hyp_case['method_h_back_abs'] = 'CFD'
        hyp_case['method_h_back_tube'] = 'CFD'
        hyp_case['method_h_rad_back_tube'] = 'CFD'

        import_1Dresults_into_Pyfluent(Inputs_PyFluent_case, df_one, res, T_fluid_in0, method)
        import_Inputs_PyFluent_to_FluentCase(tui, Inputs_PyFluent_case, 'init')
        Inputs_PyFluent_case.to_csv(file_path_result_PyFluent, sep=';')

        # L_abs = panelSpecs_case['main']['L_abs']
        # c = bht.speed_natural_convection(T_fluid_in0, T_amb, theta,L_abs)

        if stepConditions['u'] >= 0.1:
            bc_to_init = "cd2_wind"
            c = stepConditions['u']
            cy =   c/np.sqrt(2)
            cz = - c/np.sqrt(2)

            cd2_wind = solver.setup.boundary_conditions.velocity_inlet['cd2_wind'].get_state()
            cd2_wind['vmag']['constant'] = c
            solver.setup.boundary_conditions.velocity_inlet['cd2_wind'].set_state(cd2_wind)
            jg.standard_initialization(tui, bc_to_init, 0, 0, cy, cz)

        else:
            bc_to_init = "cd_fc"
            # c = bht.speed_natural_convection(stepConditions['T_fluid_in0'], stepConditions['T_amb'], math.radians(45), Inputs_PyFluent_case.loc['L_PV', 'value'])
            c = 0.1
            cy = - c/np.sqrt(2)
            cz =   c/np.sqrt(2)
            jg.standard_initialization(tui, bc_to_init, 0, 0, 0, c)

        save_BC(solver, folder_path_case)

        while big_it < nb_big_it :

            time_start = time.time()

            # Save 1D
            save_simu_1D(df_one, res, df_one_fp, res_fp, df_one_per_part_fp)

            # CFD
            solver.solution.run_calculation.iterate(number_of_iterations=nb_it)
            convert_CFD_results_into_phis(tui, folder_path_case, file_path_result_CFD, file_path_result_PyFluent, Inputs_PyFluent_case, panelSpecs_case, hyp_case, big_it) # on remplit les phis

            save_residuals_data_casdat_from_CFD(tui, folder_path_case, case, big_it)

            try:

                # 1D
                df_one, res = ty.simu_one_steady_state_all_he(panelSpecs_case, stepConditions, hyp_case)

                T_out_prev = T_out_now
                T_out_now = df_one['T_fluid_out'].loc[0]

                T_m = df_one['T_fluid_mean'].loc[0]
                update_air_properties(solver, p_op, (T_amb + T_m)/2)

                # On alimente la CFD 
                if big_it < nb_big_it - 1:
                    import_1Dresults_into_Pyfluent(Inputs_PyFluent_case, df_one, res, T_fluid_in0, method)
                    import_Inputs_PyFluent_to_FluentCase(tui, Inputs_PyFluent_case, 'update')
                    Inputs_PyFluent_case.to_csv(file_path_result_PyFluent, sep=';')

                big_it = big_it + 1
                hyp_case['big_it'] = big_it
                file_path_result_CFD, file_path_result_PyFluent, df_one_per_part_fp, df_one_fp, res_fp = create_save_paths(hyp_case, hyp_case['CFD_ht_path'], hyp_case['big_it'])

            except Exception as e:
                print(e)
                time_end = time.time()
                log_df.loc[big_it] = {'iteration': big_it, 'time': time_end - time_start, 'log' : f'1D model error at big it : {big_it}', 'errors': e}
                break

            time_end = time.time()
            save_simu_1D(df_one, res, df_one_fp, res_fp, df_one_per_part_fp)

            if ((big_it-1) >= 2) and (abs(T_out_now-T_out_prev) < T_out_diff_convergence_control): 
                log_df.loc[big_it] = {'iteration': big_it-1, 'time': time_end - time_start, 'log' : f'temperature profile converged at big it : {big_it}', 'errors': 'None'}
                break
            else:
                log_df.loc[big_it] = {'iteration': big_it-1, 'time': time_end - time_start, 'log' : 'continue' if big_it < nb_big_it else 'big it limit reached', 'errors': 'None'}

        log_df.to_csv(os.path.join(folder_path_case,f'log.csv'), sep=';')