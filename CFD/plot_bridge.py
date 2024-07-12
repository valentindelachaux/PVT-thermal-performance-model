## IMPORTS
import os
import sys
import time
import math
from datetime import datetime
from io import StringIO

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from scipy.stats import linregress
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

sys.path.append(r'D:\seagu_OneDrive\Documents\GitHub\parallel-flow-distribution-pressure-loss\ansys')
sys.path.append(r'D:\seagu_OneDrive\Documents\GitHub\PVT-thermal-performance-model')
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


## FUNCTIONS
def get_data(plot_hyp, panelSpecs, hyp, stepConditions) : 
    method = plot_hyp['method']
    nb_it = plot_hyp['nb_it']
    folder_path	= plot_hyp['folder_path']
    folder_mesh_base	= plot_hyp['folder_mesh']
    folder_name_base = plot_hyp['folder_name']
    new_save = plot_hyp['new_save']

    if new_save == True :
        if method == 'case' :
            no_case = plot_hyp['no_case']
            no_mesh = plot_hyp['no_mesh']

            folder_mesh = folder_mesh_base+f'{no_mesh+1}'
            folder_path_mesh = os.path.join(folder_path, folder_mesh)

            ht_tot_list = []
            ht_rad_list = []
            ht_conv_list = []
            CFD_list = []
            df_one_list = []
            slices_df_list = []
            PyFluent_list = []

            folder_name = folder_name_base+f'{no_case}'
            folder_path_case = os.path.join(folder_path_mesh, folder_name)

            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test')

            for iteration in range(nb_it) :
                file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
                file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
                file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
                file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

                df_PyFluent = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
                df_PyFluent['iteration'] = iteration
                PyFluent_list.append(df_PyFluent)

                nb_hx = int(apb.get_value('nb_hx', 'named_expression', df_PyFluent))

                ht_tot = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_{iteration}.csv'), sep=',')
                ht_tot['iteration'] = iteration
                ht_tot_list.append(ht_tot)

                ht_rad = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_{iteration}.csv'), sep=',')
                ht_rad['iteration'] = iteration
                ht_rad_list.append(ht_rad)
                
                ht_conv = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_{iteration}.csv'), sep=';')
                ht_conv['iteration'] = iteration
                ht_conv_list.append(ht_conv)

                df_CFD = pd.read_csv(file_path_result_CFD, sep=';')
                df_CFD['iteration'] = iteration
                df_CFD.index = [f'part{i}' for i in range(1, nb_hx+3)]
                CFD_list.append(df_CFD)

                df_one = pd.read_csv(file_path_df_one, sep=';')
                df_one['iteration'] = iteration
                df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
                df_one_list.append(df_one)

                slices_df = pd.read_csv(file_path_slices_df, sep=';')
                slices_df['iteration'] = iteration
                slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
                slices_df_list.append(slices_df)

            file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration+1}.csv'
            file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration+1}.csv'

            df_one = pd.read_csv(file_path_df_one, sep=';')
            df_one['iteration'] = iteration+1
            df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
            df_one_list.append(df_one)

            slices_df = pd.read_csv(file_path_slices_df, sep=';')
            slices_df['iteration'] = iteration+1
            slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
            slices_df_list.append(slices_df)
            
            return(ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list)

        elif method == 'mesh' :
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

            ht_tot_mesh_case_list = []
            ht_rad_mesh_case_list = []
            ht_conv_mesh_case_list = []
            CFD_mesh_case_list = []
            df_one_mesh_case_list = []
            slices_df_mesh_case_list = []
            PyFluent_mesh_case_list = []

            for mesh in range(nb_mesh) :

                folder_mesh = folder_mesh_base+f'{mesh+1}'
                folder_path_mesh = os.path.join(folder_path, folder_mesh)

                ht_tot_case_list = []
                ht_rad_case_list = []
                ht_conv_case_list = []
                CFD_case_list = []
                df_one_case_list = []
                slices_df_case_list = []
                PyFluent_case_list = []

                for case in range(nb_cases) :

                    ht_tot_list = []
                    ht_rad_list = []
                    ht_conv_list = []
                    CFD_list = []
                    df_one_list = []
                    slices_df_list = []
                    PyFluent_list = []

                    folder_name = folder_name_base+f'{case}'
                    folder_path_case = os.path.join(folder_path_mesh, folder_name)

                    hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test')

                    for iteration in range(nb_it) :
                        file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
                        file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
                        file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
                        file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

                        df_PyFluent = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
                        df_PyFluent['iteration'] = iteration
                        PyFluent_list.append(df_PyFluent)

                        nb_hx = int(apb.get_value('nb_hx', 'named_expression', df_PyFluent))

                        ht_tot = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_{iteration}.csv'), sep=',')
                        ht_tot['iteration'] = iteration
                        ht_tot_list.append(ht_tot)

                        ht_rad = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_{iteration}.csv'), sep=',')
                        ht_rad['iteration'] = iteration
                        ht_rad_list.append(ht_rad)
                        
                        ht_conv = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_{iteration}.csv'), sep=';')
                        ht_conv['iteration'] = iteration
                        ht_conv_list.append(ht_conv)
                        df_CFD = pd.read_csv(file_path_result_CFD, sep=';')
                        df_CFD['iteration'] = iteration
                        df_CFD.index = [f'part{i}' for i in range(1, nb_hx+3)]
                        CFD_list.append(df_CFD)

                        df_one = pd.read_csv(file_path_df_one, sep=';')
                        df_one['iteration'] = iteration
                        df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
                        df_one_list.append(df_one)

                        slices_df = pd.read_csv(file_path_slices_df, sep=';')
                        slices_df['iteration'] = iteration
                        slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
                        slices_df_list.append(slices_df)

                    file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration+1}.csv'
                    file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration+1}.csv'

                    df_one = pd.read_csv(file_path_df_one, sep=';')
                    df_one['iteration'] = iteration+1
                    df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
                    df_one_list.append(df_one)

                    slices_df = pd.read_csv(file_path_slices_df, sep=';')
                    slices_df['iteration'] = iteration+1
                    slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
                    slices_df_list.append(slices_df)

                    ht_tot_case_list.append(ht_tot_list)
                    ht_rad_case_list.append(ht_rad_list)
                    ht_conv_case_list.append(ht_conv_list)
                    CFD_case_list.append(CFD_list)
                    df_one_case_list.append(df_one_list)
                    slices_df_case_list.append(slices_df_list)
                    PyFluent_case_list.append(PyFluent_list)

                ht_tot_mesh_case_list.append(ht_tot_case_list)
                ht_rad_mesh_case_list.append(ht_rad_case_list)
                ht_conv_mesh_case_list.append(ht_conv_case_list)
                CFD_mesh_case_list.append(CFD_case_list)
                df_one_mesh_case_list.append(df_one_case_list)
                slices_df_mesh_case_list.append(slices_df_case_list)
                PyFluent_mesh_case_list.append(PyFluent_case_list)

            return ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list

        elif method == 'ref' :

            ht_tot_AR_list = []
            ht_rad_AR_list = []
            ht_conv_AR_list = []
            CFD_AR_list = []
            df_one_AR_list = []
            slices_df_AR_list = []
            PyFluent_AR_list = []

            folder_name = folder_name_base+f'ref'
            folder_path_case = os.path.join(folder_path, folder_name)

            # Définir le chemin du dossier pour les fichiers de résultats
            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test')

            # for iteration in range(1, limit_big_it+1) :
            for iteration in range(nb_it) :
                file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
                file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
                file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
                file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

                df_PyFluent = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
                df_PyFluent['iteration'] = iteration
                PyFluent_AR_list.append(df_PyFluent)

                nb_hx = int(apb.get_value('nb_hx', 'named_expression', df_PyFluent))

                ht_tot = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_{iteration}.csv'), sep=',')
                ht_tot['iteration'] = iteration
                ht_tot_AR_list.append(ht_tot)

                ht_rad = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_{iteration}.csv'), sep=',')
                ht_rad['iteration'] = iteration
                ht_rad_AR_list.append(ht_rad)
                
                ht_conv = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_{iteration}.csv'), sep=';')
                ht_conv['iteration'] = iteration
                ht_conv_AR_list.append(ht_conv)

                df_CFD = pd.read_csv(file_path_result_CFD, sep=';')
                df_CFD['iteration'] = iteration
                df_CFD.index = [f'part{i}' for i in range(1, nb_hx+3)]
                CFD_AR_list.append(df_CFD)

                df_one = pd.read_csv(file_path_df_one, sep=';')
                df_one['iteration'] = iteration
                df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
                df_one_AR_list.append(df_one)

                slices_df = pd.read_csv(file_path_slices_df, sep=';')
                slices_df['iteration'] = iteration
                slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
                slices_df_AR_list.append(slices_df)


            file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration+1}.csv'
            file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration+1}.csv'

            df_one = pd.read_csv(file_path_df_one, sep=';')
            df_one['iteration'] = iteration+1
            df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
            df_one_AR_list.append(df_one)

            slices_df = pd.read_csv(file_path_slices_df, sep=';')
            slices_df['iteration'] = iteration+1
            slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
            slices_df_AR_list.append(slices_df)


            iteration = 0
            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test_CFD_uniform')

            file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
            file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
            file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
            file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

            CFD_uniform = pd.read_csv(file_path_result_CFD, sep=';')
            CFD_uniform['iteration'] = iteration
            CFD_uniform.index = [f'part{i}' for i in range(1, nb_hx+3)]


            df_one_uniform = pd.read_csv(file_path_df_one, sep=';')
            df_one_uniform['iteration'] = iteration
            df_one_uniform.index = [f'part{i}' for i in range(1, nb_hx+3)]

            slices_df_uniform = pd.read_csv(file_path_slices_df, sep=';')
            slices_df_uniform['iteration'] = iteration
            slices_df_uniform.index = [f'part{i}' for i in range(1, nb_hx+3)]

            df_PyFluent_uniform = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
            df_PyFluent_uniform['iteration'] = iteration

            ht_tot_uniform = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_uniform_{iteration}.csv'), sep=',')
            ht_tot_uniform['iteration'] = iteration

            ht_rad_uniform = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_uniform_{iteration}.csv'), sep=',')
            ht_rad_uniform['iteration'] = iteration

            ht_conv_uniform = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_uniform_{iteration}.csv'), sep=';')
            ht_conv_uniform['iteration'] = iteration



            iteration = 0
            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test_1D')

            file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
            file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
            file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

            df_one_1D = pd.read_csv(file_path_df_one, sep=';')
            df_one_1D['iteration'] = iteration
            df_one_1D.index = [f'part{i}' for i in range(1, nb_hx+3)]

            slices_df_1D = pd.read_csv(file_path_slices_df, sep=';')
            slices_df_1D['iteration'] = iteration
            slices_df_1D.index = [f'part{i}' for i in range(1, nb_hx+3)]

            df_PyFluent_1D = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
            df_PyFluent_1D['iteration'] = iteration

        return ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D


    else :
        if method == 'case' :
            no_case = plot_hyp['no_case']
            no_mesh = plot_hyp['no_mesh']

            folder_mesh = folder_mesh_base+f'{no_mesh+1}'
            folder_path_mesh = os.path.join(folder_path, folder_mesh)

            ht_tot_list = []
            ht_rad_list = []
            ht_conv_list = []
            CFD_list = []
            df_one_list = []
            slices_df_list = []
            PyFluent_list = []

            folder_name = folder_name_base+f'{no_case}'
            folder_path_case = os.path.join(folder_path_mesh, folder_name)

            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test')

            for iteration in range(nb_it) :
                file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
                file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
                file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
                file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

                df_PyFluent = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
                df_PyFluent['iteration'] = iteration
                PyFluent_list.append(df_PyFluent)

                nb_hx = int(apb.get_value('nb_hx', 'named_expression', df_PyFluent))

                ht_tot = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_{iteration}.csv'), sep=',')
                ht_tot['iteration'] = iteration
                ht_tot_list.append(ht_tot)

                ht_rad = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_{iteration}.csv'), sep=',')
                ht_rad['iteration'] = iteration
                ht_rad_list.append(ht_rad)
                
                ht_conv = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_{iteration}.csv'), sep=';')
                ht_conv['iteration'] = iteration
                ht_conv_list.append(ht_conv)

                df_CFD = pd.read_csv(file_path_result_CFD, sep=';')
                df_CFD['iteration'] = iteration
                df_CFD.index = [f'part{i}' for i in range(1, nb_hx+3)]
                CFD_list.append(df_CFD)

                df_one = pd.read_csv(file_path_df_one, sep=';')
                df_one['iteration'] = iteration
                df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
                df_one_list.append(df_one)

                slices_df = pd.read_csv(file_path_slices_df, sep=';')
                slices_df['iteration'] = iteration
                slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
                slices_df_list.append(slices_df)
            
            return(ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list)

        elif method == 'mesh' :
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

            ht_tot_mesh_case_list = []
            ht_rad_mesh_case_list = []
            ht_conv_mesh_case_list = []
            CFD_mesh_case_list = []
            df_one_mesh_case_list = []
            slices_df_mesh_case_list = []
            PyFluent_mesh_case_list = []

            for mesh in range(nb_mesh) :

                folder_mesh = folder_mesh_base+f'{mesh+1}'
                folder_path_mesh = os.path.join(folder_path, folder_mesh)

                ht_tot_case_list = []
                ht_rad_case_list = []
                ht_conv_case_list = []
                CFD_case_list = []
                df_one_case_list = []
                slices_df_case_list = []
                PyFluent_case_list = []

                for case in range(nb_cases) :

                    ht_tot_list = []
                    ht_rad_list = []
                    ht_conv_list = []
                    CFD_list = []
                    df_one_list = []
                    slices_df_list = []
                    PyFluent_list = []

                    folder_name = folder_name_base+f'{case}'
                    folder_path_case = os.path.join(folder_path_mesh, folder_name)

                    hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test')

                    for iteration in range(nb_it) :
                        file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
                        file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
                        file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
                        file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

                        df_PyFluent = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
                        df_PyFluent['iteration'] = iteration
                        PyFluent_list.append(df_PyFluent)

                        nb_hx = int(apb.get_value('nb_hx', 'named_expression', df_PyFluent))

                        ht_tot = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_{iteration}.csv'), sep=',')
                        ht_tot['iteration'] = iteration
                        ht_tot_list.append(ht_tot)

                        ht_rad = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_{iteration}.csv'), sep=',')
                        ht_rad['iteration'] = iteration
                        ht_rad_list.append(ht_rad)
                        
                        ht_conv = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_{iteration}.csv'), sep=';')
                        ht_conv['iteration'] = iteration
                        ht_conv_list.append(ht_conv)
                        df_CFD = pd.read_csv(file_path_result_CFD, sep=';')
                        df_CFD['iteration'] = iteration
                        df_CFD.index = [f'part{i}' for i in range(1, nb_hx+3)]
                        CFD_list.append(df_CFD)

                        df_one = pd.read_csv(file_path_df_one, sep=';')
                        df_one['iteration'] = iteration
                        df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
                        df_one_list.append(df_one)

                        slices_df = pd.read_csv(file_path_slices_df, sep=';')
                        slices_df['iteration'] = iteration
                        slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
                        slices_df_list.append(slices_df)
                    
                    ht_tot_case_list.append(ht_tot_list)
                    ht_rad_case_list.append(ht_rad_list)
                    ht_conv_case_list.append(ht_conv_list)
                    CFD_case_list.append(CFD_list)
                    df_one_case_list.append(df_one_list)
                    slices_df_case_list.append(slices_df_list)
                    PyFluent_case_list.append(PyFluent_list)

                ht_tot_mesh_case_list.append(ht_tot_case_list)
                ht_rad_mesh_case_list.append(ht_rad_case_list)
                ht_conv_mesh_case_list.append(ht_conv_case_list)
                CFD_mesh_case_list.append(CFD_case_list)
                df_one_mesh_case_list.append(df_one_case_list)
                slices_df_mesh_case_list.append(slices_df_case_list)
                PyFluent_mesh_case_list.append(PyFluent_case_list)

            return ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list

        elif method == 'ref' :

            ht_tot_AR_list = []
            ht_rad_AR_list = []
            ht_conv_AR_list = []
            CFD_AR_list = []
            df_one_AR_list = []
            slices_df_AR_list = []
            PyFluent_AR_list = []

            folder_name = folder_name_base+f'ref'
            folder_path_case = os.path.join(folder_path, folder_name)

            # Définir le chemin du dossier pour les fichiers de résultats
            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test')

            # for iteration in range(1, limit_big_it+1) :
            for iteration in range(nb_it) :
                file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
                file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
                file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
                file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

                df_PyFluent = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
                df_PyFluent['iteration'] = iteration
                PyFluent_AR_list.append(df_PyFluent)

                nb_hx = int(apb.get_value('nb_hx', 'named_expression', df_PyFluent))

                ht_tot = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_{iteration}.csv'), sep=',')
                ht_tot['iteration'] = iteration
                ht_tot_AR_list.append(ht_tot)

                ht_rad = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_{iteration}.csv'), sep=',')
                ht_rad['iteration'] = iteration
                ht_rad_AR_list.append(ht_rad)
                
                ht_conv = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_{iteration}.csv'), sep=';')
                ht_conv['iteration'] = iteration
                ht_conv_AR_list.append(ht_conv)

                df_CFD = pd.read_csv(file_path_result_CFD, sep=';')
                df_CFD['iteration'] = iteration
                df_CFD.index = [f'part{i}' for i in range(1, nb_hx+3)]
                CFD_AR_list.append(df_CFD)

                df_one = pd.read_csv(file_path_df_one, sep=';')
                df_one['iteration'] = iteration
                df_one.index = [f'part{i}' for i in range(1, nb_hx+3)]
                df_one_AR_list.append(df_one)

                slices_df = pd.read_csv(file_path_slices_df, sep=';')
                slices_df['iteration'] = iteration
                slices_df.index = [f'part{i}' for i in range(1, nb_hx+3)]
                slices_df_AR_list.append(slices_df)


            iteration = 0
            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test_CFD_uniform')

            file_path_result_CFD = hyp['CFD_ht_path']+f'_{iteration}.csv' 
            file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
            file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
            file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

            CFD_uniform = pd.read_csv(file_path_result_CFD, sep=';')
            CFD_uniform['iteration'] = iteration
            CFD_uniform.index = [f'part{i}' for i in range(1, nb_hx+3)]


            df_one_uniform = pd.read_csv(file_path_df_one, sep=';')
            df_one_uniform['iteration'] = iteration
            df_one_uniform.index = [f'part{i}' for i in range(1, nb_hx+3)]

            slices_df_uniform = pd.read_csv(file_path_slices_df, sep=';')
            slices_df_uniform['iteration'] = iteration
            slices_df_uniform.index = [f'part{i}' for i in range(1, nb_hx+3)]

            df_PyFluent_uniform = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
            df_PyFluent_uniform['iteration'] = iteration

            ht_tot_uniform = pd.read_csv(os.path.join(folder_path_case,f'ht_tot_report_uniform_{iteration}.csv'), sep=',')
            ht_tot_uniform['iteration'] = iteration

            ht_rad_uniform = pd.read_csv(os.path.join(folder_path_case,f'ht_rad_report_uniform_{iteration}.csv'), sep=',')
            ht_rad_uniform['iteration'] = iteration

            ht_conv_uniform = pd.read_csv(os.path.join(folder_path_case,f'ht_conv_report_uniform_{iteration}.csv'), sep=';')
            ht_conv_uniform['iteration'] = iteration



            iteration = 0
            hyp['CFD_ht_path'] = os.path.join(folder_path_case, 'test_1D')

            file_path_df_one = hyp['CFD_ht_path']+'_df_one' + f'_{iteration}.csv'
            file_path_slices_df = hyp['CFD_ht_path']+'_slices_df' + f'_{iteration}.csv'
            file_path_Inputs_PyFluent = hyp['CFD_ht_path']+'_PyFluent' + f'_{iteration}.csv'

            df_one_1D = pd.read_csv(file_path_df_one, sep=';')
            df_one_1D['iteration'] = iteration
            df_one_1D.index = [f'part{i}' for i in range(1, nb_hx+3)]

            slices_df_1D = pd.read_csv(file_path_slices_df, sep=';')
            slices_df_1D['iteration'] = iteration
            slices_df_1D.index = [f'part{i}' for i in range(1, nb_hx+3)]

            df_PyFluent_1D = pd.read_csv(file_path_Inputs_PyFluent, sep=';')
            df_PyFluent_1D['iteration'] = iteration

            return ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D

        else : 
            raise ValueError('method should be either mesh, case or ref')

def calculate_Qdot(plot_hyp, panelSpecs, hyp, stepConditions, mesh = 0, case = 0, iteration = 0):
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

    method = plot_hyp['method']
    nb_it = plot_hyp['nb_it']

    if method == 'case':
        no_case = plot_hyp['no_case']
        no_mesh = plot_hyp['no_mesh']

        ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
        nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_list[0]))

        ht_tot = ht_tot_list[iteration]
        ht_rad = ht_rad_list[iteration]
        ht_conv = ht_conv_list[iteration]

        Qdot_tube_back = []
        Qdot_top = []
        Qdot_top_rad = []
        Qdot_tube_fluid = []
        Qdot_PV_sky = []

        for i in range(1, nb_hx + 3):
            if i == 3 or i == 5:
                Qdot_tube_back.append(4.75 * ht_conv[ht_conv['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())
                Qdot_tube_fluid.append(-4.75 * ht_tot[ht_tot['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())
                Qdot_top.append(4.75 * ht_conv[ht_conv['Component'].isin(parts_top[i - 1])]['conv_ht'].sum())
                
                Area_part = pr.top_area_tube_contact_PV(panelSpecs[f'part{i}'])
                Area_panel = pr.top_area_tube_contact_PV(panelSpecs['main'])
                
                Qdot_top_rad.append((Area_part / Area_panel) * Qdot_top[-1])
                Qdot_PV_sky.append(
                    4.75 * ht_tot[ht_tot['Component'].isin(parts_top[i - 1])]['ht'].sum() + 
                    (Area_part / Area_panel) * 4.75 * ht_tot[ht_tot['Component'].isin(PV)]['ht'].sum()
                )
            else:
                Qdot_tube_back.append(4.75 * ht_conv[ht_conv['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())
                Qdot_tube_fluid.append(-4.75 * ht_tot[ht_tot['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())
                Qdot_top.append(1e-6)
                Qdot_top_rad.append(1e-6)
                Qdot_PV_sky.append(1e-6)

        return Qdot_tube_fluid, Qdot_top, Qdot_top_rad, Qdot_tube_back, Qdot_PV_sky
    
    elif method == 'mesh' :
        ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
        nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_mesh_case_list[mesh][case][0]))

        ht_tot = ht_tot_mesh_case_list[mesh][case][iteration]
        ht_rad = ht_rad_mesh_case_list[mesh][case][iteration]
        ht_conv = ht_conv_mesh_case_list[mesh][case][iteration]
        
        Qdot_tube_back = []
        Qdot_top = []
        Qdot_top_rad = []
        Qdot_tube_fluid = []
        Qdot_PV_sky = []

        for i in range(1, nb_hx + 3):
            if i == 3 or i == 5:
                Qdot_tube_back.append(4.75 * ht_conv[ht_conv['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())
                Qdot_tube_fluid.append(-4.75 * ht_tot[ht_tot['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())
                Qdot_top.append(4.75 * ht_conv[ht_conv['Component'].isin(parts_top[i - 1])]['conv_ht'].sum())
                
                Area_part = pr.top_area_tube_contact_PV(panelSpecs[f'part{i}'])
                Area_panel = pr.top_area_tube_contact_PV(panelSpecs['main'])
                
                Qdot_top_rad.append((Area_part / Area_panel) * Qdot_top[-1])
                Qdot_PV_sky.append(
                    4.75 * ht_tot[ht_tot['Component'].isin(parts_top[i - 1])]['ht'].sum() + 
                    (Area_part / Area_panel) * 4.75 * ht_tot[ht_tot['Component'].isin(PV)]['ht'].sum()
                )
            else:
                Qdot_tube_back.append(4.75 * ht_conv[ht_conv['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())
                Qdot_tube_fluid.append(-4.75 * ht_tot[ht_tot['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())
                Qdot_top.append(1e-6)
                Qdot_top_rad.append(1e-6)
                Qdot_PV_sky.append(1e-6)

        return Qdot_tube_fluid, Qdot_top, Qdot_top_rad, Qdot_tube_back, Qdot_PV_sky
    elif method == 'ref' :
        ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
        nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_AR_list[0]))

        ht_tot_AR = ht_tot_AR_list[iteration]
        ht_rad_AR = ht_rad_AR_list[iteration]
        ht_conv_AR = ht_conv_AR_list[iteration]

        ht_tot_uniform = ht_tot_uniform[iteration]
        ht_rad_uniform = ht_rad_uniform[iteration]
        ht_conv_uniform = ht_conv_uniform[iteration]

        Qdot_tube_back_AR = []
        Qdot_top_AR = []
        Qdot_top_rad_AR = []
        Qdot_tube_fluid_AR = []
        Qdot_PV_sky_AR = []

        Qdot_tube_back_uniform = []
        Qdot_top_uniform = []
        Qdot_top_rad_uniform = []
        Qdot_tube_fluid_uniform = []
        Qdot_PV_sky_uniform = []

        for i in range(1, nb_hx + 3):
            if i == 3 or i == 5:
                Qdot_tube_back_uniform.append(4.75 * ht_conv_uniform[ht_conv_uniform['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())
                Qdot_tube_back_AR.append(4.75 * ht_conv_AR[ht_conv_AR['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())

                Qdot_tube_fluid_uniform.append(-4.75 * ht_tot_uniform[ht_tot_uniform['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())
                Qdot_tube_fluid_AR.append(-4.75 * ht_tot_AR[ht_tot_AR['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())

                Qdot_top_uniform.append(4.75 * ht_conv_uniform[ht_conv_uniform['Component'].isin(parts_top[i - 1])]['conv_ht'].sum())
                Qdot_top_AR.append(4.75 * ht_conv_AR[ht_conv_AR['Component'].isin(parts_top[i - 1])]['conv_ht'].sum())

                Area_part = pr.top_area_tube_contact_PV(panelSpecs[f'part{i}'])
                Area_panel = pr.top_area_tube_contact_PV(panelSpecs['main'])
                
                Qdot_top_rad_uniform.append((Area_part / Area_panel) * Qdot_top[-1]) # to be checked
                Qdot_top_rad_AR.append((Area_part / Area_panel) * Qdot_top[-1])

                Qdot_PV_sky_uniform.append(
                    4.75 * ht_tot_uniform[ht_tot_uniform['Component'].isin(parts_top[i - 1])]['ht'].sum() + 
                    (Area_part / Area_panel) * 4.75 * ht_tot_uniform[ht_tot_uniform['Component'].isin(PV)]['ht'].sum()
                )
                Qdot_PV_sky_AR.append(
                    4.75 * ht_tot_AR[ht_tot_AR['Component'].isin(parts_top[i - 1])]['conv_ht'].sum() + 
                    (Area_part / Area_panel) * 4.75 * ht_tot_AR[ht_tot_AR['Component'].isin(PV)]['conv_ht'].sum()
                )
            else:
                Qdot_tube_back_uniform.append(4.75 * ht_conv_uniform[ht_conv_uniform['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())
                Qdot_tube_back_AR.append(4.75 * ht_conv_AR[ht_conv_AR['Component'].isin(parts_tube_back[i - 1])]['conv_ht'].sum())
                Qdot_tube_fluid_uniform.append(-4.75 * ht_tot_uniform[ht_tot_uniform['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())
                Qdot_tube_fluid_AR.append(-4.75 * ht_tot_AR[ht_tot_AR['Component'].isin(parts_tube_fluid[i - 1])]['ht'].sum())
                Qdot_top_uniform.append(1e-6)
                Qdot_top_AR.append(1e-6)
                Qdot_top_rad_uniform.append(1e-6)
                Qdot_top_rad_AR.append(1e-6)
                Qdot_PV_sky_uniform.append(1e-6)
                Qdot_PV_sky_AR.append(1e-6)

        return Qdot_tube_fluid_AR, Qdot_top_AR, Qdot_top_rad_AR, Qdot_tube_back_AR, Qdot_PV_sky_AR, Qdot_tube_fluid_uniform, Qdot_top_uniform, Qdot_top_rad_uniform, Qdot_tube_back_uniform, Qdot_PV_sky_uniform
    
    else : 
        raise ValueError('method should be either mesh, case or ref')

def plot_CFD_last_it(Qdot, plot_hyp, panelSpecs, hyp, stepConditions) :
        method = plot_hyp['method']
        nb_it = plot_hyp['nb_it']

        if method == 'case' :
            ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_list[0]))
            no_case = plot_hyp['no_case']
            no_mesh = plot_hyp['no_mesh']

            fig_comparison = go.Figure()

            sim_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            fig_comparison_list = []
            bar_width = 0.8
            bar_positions = np.arange(1, nb_hx + 3)

            values = []

            for part in range(1, nb_hx + 3):
                iteration = nb_it - 1
                Qdot_tube_fluid, Qdot_top, Qdot_top_rad, Qdot_tube_back, Qdot_PV_sky = calculate_Qdot(plot_hyp, panelSpecs, hyp, stepConditions, no_mesh, no_case, iteration)
                if Qdot == 'Qdot_tube_fluid':
                    values.append(Qdot_tube_fluid[part - 1])
                elif Qdot == 'Qdot_top_conv':
                    values.append(Qdot_top[part - 1])
                elif Qdot == 'Qdot_top_rad':
                    values.append(Qdot_top_rad[part - 1])
                elif Qdot == 'Qdot_tube_back':
                    values.append(Qdot_tube_back[part - 1])
                elif Qdot == 'Qdot_PV_sky':
                    values.append(Qdot_PV_sky[part - 1])
                else:
                    raise ValueError('Qdot should be either Qdot_tube_fluid, Qdot_top_conv, Qdot_top_rad, Qdot_tube_back or Qdot_PV_sky')

            fig_comparison.add_trace(go.Bar(
                x=bar_positions,
                y=values,
                name=f'Simulation M{no_mesh + 1}',
                marker_color=sim_colors[no_mesh],
                width=bar_width,
                opacity=0.8
            ))
            fig_comparison.update_layout(
                title=f'Comparaison des simulations pour {Qdot} - Cas {no_case}',
                xaxis_title='Parties',
                yaxis_title= f'{Qdot} [W]',
                xaxis=dict(
                    tickmode='array',
                    tickvals=np.arange(1, nb_hx + 3),
                    ticktext=[f'Part {i}' for i in range(1, nb_hx + 3)]
                ),
                barmode='group',
                legend_title='Simulation'
            )

            fig_comparison.show()

        elif method == 'mesh' :
            ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_mesh_case_list[0][0][0]))
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

            sim_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            fig_comparison_list = []
            bar_width = 0.8 / nb_mesh
            bar_positions = np.arange(1, nb_hx + 3)

            for case in range(nb_cases):

                fig_comparison = go.Figure()

                for mesh in range(nb_mesh):

                    values = []

                    for part in range(1, nb_hx + 3):
                        iteration = nb_it - 1
                        Qdot_tube_fluid, Qdot_top, Qdot_top_rad, Qdot_tube_back, Qdot_PV_sky = calculate_Qdot(plot_hyp, panelSpecs, hyp, stepConditions, mesh, case, iteration)
                        if Qdot == 'Qdot_tube_fluid':
                            values.append(Qdot_tube_fluid[part - 1])
                        elif Qdot == 'Qdot_top_conv':
                            values.append(Qdot_top[part - 1])
                        elif Qdot == 'Qdot_top_rad':
                            values.append(Qdot_top_rad[part - 1])
                        elif Qdot == 'Qdot_tube_back':
                            values.append(Qdot_tube_back[part - 1])
                        elif Qdot == 'Qdot_PV_sky':
                            values.append(Qdot_PV_sky[part - 1])
                        else:
                            raise ValueError('Qdot should be either Qdot_tube_fluid, Qdot_top_conv, Qdot_top_rad, Qdot_tube_back or Qdot_PV_sky')

                    fig_comparison.add_trace(go.Bar(
                        x=bar_positions + (mesh - (nb_mesh - 1) / 2) * bar_width,
                        y=values,
                        name=f'Simulation M{mesh + 1}',
                        marker_color=sim_colors[mesh],
                        width=bar_width,
                        opacity=0.8
                    ))
                fig_comparison.update_layout(
                    title=f'Comparaison des simulations pour Qdot_tube - Cas {case}',
                    xaxis_title='Parties',
                    yaxis_title='Qdot_tube',
                    xaxis=dict(
                        tickmode='array',
                        tickvals=np.arange(1, nb_hx + 3),
                        ticktext=[f'Part {i}' for i in range(1, nb_hx + 3)]
                    ),
                    barmode='group',
                    legend_title='Simulation'
                )

                fig_comparison_list.append(fig_comparison)

            for fig_comparison in fig_comparison_list:
                fig_comparison.show()
        
        elif method == 'ref' :
            ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_AR_list[0]))

            fig_comparison = go.Figure()

            sim_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            fig_comparison_list = []
            bar_width = 0.8 / 2
            bar_positions = np.arange(1, nb_hx + 3)

            values_AR = []
            values_uniform = []
            for part in range(1, nb_hx + 3):
                iteration = nb_it - 1
                Qdot_tube_fluid_AR, Qdot_top_AR, Qdot_top_rad_AR, Qdot_tube_back_AR, Qdot_PV_sky_AR, Qdot_tube_fluid_uniform, Qdot_top_uniform, Qdot_top_rad_uniform, Qdot_tube_back_uniform, Qdot_PV_sky_uniform = calculate_Qdot(plot_hyp, panelSpecs, hyp, stepConditions, mesh = 0, case = 0, iteration = 0)
                if Qdot == 'Qdot_tube_fluid':
                    values_AR.append(Qdot_tube_fluid_AR)
                    values_uniform.append(Qdot_tube_fluid_uniform)

                elif Qdot == 'Qdot_top_conv':
                    values_AR.append(Qdot_top_AR)
                    values_uniform.append(Qdot_top_uniform)

                elif Qdot == 'Qdot_top_rad':
                    values_AR.append(Qdot_top_rad_AR)
                    values_uniform.append(Qdot_top_rad_uniform)

                elif Qdot == 'Qdot_tube_back':
                    values_AR.append(Qdot_tube_back_AR)
                    values_uniform.append(Qdot_tube_back_uniform)

                elif Qdot == 'Qdot_PV_sky':
                    values_AR.append(Qdot_PV_sky_AR)
                    values_uniform.append(Qdot_PV_sky_uniform)
                else :
                    raise ValueError('Qdot should be either Qdot_tube_fluid, Qdot_top_conv, Qdot_top_rad, Qdot_tube_back or Qdot_PV_sky')

            fig_comparison.add_trace(go.Bar(
                x=bar_positions + (0 - (2 - 1) / 2) * bar_width,
                y=values_AR,
                name='Méthode 1D',
                marker_color=sim_colors[no_mesh],
                width=bar_width,
                opacity=0.8
            ))
            fig_comparison.add_trace(go.Bar(
                x=bar_positions + (1 - (2 - 1) / 2) * bar_width,
                y=values_uniform,
                name='Méthode CFD hx uniforme',
                marker_color=sim_colors[no_mesh],
                width=bar_width,
                opacity=0.8
            ))

            fig_comparison.update_layout(
                title=f'Comparaison des simulations pour {Qdot} - Cas référence',
                xaxis_title='Parties',
                yaxis_title= f'{Qdot} [W]',
                xaxis=dict(
                    tickmode='array',
                    tickvals=np.arange(1, nb_hx + 3),
                    ticktext=[f'Part {i}' for i in range(1, nb_hx + 3)]
                ),
                barmode='group',
                legend_title='Simulation'
            )

            fig_comparison.show()
        
        else :
            raise ValueError('method should be either mesh, case or ref')

def plot_1D_last_it(Qdot, plot_hyp, panelSpecs, hyp, stepConditions) :
        method = plot_hyp['method']
        nb_it = plot_hyp['nb_it']

        if method == 'case' :
            ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_list[0]))
            no_case = plot_hyp['no_case']
            no_mesh = plot_hyp['no_mesh']

            fig_comparison_list = []

            sim_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            bar_width = 0.8
            bar_positions = np.arange(1, nb_hx + 4)

            fig_comparison = go.Figure()

            values = []

            for part in range(1, nb_hx + 3):
                iteration = nb_it - 1
                value = df_one_list[iteration][Qdot].loc[f'part{part}']
                values.append(value)

            fig_comparison.add_trace(go.Bar(
                x=bar_positions,
                y=values,
                name=f'Simulation M{no_mesh + 1}',
                marker_color=sim_colors[no_mesh],
                width=bar_width,
                opacity=0.8
            ))

            fig_comparison.update_layout(
                title=f'Comparaison des simulations pour {Qdot} - Cas {no_case}',
                xaxis_title='Parties',
                yaxis_title=f'{Qdot} [W]',
                barmode='group',
                legend_title='Simulation'
            )

            fig_comparison.show()

        elif method == 'mesh' :
            ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_mesh_case_list[0][0][0]))
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

            fig_comparison_list = []

            sim_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            bar_width = 0.8 / nb_mesh
            bar_positions = np.arange(1, nb_hx + 4)

            for case in range(nb_cases):

                fig_comparison = go.Figure()

                for mesh in range(nb_mesh):

                    values = []

                    for part in range(1, nb_hx + 3):
                        iteration = nb_it - 1
                        value = df_one_mesh_case_list[mesh][case][iteration]['Qdot_tube_back_conv'].loc[f'part{part}']
                        values.append(value)

                    fig_comparison.add_trace(go.Bar(
                        x=bar_positions + mesh * bar_width - 0.4,
                        y=values,
                        name=f'Simulation M{mesh + 1}',
                        marker_color=sim_colors[mesh],
                        width=bar_width,
                        opacity=0.8
                    ))

                fig_comparison.update_layout(
                    title=f'Comparaison des simulations pour {Qdot} - Cas {case}',
                    xaxis_title='Parties',
                    yaxis_title=f'{Qdot} [W]',
                    barmode='group',
                    legend_title='Simulation'
                )

                fig_comparison_list.append(fig_comparison)

            for fig_comparison in fig_comparison_list:
                fig_comparison.show()

        elif method == 'ref' :
            ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_AR_list[0]))
            
            fig_comparison_list = []

            sim_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            bar_width = 0.8 / 2
            bar_positions = np.arange(1, nb_hx + 4)

            fig_comparison = go.Figure()

            values_AR = []
            values_1D = []

            for part in range(1, nb_hx + 3):
                iteration = nb_it - 1
                value_AR = df_one_AR_list[iteration][Qdot].loc[f'part{part}']
                values_AR.append(value_AR)
                value_1D = df_one_1D[Qdot].loc[f'part{part}']
                values_1D.append(value_1D)

            fig_comparison.add_trace(go.Bar(
                x=bar_positions + 0 * bar_width - 0.4,
                y=values_AR,
                name='Méthode AR',
                marker_color=sim_colors[0],
                width=bar_width,
                opacity=0.8
            ))

            fig_comparison.add_trace(go.Bar(
                x=bar_positions + 1 * bar_width - 0.4,
                y=values_1D,
                name='Méthode 1D',
                marker_color=sim_colors[1],
                width=bar_width,
                opacity=0.8
            ))


            fig_comparison.update_layout(
                title=f'Comparaison des simulations pour {Qdot} - Cas référence',
                xaxis_title='Parties',
                yaxis_title=f'{Qdot} [W]',
                barmode='group',
                legend_title='Simulation'
            )

            fig_comparison.show()

        else : 
            raise ValueError('method should be either mesh, case or ref')

def plot_big_it(Qdot, plot_hyp, panelSpecs, hyp, stepConditions) :
        method = plot_hyp['method']
        nb_it = plot_hyp['nb_it']

        if method == 'case' :
            ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_list[0]))
            no_case = plot_hyp['no_case']
            no_mesh = plot_hyp['no_mesh']

            part_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

            values_cfd = []
            values_1d = [[] for _ in range(nb_hx + 3)]
            values = [[] for _ in range(nb_hx + 3)]
            iterations_cfd = []
            iterations_1d = [[] for _ in range(nb_hx + 3)]
            iterations = [[] for _ in range(nb_hx + 3)]

            # Créer une figure pour les valeurs de Qdot_top
            fig = go.Figure()

            # Boucle sur chaque itération pour calculer les valeurs et les afficher
            for iteration in range(nb_it):
                Qdot_tube_fluid, Qdot_top, Qdot_top_rad, Qdot_tube_back, Qdot_PV_sky = calculate_Qdot(plot_hyp, panelSpecs, hyp, stepConditions, no_mesh, no_case, iteration)
                if Qdot == 'Qdot_tube_fluid':
                    value_cfd = Qdot_tube_fluid
                elif Qdot == 'Qdot_top_conv':
                    value_cfd = Qdot_top
                elif Qdot == 'Qdot_top_rad':
                    value_cfd = Qdot_top_rad
                elif Qdot == 'Qdot_tube_back':
                    value_cfd = Qdot_tube_back
                elif Qdot == 'Qdot_PV_sky':
                    value_cfd = Qdot_PV_sky
                else:
                    raise ValueError('Qdot not implemented')
                values_cfd.append(value_cfd)
                iterations_cfd.append(iteration + 1)

                # Récupérer les valeurs pour Qdot_top_conv pour 1D
                for part in range(1, nb_hx + 3):
                    value_1d = df_one_list[iteration][Qdot].loc[f'part{part}']
                    values_1d[part-1].append(value_1d)
                    iterations_1d[part-1].append(iteration + 0.5)
                    values[part-1].append(value_1d)
                    values[part-1].append(value_cfd[part-1])
                    iterations[part-1].append(iteration + 0.5)
                    iterations[part-1].append(iteration + 1)

                # Ajouter les valeurs CFD et 1D à la figure
                for part in range(1, nb_hx + 3):
                    # Tracer les valeurs CFD
                    fig.add_trace(go.Scatter(x=iterations_cfd, 
                                                y=[qt[part-1] for qt in values_cfd],
                                                mode='markers',
                                                name=f'Part {part} CFD',
                                                line=dict(width=2, color=part_colors[part-1]),
                                                marker=dict(symbol='cross'),
                                                hovertemplate='CFD Itération: %{x}'+f'<br>{Qdot} :'+ '%{y:.2f} W<br>Part '+f'{part}',))

                    # Tracer les valeurs 1D
                    fig.add_trace(go.Scatter(x=iterations_1d[part-1], 
                                                y=values_1d[part-1],
                                                mode='markers',
                                                name=f'Part {part} 1D',
                                                line=dict(width=2, color=part_colors[part-1], dash='dot'),
                                                marker=dict(symbol='x'),
                                                hovertemplate='1D Itération: %{x}'+f'<br>{Qdot} :'+ '%{y:.2f} W<br>Part '+f'{part}',))
                    
                    fig.add_trace(go.Scatter(x=iterations[part-1], 
                                                y=values[part-1],
                                                mode='lines',
                                                name=f'Part {part}',
                                                line=dict(width=2, color=part_colors[part-1], dash='dot')))

                # Mise en forme du titre et des axes
                fig.update_layout(title=f'{Qdot} - Cas {no_case} - Simu M{no_mesh+1}',
                                    xaxis_title='Itération',
                                    yaxis_title=f'{Qdot} [W]')

            # Affichage de la figure
            fig.show()

        elif method == 'mesh' :
            ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_mesh_case_list[0][0][0]))
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

            part_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

            for case in range(nb_cases):
                
                fig = go.Figure()

                for mesh in range(nb_mesh):
                    values_cfd = []
                    values_1d = [[] for _ in range(nb_hx + 3)]
                    values = [[] for _ in range(nb_hx + 3)]
                    iterations_cfd = []
                    iterations_1d = [[] for _ in range(nb_hx + 3)]
                    iterations_combined = [[] for _ in range(nb_hx + 3)]

                    for iteration in range(nb_it):

                        Qdot_tube_fluid, Qdot_top, Qdot_top_rad, Qdot_tube_back, Qdot_PV_sky = calculate_Qdot(plot_hyp, panelSpecs, hyp, stepConditions, mesh, case, iteration)
                        if Qdot == 'Qdot_tube_fluid':
                            value_cfd = Qdot_tube_fluid
                        elif Qdot == 'Qdot_top_conv':
                            value_cfd = Qdot_top
                        elif Qdot == 'Qdot_top_rad':
                            value_cfd = Qdot_top_rad
                        elif Qdot == 'Qdot_tube_back':
                            value_cfd = Qdot_tube_back
                        elif Qdot == 'Qdot_PV_sky':
                            value_cfd = Qdot_PV_sky
                        else:
                            raise ValueError('Qdot should be either Qdot_tube_fluid, Qdot_top_conv, Qdot_top_rad, Qdot_tube_back or Qdot_PV_sky')
                        values_cfd.append(value_cfd)
                        iterations_cfd.append(iteration + 1)

                        for part in range(1, nb_hx + 3):
                            value_1d = df_one_mesh_case_list[mesh][case][iteration][Qdot].loc[f'part{part}']
                            values_1d[part-1].append(value_1d)
                            iterations_1d[part-1].append(iteration + 0.5)

                            values[part-1].append(value_1d)
                            values[part-1].append(value_cfd[part-1])
                            iterations_combined[part-1].append(iteration + 0.5)
                            iterations_combined[part-1].append(iteration + 1)
                            
                    for part in range(1, nb_hx + 3):
                        fig.add_trace(go.Scatter(x=iterations_cfd, 
                                                    y=[-qt[part-1] for qt in values_cfd],
                                                    mode='markers',
                                                    name=f'Simu M{mesh+1}-Part{part} CFD',
                                                    line=dict(width=2, color=part_colors[part-1]),
                                                    marker=dict(symbol='cross'),
                                                    hovertemplate=f'CFD Simu M{mesh+1}<br>Itération:'+' %{x}'+f'<br>{Qdot} :'+ '%{y}<extra></extra>'))

                        fig.add_trace(go.Scatter(x=iterations_1d[part-1], 
                                                    y=values_1d[part-1],
                                                    mode='markers',
                                                    name=f'Simu M{mesh+1}-Part{part} 1D',
                                                    line=dict(width=2, color=part_colors[part-1], dash='dot'),
                                                    marker=dict(symbol='x'),
                                                    hovertemplate=f'1D Simu M{mesh+1}<br>Itération:' +' %{x}'+f'<br>{Qdot} :'+ '%{y}<extra></extra>'))
                        
                        fig.add_trace(go.Scatter(x=iterations_combined[part-1], 
                                                    y=values[part-1],
                                                    mode='lines',
                                                    name=f'Part {part} Combiné',
                                                    line=dict(width=2, color=part_colors[part-1], dash ='dot')))

                fig.update_layout(title=f'{Qdot} - Cas {case}',
                                    xaxis_title='Itération',
                                    yaxis_title=f'{Qdot} [W]',
                                    legend_title='Partie et Simulation')

                fig.show()

        elif method == 'ref' :
            ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_AR_list[0]))
            part_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

            values_cfd = []
            values_1d = [[] for _ in range(nb_hx + 3)]
            values = [[] for _ in range(nb_hx + 3)]
            iterations_cfd = []
            iterations_1d = [[] for _ in range(nb_hx + 3)]
            iterations = [[] for _ in range(nb_hx + 3)]

            # Créer une figure pour les valeurs de Qdot_top
            fig = go.Figure()

            # Boucle sur chaque itération pour calculer les valeurs et les afficher
            for iteration in range(nb_it):
                iteration = nb_it - 1
                Qdot_tube_fluid_AR, Qdot_top_AR, Qdot_top_rad_AR, Qdot_tube_back_AR, Qdot_PV_sky_AR, Qdot_tube_fluid_uniform, Qdot_top_uniform, Qdot_top_rad_uniform, Qdot_tube_back_uniform, Qdot_PV_sky_uniform = calculate_Qdot(plot_hyp, panelSpecs, hyp, stepConditions, mesh = 0, case = 0, iteration = 0)
                if Qdot == 'Qdot_tube_fluid':
                    value_cfd = Qdot_tube_fluid_AR
                elif Qdot == 'Qdot_top_conv':
                    value_cfd = Qdot_top_AR
                elif Qdot == 'Qdot_top_rad':
                    value_cfd = Qdot_top_rad_AR
                elif Qdot == 'Qdot_tube_back':
                    value_cfd = Qdot_tube_back_AR
                elif Qdot == 'Qdot_PV_sky':
                    value_cfd = Qdot_PV_sky_AR
                else:
                    raise ValueError('Qdot should be either Qdot_tube_fluid, Qdot_top_conv, Qdot_top_rad, Qdot_tube_back or Qdot_PV_sky')
                values_cfd.append(value_cfd)
                iterations_cfd.append(iteration + 1)

                # Récupérer les valeurs pour Qdot_top_conv pour 1D
                for part in range(1, nb_hx + 3):
                    value_1d = df_one_list[iteration][Qdot].loc[f'part{part}']
                    values_1d[part-1].append(value_1d)
                    iterations_1d[part-1].append(iteration + 0.5)
                    values[part-1].append(value_1d)
                    values[part-1].append(value_cfd[part-1])
                    iterations[part-1].append(iteration + 0.5)
                    iterations[part-1].append(iteration + 1)

            # Ajouter les valeurs CFD et 1D à la figure
            for part in range(1, nb_hx + 3):
                # Tracer les valeurs CFD
                fig.add_trace(go.Scatter(x=iterations_cfd, 
                                            y=[qt[part-1] for qt in values_cfd],
                                            mode='markers',
                                            name=f'Part {part} CFD',
                                            line=dict(width=2, color=part_colors[part-1]),
                                            marker=dict(symbol='cross'),
                                            hovertemplate='CFD Itération: %{x}'+f'<br>{Qdot} :'+ '%{y:.2f} W<br>Part %{text}',))

                # Tracer les valeurs 1D
                fig.add_trace(go.Scatter(x=iterations_1d[part-1], 
                                            y=values_1d[part-1],
                                            mode='markers',
                                            name=f'Part {part} 1D',
                                            line=dict(width=2, color=part_colors[part-1], dash='dot'),
                                            marker=dict(symbol='x'),
                                            hovertemplate='1D Itération: %{x}'+f'<br>{Qdot} :'+ '%{y:.2f} W<br>Part %{text}',))
                
                fig.add_trace(go.Scatter(x=iterations[part-1], 
                                            y=values[part-1],
                                            mode='lines',
                                            name=f'Part {part}',
                                            line=dict(width=2, color=part_colors[part-1], dash='dot')))

            # Mise en forme du titre et des axes
            fig.update_layout(title=f'{Qdot} - Cas référence',
                                xaxis_title='Itération',
                                yaxis_title=f'{Qdot} [W]')

            # Affichage de la figure
            fig.show()

        else :
            raise ValueError('method should be either mesh, case or ref')

def T_fluid_fct(y_values, T_fluid_in, a_f, b_f):
        return (T_fluid_in + b_f / a_f) * np.exp(a_f * y_values) - b_f / a_f

def plot_profile_temp(plot_hyp, panelSpecs, hyp, stepConditions) :
        method = plot_hyp['method']
        nb_it = plot_hyp['nb_it']

        if method == 'case' :
            ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_list[0]))
            no_case = plot_hyp['no_case']
            no_mesh = plot_hyp['no_mesh']

            step = apb.get_value('L_1', 'named_expression', PyFluent_list[0]) / 100
            temperature_profiles_it = []

            cmap = px.colors.sequential.Viridis
            fig = go.Figure()

            for iteration in range(nb_it):
                y_values_tot_iter = [0]
                T_fluid_values_tot_iter = [apb.get_value('T_fluid_in_1','named_expression',PyFluent_list[iteration])]
                temperature_profiles = []
                
                for i in range(1, nb_hx + 1):
                    T_fluid_in = apb.get_value(f'T_fluid_in_{i}','named_expression',PyFluent_list[iteration])
                    L = apb.get_value(f'L_{i}', 'named_expression', PyFluent_list[iteration])
                    a_f = apb.get_value(f'a_f_{i}', 'named_expression', PyFluent_list[iteration])
                    b_f = apb.get_value(f'b_f_{i}', 'named_expression', PyFluent_list[iteration])
                    
                    y_values = np.linspace(y_values_tot_iter[-1], y_values_tot_iter[-1] + L, int(L / step))
                    if i == 3:
                        T_fluid_mean_value = apb.T_fluid_mean(T_fluid_in, L, a_f, b_f)
                        T_fluid_values = T_fluid_mean_value * np.ones(len(y_values))
                    else:
                        T_fluid_values = T_fluid_fct(y_values - y_values_tot_iter[-1], T_fluid_in, a_f, b_f)

                    y_values_tot_iter = np.concatenate((y_values_tot_iter, y_values[1:]))
                    T_fluid_values_tot_iter = np.concatenate((T_fluid_values_tot_iter, T_fluid_values[1:]))
                    temperature_profiles.append((y_values, T_fluid_values, f'part{i}'))
                
                temperature_profiles_it.append(temperature_profiles)
                color = cmap[iteration % len(cmap)]
                for y_values, T_fluid_values, part in temperature_profiles_it[iteration]:
                    fig.add_trace(go.Scatter(
                        x=y_values,
                        y=T_fluid_values,
                        mode='lines',
                        name=part + f' (iteration {iteration})',
                        line=dict(color=color),
                        hovertemplate = 'y = %{y:.2f} m <br>Température = %{y:.2f} K'
                    ))
            
            fig.update_layout(
                title='Profils de température pour chaque partie',
                xaxis_title='y (m)',
                yaxis_title='Température (K)',
                legend_title='Partie',
                template='plotly_white'
            )

            fig.show()

        elif method == 'mesh' :
            ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_mesh_case_list[0][0][0]))
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

            fig_list = []
            for case in range(nb_cases) :
                fig = go.Figure()
                for mesh in range(nb_mesh) :
                    
                    step = apb.get_value('L_1', 'named_expression', PyFluent_mesh_case_list[0][0][0]) / 100
                    temperature_profiles_it = []

                    cmap = px.colors.sequential.Viridis

                    for iteration in range(nb_it):
                        y_values_tot_iter = [0]
                        T_fluid_values_tot_iter = [apb.get_value('T_fluid_in_1','named_expression',PyFluent_mesh_case_list[mesh][case][iteration])]
                        temperature_profiles = []
                        
                        for i in range(1, nb_hx + 1):
                            T_fluid_in = apb.get_value(f'T_fluid_in_{i}','named_expression',PyFluent_mesh_case_list[mesh][case][iteration])
                            L = apb.get_value(f'L_{i}', 'named_expression', PyFluent_mesh_case_list[mesh][case][iteration])
                            a_f = apb.get_value(f'a_f_{i}', 'named_expression', PyFluent_mesh_case_list[mesh][case][iteration])
                            b_f = apb.get_value(f'b_f_{i}', 'named_expression', PyFluent_mesh_case_list[mesh][case][iteration])

                            y_values = np.linspace(y_values_tot_iter[-1], y_values_tot_iter[-1] + L, int(L / step))
                            if i == 3:
                                T_fluid_mean_value = apb.T_fluid_mean(T_fluid_in, L, a_f, b_f)
                                T_fluid_values = T_fluid_mean_value * np.ones(len(y_values))
                            else:
                                T_fluid_values = T_fluid_fct(y_values - y_values_tot_iter[-1], T_fluid_in, a_f, b_f)

                            y_values_tot_iter = np.concatenate((y_values_tot_iter, y_values[1:]))
                            T_fluid_values_tot_iter = np.concatenate((T_fluid_values_tot_iter, T_fluid_values[1:]))
                            temperature_profiles.append((y_values, T_fluid_values, f'part{i}'))
                        
                        temperature_profiles_it.append(temperature_profiles)
                        color = cmap[iteration % len(cmap)]
                        for y_values, T_fluid_values, part in temperature_profiles_it[iteration]:
                            fig.add_trace(go.Scatter(
                                x=y_values,
                                y=T_fluid_values,
                                mode='lines',
                                name=part + f' (iteration {iteration}) <br> Simu M{mesh + 1}',
                                line=dict(color=color),
                                hovertemplate = 'y = %{y:.2f} m <br>Température = %{y:.2f} K'
                            ))
                    
                    fig.update_layout(
                        title=f'Profils de température - case {case}',
                        xaxis_title='y (m)',
                        yaxis_title='Température (K)',
                        legend_title='Partie',
                        template='plotly_white'
                    )
                fig_list.append(fig)

            for fig in fig_list :
                fig.show()

        elif method == 'ref' :
            ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_AR_list[0]))

            step = apb.get_value('L_1', 'named_expression', PyFluent_AR_list[0]) / 100
            temperature_profiles_it = []

            cmap = px.colors.sequential.Viridis
            fig = go.Figure()

            for iteration in range(nb_it):
                y_values_tot_iter = [0]
                T_fluid_values_tot_iter = [apb.get_value('T_fluid_in_1','named_expression',PyFluent_AR_list[iteration])]
                temperature_profiles = []
                
                for i in range(1, nb_hx + 1):
                    T_fluid_in = apb.get_value(f'T_fluid_in_{i}','named_expression',PyFluent_AR_list[iteration])
                    L = apb.get_value(f'L_{i}', 'named_expression', PyFluent_AR_list[iteration])
                    a_f = apb.get_value(f'a_f_{i}', 'named_expression', PyFluent_AR_list[iteration])
                    b_f = apb.get_value(f'b_f_{i}', 'named_expression', PyFluent_AR_list[iteration])
                    
                    y_values = np.linspace(y_values_tot_iter[-1], y_values_tot_iter[-1] + L, int(L / step))
                    if i == 3:
                        T_fluid_mean_value = apb.T_fluid_mean(T_fluid_in, L, a_f, b_f)
                        T_fluid_values = T_fluid_mean_value * np.ones(len(y_values))
                    else:
                        T_fluid_values = T_fluid_fct(y_values - y_values_tot_iter[-1], T_fluid_in, a_f, b_f)

                    y_values_tot_iter = np.concatenate((y_values_tot_iter, y_values[1:]))
                    T_fluid_values_tot_iter = np.concatenate((T_fluid_values_tot_iter, T_fluid_values[1:]))
                    temperature_profiles.append((y_values, T_fluid_values, f'part{i}'))
                
                temperature_profiles_it.append(temperature_profiles)
                color = cmap[iteration % len(cmap)]
                for y_values, T_fluid_values, part in temperature_profiles_it[iteration]:
                    fig.add_trace(go.Scatter(
                        x=y_values,
                        y=T_fluid_values,
                        mode='lines',
                        name=part + f' (iteration {iteration})<br> Méthode AR',
                        line=dict(color=color),
                        hovertemplate = 'y = %{y:.2f} m <br>Température = %{y:.2f} K'
                    ))
                

            y_values_tot_iter = [0]
            T_fluid_values_tot_iter = [apb.get_value('T_fluid_in_1','named_expression',df_PyFluent_uniform[0])]
            temperature_profiles = []
            
            for i in range(1, nb_hx + 1):
                T_fluid_in = apb.get_value(f'T_fluid_in_{i}','named_expression',df_PyFluent_uniform[0])
                L = apb.get_value(f'L_{i}', 'named_expression', df_PyFluent_uniform[0])
                a_f = apb.get_value(f'a_f_{i}', 'named_expression', df_PyFluent_uniform[0])
                b_f = apb.get_value(f'b_f_{i}', 'named_expression', df_PyFluent_uniform[0])
                
                y_values = np.linspace(y_values_tot_iter[-1], y_values_tot_iter[-1] + L, int(L / step))
                if i == 3:
                    T_fluid_mean_value = apb.T_fluid_mean(T_fluid_in, L, a_f, b_f)
                    T_fluid_values = T_fluid_mean_value * np.ones(len(y_values))
                else:
                    T_fluid_values = T_fluid_fct(y_values - y_values_tot_iter[-1], T_fluid_in, a_f, b_f)

                y_values_tot_iter = np.concatenate((y_values_tot_iter, y_values[1:]))
                T_fluid_values_tot_iter = np.concatenate((T_fluid_values_tot_iter, T_fluid_values[1:]))
                temperature_profiles.append((y_values, T_fluid_values, f'part{i}'))
            
            temperature_profiles_it.append(temperature_profiles)
            color = cmap[nb_it+1 % len(cmap)]
            for y_values, T_fluid_values, part in temperature_profiles_it[0]:
                fig.add_trace(go.Scatter(
                    x=y_values,
                    y=T_fluid_values,
                    mode='lines',
                    name=part + '<br> Méthode CFD',
                    line=dict(color=color),
                    hovertemplate = 'y = %{y:.2f} m <br>Température = %{y:.2f} K'
                ))


            y_values_tot_iter = [0]
            T_fluid_values_tot_iter = [apb.get_value('T_fluid_in_1','named_expression',df_PyFluent_1D[0])]
            temperature_profiles = []
            
            for i in range(1, nb_hx + 1):
                T_fluid_in = apb.get_value(f'T_fluid_in_{i}','named_expression',df_PyFluent_1D[0])
                L = apb.get_value(f'L_{i}', 'named_expression', df_PyFluent_1D[0])
                a_f = apb.get_value(f'a_f_{i}', 'named_expression', df_PyFluent_1D[0])
                b_f = apb.get_value(f'b_f_{i}', 'named_expression', df_PyFluent_1D[0])
                
                y_values = np.linspace(y_values_tot_iter[-1], y_values_tot_iter[-1] + L, int(L / step))
                if i == 3:
                    T_fluid_mean_value = apb.T_fluid_mean(T_fluid_in, L, a_f, b_f)
                    T_fluid_values = T_fluid_mean_value * np.ones(len(y_values))
                else:
                    T_fluid_values = T_fluid_fct(y_values - y_values_tot_iter[-1], T_fluid_in, a_f, b_f)

                y_values_tot_iter = np.concatenate((y_values_tot_iter, y_values[1:]))
                T_fluid_values_tot_iter = np.concatenate((T_fluid_values_tot_iter, T_fluid_values[1:]))
                temperature_profiles.append((y_values, T_fluid_values, f'part{i}'))
            
            temperature_profiles_it.append(temperature_profiles)
            color = cmap[nb_it+1 % len(cmap)]
            for y_values, T_fluid_values, part in temperature_profiles_it[0]:
                fig.add_trace(go.Scatter(
                    x=y_values,
                    y=T_fluid_values,
                    mode='lines',
                    name=part + '<br> Méthode 1D',
                    line=dict(color=color),
                    hovertemplate = 'y = %{y:.2f} m <br>Température = %{y:.2f} K'
                ))
                
            
            fig.update_layout(
                title='Profils de température - cas de référence',
                xaxis_title='y (m)',
                yaxis_title='Température (K)',
                legend_title='Partie',
                template='plotly_white'
            )

            fig.show()

        else :
            raise ValueError('method should be either mesh, case or ref')

def plot_1D_DeltaT(Qdot, plot_hyp, panelSpecs, hyp, stepConditions) : ## A MODIFIER avec la nouvelle implémentation
        method = plot_hyp['method']
        nb_it = plot_hyp['nb_it']

        if method == 'case' :
            print('Must be mesh method')

        elif method == 'mesh' :
            ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_mesh_case_list[0][0][0]))
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

            part_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

            for mesh in range(nb_mesh):
                fig_top = go.Figure()

                for part in range(1, nb_hx+1): 
                    tube_x = []
                    tube_y = []
                    
                    for case in range(nb_cases):
                        tube_value = df_one_mesh_case_list[mesh][case][-1][Qdot].loc[f'part{part}']
                        T_in = stepConditions[case]['T_fluid_in0']
                        T_out = apb.get_value('T_fluid_in_6', 'named_expression', PyFluent_mesh_case_list[mesh][case][-1])
                        Tmean = (T_in - T_out) / 2
                        delta_temp = Tmean - stepConditions[case]['T_amb']

                        tube_x.append(delta_temp)
                        tube_y.append(tube_value)

                    a, b, r_value, p_value, std_err = linregress(tube_x, tube_y)
                    line_y = np.array(tube_x) * a + b
                    r_value = r_value ** 2
                    
                    fig_top.add_trace(go.Scatter(x=tube_x, y=tube_y,
                                                mode='lines+markers',
                                                name=f'Partie {part}',
                                                line=dict(width=2, color=part_colors[part-1]),
                                                hovertemplate=f'Partie {part}<br>Cas {case + 1}<br>%{{y:.2f}}'))
                    
                    fig_top.add_trace(go.Scatter(x=tube_x, y=line_y,
                                                mode='lines',
                                                name=f'{Qdot} = {a:.2f}*DT + {b:.2f} <br> r²={r_value:.2f}',
                                                line=dict(width=2, color=part_colors[part-1], dash='dash')))

                fig_top.update_layout(title=f'{Qdot} - Simu M{mesh + 1}',
                                    xaxis_title='T_mean-T_amb [K]',
                                    yaxis_title=f'{Qdot} [W]',
                                    legend_title='Partie et Régression Linéaire')

                fig_top.show()

        elif method == 'ref' :
            print('Must be mesh method')

        else :
            raise ValueError('method should be either mesh, case or ref')






def plot_template(Qdot, plot_hyp, panelSpecs, hyp, stepConditions) :
        method = plot_hyp['method']
        nb_it = plot_hyp['nb_it']

        if method == 'case' :
            ht_tot_list, ht_rad_list, ht_conv_list, CFD_list, df_one_list, slices_df_list, PyFluent_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_list[0]))
            no_case = plot_hyp['no_case']
            no_mesh = plot_hyp['no_mesh']

        elif method == 'mesh' :
            ht_tot_mesh_case_list, ht_rad_mesh_case_list, ht_conv_mesh_case_list, CFD_mesh_case_list, df_one_mesh_case_list, slices_df_mesh_case_list, PyFluent_mesh_case_list = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_mesh_case_list[0][0][0]))
            nb_mesh = plot_hyp['nb_mesh']
            nb_cases = plot_hyp['nb_cases']

        elif method == 'ref' :
            ht_tot_AR_list, ht_rad_AR_list, ht_conv_AR_list, CFD_AR_list, df_one_AR_list, slices_df_AR_list, PyFluent_AR_list, ht_tot_uniform, ht_rad_uniform, ht_conv_uniform, CFD_uniform, df_one_uniform, slices_df_uniform, df_PyFluent_uniform, df_one_1D, slices_df_1D, df_PyFluent_1D = get_data(plot_hyp, panelSpecs, hyp, stepConditions)
            nb_hx = int(apb.get_value('nb_hx', 'named_expression', PyFluent_AR_list[0]))

        else :
            raise ValueError('method should be either mesh, case or ref')