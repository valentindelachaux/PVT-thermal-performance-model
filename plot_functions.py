import math
from datetime import datetime
import openpyxl as opxl
from openpyxl.utils.dataframe import dataframe_to_rows
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import model as ty

import networkx as nx
import plotly.graph_objects as go

def plot_model_tuv(u_list,df_res,condi,popt_mod_list,popt_tuv_list,color_list):

    def lin(x,a,b):
            return a*x+b
    
    fig = go.Figure()

    for i in range(len(u_list)):
        # Add traces

        fig = fig.add_trace(go.Scatter(x=-df_res.loc[df_res["u"]==u_list[i]]['-(T_m - T_a)'], y=df_res.loc[df_res["u"]==u_list[i]]['Qdot / AG'],
                            mode = 'markers',
                            marker=dict(color=color_list[2*i]),
                            name='Model 1D - u = '+str(u_list[i])+' m/s')
                            )
        
        fig = fig = fig.add_trace(go.Scatter(x=condi.loc[condi["u"]==u_list[i]]['T_m - T_a'], y=condi.loc[condi["u"]==u_list[i]]['Qdot / AG'],
                            mode = 'markers',
                            marker=dict(color=color_list[2*i+1]),
                            name='TUV - u = '+str(u_list[i])+' m/s')
                            )


        fig = fig.add_trace(go.Scatter(x=-df_res.loc[df_res["u"]==u_list[i]]['-(T_m - T_a)'], y=lin(-df_res.loc[df_res["u"]==u_list[i]]['-(T_m - T_a)'],*popt_mod_list[i]),
                            mode = 'lines',
                            line=dict(color=color_list[2*i], width=1,
                                dash='dashdot'),
                            name='Linear fit model 1D - u = '+str(u_list[i])+' m/s')
                            )
        
        fig = fig.add_trace(go.Scatter(x=condi.loc[condi["u"]==u_list[i]]['T_m - T_a'], y=lin(condi.loc[condi["u"]==u_list[i]]['T_m - T_a'],*popt_tuv_list[i]),
                            mode = 'lines',
                            line=dict(color=color_list[2*i+1], width=1,
                                dash='dashdot'),
                            name='Linear fit TUV - u = '+str(u_list[i])+' m/s')
                            )
    # Set x-axis title
    fig = fig.update_xaxes(title_text="T_m - T_amb")
    fig.update_yaxes(title_text="Power related to gross (W/m2 coll.)")



    fig = fig.update_layout(
        autosize=False,
        width=1200,
        height=700,
        margin=dict(
            l=0,
            r=0,
            b=50,
            t=50,
            pad=1
        ),
    )

    return fig

def plot_graph(df_gr,par,case):

    df_temp0 = df_gr[["T_PV","T_PV_Base_mean","T_PV_absfin_mean","T_abs_mean","T_Base_mean","T_absfin_mean","T_ins_mean","T_ins_tube_mean","T_ins_absfin_mean","T_tube_mean","T_fluid_mean"]].copy()
    df_tr0 = df_gr[["Q_S","Q_top_conv","Q_top_rad","Q_PV_plate","Q_PV_Base","Q_PV_absfin","Q_absfins_Base","Q_tube_fluid","Q_ins_tube_back_conv","Q_ins_tube_back_rad","Q_ins_absfin_back_conv","Q_ins_absfin_back_rad","Q_tube_back_conv","Q_tube_back_rad","Q_absfin_back"]].copy()
    df_temp = df_gr[["T_PV","T_Base_mean","T_absfin_mean","T_ins_tube_mean","T_ins_absfin_mean","T_tube_mean","T_fluid_mean","T_amb"]].copy()
    df_tr = df_gr[["Q_S","Q_top_conv","Q_top_rad","Q_PV_Base","Q_PV_absfin","Q_absfins_Base","Q_ins_tube_back_conv","Q_ins_tube_back_rad","Q_ins_absfin_back_conv","Q_ins_absfin_back_rad","Q_Base_tube","Q_tube_fluid"]].copy()
    df_tr0["Q_S_calc"] = df_tr0["Q_PV_plate"] + df_tr0["Q_top_conv"] + df_tr0["Q_top_rad"]
    # Multiplie par le nombre de canaux !
    # df_tr = par["N_harp"]*df_tr

    tr_labels = []
    for i in range(len(list(df_tr.keys()))):
        tr_labels.append(round(df_tr[list(df_tr.keys())[i]][case],3))

    T_labels = {}
    for i in range(len(list(df_temp.keys()))):
        T_labels[list(df_temp.keys())[i]] = list(df_temp.keys())[i]+' '+str(round(df_temp[list(df_temp.keys())[i]][case]-273.15,3))+'°C'

    T_labels['T_PV_mean'] = T_labels['T_PV']
    T_labels["Sun"] = "Sun 0 W/m2"

    # Create an empty directed graph
    G = nx.DiGraph()

    # Add the "T_PV_mean" node
    G.add_node("T_PV_mean")

    n0 = "Sun"
    n1 = "T_PV_mean"
    n4 = "T_Base_mean"
    n5 = "T_absfin_mean"
    n6 = "T_tube_mean"
    n7 = "T_fluid_mean"
    n8 = "T_ins_tube_mean"
    n9 = "T_ins_absfin_mean"
    n10 = "T_amb"

    G.add_edge(n0,n1,tr=tr_labels[0])
    G.add_edge(n1,n10,tr=str(tr_labels[1])+'+'+str(tr_labels[2])+'='+str(round(tr_labels[1]+tr_labels[2],3)))
    # G.add_edge(n1,n10,key=2,tr=tr_labels[2])
    G.add_edge(n1,n4,tr=tr_labels[3])
    G.add_edge(n1,n5,tr=tr_labels[4])
    G.add_edge(n5,n4,tr=tr_labels[5])
    G.add_edge(n5,n9,tr='')
    G.add_edge(n9,n10,key=1,tr=str(tr_labels[8])+'+'+str(tr_labels[9])+'='+str(round(tr_labels[8]+tr_labels[9],3)))
    G.add_edge(n4,n6,tr=tr_labels[10])
    G.add_edge(n6,n7,tr=tr_labels[11])
    G.add_edge(n6,n8,tr='')
    G.add_edge(n8,n10,tr=str(tr_labels[6])+'+'+str(tr_labels[7])+'='+str(round(tr_labels[6]+tr_labels[7],3)))

    pos = {'T_PV_mean': np.array([0.63147258, 1.01149048]),
    'T_amb': np.array([-0.43238774,  0.60724458]),
    'T_Base_mean': np.array([0.63147258, 0.5303845 ]),
    'T_absfin_mean': np.array([0.21924969, 0.5303845 ]),
    'T_ins_absfin_mean': np.array([-0.08075031,  0.29501669]),
    'T_tube_mean': np.array([ 0.63147258, -0.10697755]),
    'T_fluid_mean': np.array([ 0.63147258, -0.56201014]),
    'T_ins_tube_mean': np.array([-0.28075031, -0.10697755]),
    'Sun': np.array([0.63147258, 1.3])}
    fig, ax = plt.subplots(figsize=(15,10))

    nx.draw_networkx_nodes(G, pos,ax=ax)

    nodeDic = {n: T_labels[n] for n in G.nodes}
    edgeDic = {e: G.get_edge_data(*e)["tr"] for e in G.edges}  

    nx.draw_networkx_labels(G,pos, labels=nodeDic, ax=ax)
    nx.draw_networkx_edges(G, pos,ax=ax)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edgeDic,ax=ax)
    # nx.draw_networkx_edge_labels(G, pos, font_size=10, label_pos=0.5)
    # nx.draw_networkx_labels(G, pos_save,font_size=6)

    ax.axis('off')
    print(df_gr["T_fluid_in"][case]-273.15)
    plt.show()