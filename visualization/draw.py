# Initialization

import math
from datetime import datetime
import openpyxl as opxl
from openpyxl.utils.dataframe import dataframe_to_rows
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import model as ty
import proc as pr
import matplotlib.ticker as mtick
import sklearn.metrics

from IPython.core.display import HTML

import heat_transfer as bht

import fluids as fds
import ht 

import general as gen

import os

import scipy.integrate as integrate
import scipy.optimize as sco

import networkx as nx

import plotly.graph_objects as go

import plot_functions_here as pfun

import pyvista as pv

def draw_pv(componentSpecs):
    # Rectangle dimensions in millimeters (L_pan is now the vertical dimension)
    L_pan = componentSpecs['pv']['L_pan'] * 18500
    w_pan = componentSpecs['pv']['w_pan'] * 9860

    fig, ax = plt.subplots(figsize=(6, 6))  # Adjust figsize to control the aspect ratio
    # Create a rectangle patch with swapped dimensions for correct orientation
    rectangle = Rectangle((0, 0), w_pan, L_pan, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rectangle)

    # Set limits for x and y axes with added margins
    margin = 200  # Increase margin if needed
    ax.set_xlim(-margin, w_pan + margin)
    ax.set_ylim(-1.5 * margin, L_pan + margin)  # Increase bottom margin for L_pan label

    # Draw arrows and labels for dimensions
    # w_pan label (horizontal dimension)
    ax.annotate('', xy=(0, -margin / 2), xytext=(w_pan, -margin / 2), arrowprops=dict(arrowstyle='<->'))
    ax.text(w_pan / 2, -margin, f'w_pan = {w_pan} mm', ha='center', va='top')

    # L_pan label (vertical dimension)
    ax.annotate('', xy=(-margin / 2, 0), xytext=(-margin / 2, L_pan), arrowprops=dict(arrowstyle='<->'))
    ax.text(-margin, L_pan / 2, f'L_pan = {L_pan} mm', ha='right', va='center', rotation=90)

    plt.axis('off')  # Turn off axis
    plt.tight_layout()  # Adjust the padding between and around subplots
    plt.show()