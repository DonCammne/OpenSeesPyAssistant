# Module with simple geometry templates
#   Carmine Schipani, 2021

# Import libraries
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from abc import abstractmethod
from copy import copy, deepcopy
import openseespy.postprocessing.internal_plotting_functions as ipltf
import openseespy.postprocessing.Get_Rendering as opsplt
from OpenSeesPyAssistant.Section import *
from OpenSeesPyAssistant.DataManagement import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *
from OpenSeesPyAssistant.Constants import *
from OpenSeesPyAssistant.Fibers import *
from OpenSeesPyAssistant.Connections import *
from OpenSeesPyAssistant.FunctionalFeatures import *
from OpenSeesPyAssistant.MemberModel import *
from OpenSeesPyAssistant.AnalysisAndPostProcessing import *


def Initialize2DModel(data_dir = "Results"):
    # Clear all
    wipe()

    # Build model (2D - 3 DOF/node)
    model('basic', '-ndm', 2, '-ndf', 3)

    # Main Results Folder
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)


def DefineFrameNodes(n_hor_axis: int, n_vert_axis: int, storey_width, storey_height, half_pz_height = np.array([]),
    origin = [0, 0], first_hor_axis = 1, first_vert_axis = 1, show_plot = True):
    if n_hor_axis < 1: raise NegativeValue()
    if n_vert_axis < 1: raise NegativeValue()
    if storey_width < 0: raise NegativeValue()
    if storey_height < 0: raise NegativeValue()
    if len(origin) != 2: raise WrongDimension()
    if first_hor_axis < 1: raise NegativeValue()
    if first_vert_axis < 1: raise NegativeValue()
    if np.size(half_pz_height) != 0 and np.size(half_pz_height) != n_vert_axis: raise WrongDimension()

    if np.size(half_pz_height) == 0: half_pz_height = np.zeros(n_vert_axis)
    node_array = []
    max_n_x = n_hor_axis + first_hor_axis - 1
    max_n_y = n_vert_axis + first_vert_axis - 1
    for xx in range(n_hor_axis):
        x_axis = xx + first_hor_axis
        for yy in range(n_vert_axis):
            y_axis = yy + first_vert_axis
            node_ID = GridIDConvention(x_axis, y_axis, max_n_x, max_n_y)
            node(node_ID, origin[0]+xx*storey_width, origin[1]+yy*storey_height+half_pz_height[yy])
            if y_axis == 1 and half_pz_height[yy] != 0: print("Warning: the nodes at the base have a panel zone height")
            node_array.append(node_ID)
    
    if show_plot:
        plot_nodes(node_array, "Frame geometry template with only nodes", True)
        plt.grid()

    return node_array


def DefineFrameNodesAndElementsSteelIShape(n_hor_axis: int, n_vert_axis: int, storey_width, storey_height,
    list_col: list, list_beam: list, geo_trans_ID: int, N_G = np.array([]), t_dp = np.array([]), panel_zone = True, show_plot = True):
    if np.size(N_G) == 0: N_G = np.zeros(n_vert_axis-1)
    if np.size(t_dp) == 0: t_dp = np.zeros(n_vert_axis-1)

    if np.size(list_col) != n_vert_axis-1: raise WrongDimension()
    if np.size(list_beam) != n_vert_axis-1: raise WrongDimension()
    if np.size(N_G) != n_vert_axis-1: raise WrongDimension()
    if np.size(t_dp) != n_vert_axis-1: raise WrongDimension()
    if geo_trans_ID < 1: raise NegativeValue()


    
    half_pz_height = np.zeros(n_vert_axis)
    if panel_zone:
        for ii, beam in enumerate(list_beam):
            half_pz_height[ii+1] = beam.d/2
    
    node_array = DefineFrameNodes(n_hor_axis, n_vert_axis, storey_width, storey_height, half_pz_height, [0, 0], 1, 1, False)
    
    beam_column_pzspring = [[], [], []]
    for xx in range(n_hor_axis):
        for yy in range(n_vert_axis):
            if yy != 0:
                node_ID = node_array[xx*n_vert_axis + yy]

                # Panel Zone
                if half_pz_height[yy] == 0:
                    col_j_node_ID = node_ID
                    beam_j_node_ID = node_ID
                else:
                    tmp_pz = PanelZoneSteelIShapeSkiadopoulos2021(node_ID, list_col[yy-1], list_beam[yy-1], geo_trans_ID, t_dp[yy-1])
                    tmp_pz.CreateMember()
                    col_j_node_ID = IDConvention(node_ID, 5, 1)
                    beam_j_node_ID = IDConvention(node_ID, 8, 1)
                    beam_column_pzspring[2].append(deepcopy(tmp_pz))

                # Column
                col_i_node_ID = node_array[xx*n_vert_axis + yy - 1]
                col_mat_i = OffsetNodeIDConvention(col_i_node_ID, "vertical", "i")
                col_mat_j = OffsetNodeIDConvention(col_j_node_ID, "vertical", "j")
                tmp_col = SpringBasedElementModifiedIMKSteelIShape(col_i_node_ID, col_j_node_ID, list_col[yy-1], geo_trans_ID,
                    col_mat_i, col_mat_j, N_G[yy-1])
                tmp_col.CreateMember()
                beam_column_pzspring[1].append(deepcopy(tmp_col))

                if xx != 0:
                    # Beam
                    if half_pz_height[yy] == 0:
                        beam_i_node_ID = node_array[(xx-1)*n_vert_axis + yy]
                    else:
                        beam_i_node_ID = IDConvention(node_array[(xx-1)*n_vert_axis + yy], 2, 1)
                    beam_mat_i = OffsetNodeIDConvention(beam_i_node_ID, "horizontal", "i")
                    beam_mat_j = OffsetNodeIDConvention(beam_j_node_ID, "horizontal", "j")
                    print(beam_i_node_ID)
                    tmp_beam = SpringBasedElementModifiedIMKSteelIShape(beam_i_node_ID, beam_j_node_ID, list_beam[yy-1], geo_trans_ID,
                        beam_mat_i, beam_mat_j, N_G[yy-1])
                    tmp_beam.CreateMember()
                    beam_column_pzspring[0].append(deepcopy(tmp_beam))

    
    if show_plot:
        opsplt.plot_model("nodes", "elements")
        
    return beam_column_pzspring


def SubassemblageNodes():
    pass

def RCSSubassemblage():
    # WIP
    pass

