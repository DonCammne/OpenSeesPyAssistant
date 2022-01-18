"""
Module with geometry templates (nodes and/or elements with associated fibers, material models, etc).
Carmine Schipani, 2021
"""

# Import libraries
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import os
from copy import copy, deepcopy
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
    """
    Function that initialise the project creating the 2D model with 3 DOF per node and set up a directory for the results.

    @param data_dir (str, optional): Directory where the data will be stored.
        The function forces the user to define it just for good practice and consistency between projects.
        Defaults to "Results".
    """
    # Clear all
    wipe()

    # Build model (2D - 3 DOF/node)
    model('basic', '-ndm', 2, '-ndf', 3)

    # Main Results Folder
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)


def DefineFrameNodes(n_hor_axis: int, n_vert_axis: int, storey_width, storey_height, half_pz_height = np.array([]),
    origin = [0, 0], first_hor_axis = 1, first_vert_axis = 1, show_plot = True):
    """
    Function that declares and initialises the grid nodes of a frame. Option to offset the grid node of the panel zones
        with the master node of the panel zone being the grid one (top center one). The function can be used multiple times 
        to create more complex geometries.

    @param n_hor_axis (int): Number of horizontal axis (or piers) for the grid of the frame.
    @param n_vert_axis (int): Number of vertical axis (or floors) for the grid of the frame.
    @param storey_width (float): Width of the bays.
    @param storey_height (float): Height of the storeys.
    @param half_pz_height (np.ndarray, optional): Array of 1 dimension with half the height of the panel zone for each floor.
        The first floor should be 0 (no panel zone in the support). Defaults to np.array([]), e.g. no panel zone.
    @param origin (list, optional): List of two entry with the origin position. Defaults to [0, 0].
    @param first_hor_axis (int, optional): Number of the first pier. Defaults to 1.
    @param first_vert_axis (int, optional): Number of the first floor. Defaults to 1.
    @param show_plot (bool, optional): Option to show the plot of the nodes declared and initialised. Defaults to True.

    @exception NegativeValue: n_hor_axis needs to be a positive integer.
    @exception NegativeValue: n_vert_axis needs to be a positive integer.
    @exception NegativeValue: storey_width needs to be positive.
    @exception NegativeValue: storey_height needs to be positive.
    @exception WrongDimension: origin has a dimension of 2.
    @exception NegativeValue: first_hor_axis needs to be a positive integer.
    @exception NegativeValue: first_vert_axis needs to be a positive integer.
    @exception WrongDimension: size of half_pz_height needs to be equal to n_vert_axis, if different from 0.

    @returns list: List with the nodes declared. 
    """
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
    list_col: list, list_beam: list, geo_trans_ID: int, N_G = np.array([]), t_dp = np.array([]), L_b_col = np.array([]), L_b_beam = np.array([]),
    fix_support = True, show_plot = True, panel_zone = True):
    """
    WIP (Work In Progress). Function that declares and initialises the grid nodes of a frame and the members using steel I shape SpringBasedElements.
    WARNING: Current limit of the geometry: n_hor_axis and n_vert_axis < 10; if exceeded, there are problems with the IDs (ID limit is exceeded, ~2.2e9).
    WARNING: if the section of the columns change, the function does not account for the splacing. Each colum section is defined from floor to floor;
      if there is a change in the column section, it happens right after the panel zone (not realistic but good enough for predesign).
    WIP: Solve ID limit for large building need implementations (for example the use of a different ID convention or the use of the class IDGenerator).

    @param n_hor_axis (int): Number of horizontal axis (or piers) for the grid of the frame.
    @param n_vert_axis (int): Number of vertical axis (or floors) for the grid of the frame.
    @param storey_width (float): Width of the bays.
    @param storey_height (float): Height of the storeys.
    @param list_col (list(SteelIShape)): List with the sections of the columns for every floor.
    @param list_beam (list(SteelIShape)): List with the sections of the beams for every bay.
    @param geo_trans_ID (int): The geometric transformation (for more information, see OpenSeesPy documentation).
    @param N_G (np.ndarray, optional): Array of dimension 1 with the axial load for each column (starting at floor 2). Defaults to np.array([]), e.g. 0.
    @param t_dp (np.ndarray, optional): Array of dimension 1 with the doubler plate thickness for each bay's beam. Defaults to np.array([]), e.g. 0.
    @param L_b_col (np.ndarray, optional): Array of dimension 1 with the maxiaml unbraced lateral buckling length for each column. Defaults to np.array([]), e.g. -1.
    @param L_b_beam (np.ndarray, optional): Array of dimension 1 with the maxiaml unbraced lateral buckling length for each beam. Defaults to np.array([]), e.g. -1.
    @param fix_support (bool, optional): Option to fix the support of the frame. Defaults to True.
    @param show_plot (bool, optional): Option to show the plot of the nodes declared and initialised. Defaults to True.
    @param panel_zone (bool, optional): Option to add the panel zones in the model. Defaults to True.

    @exception WrongDimension: N_G dimension needs to be equal to n_vert_axis-1, if different from 0.
    @exception WrongDimension: t_dp dimension needs to be equal to n_vert_axis-1, if different from 0.
    @exception WrongDimension: L_b_col dimension needs to be equal to n_vert_axis-1, if different from 0.
    @exception WrongDimension: L_b_beam dimension needs to be equal to n_hor_axis-1, if different from 0.
    @exception WrongDimension: list_col dimension needs to be equal to n_vert_axis-1.
    @exception WrongDimension: list_beam dimension needs to be equal to n_vert_axis-1.
    @exception NegativeValue: geo_trans_ID needs to be a positive integer.

    @returns List: List with the element objects in the frame.
    """
    panel_zone = True
    if np.size(N_G) == 0: N_G = np.zeros(n_vert_axis-1)
    if np.size(t_dp) == 0: t_dp = np.zeros(n_vert_axis-1)
    if np.size(L_b_col) == 0: L_b_col = np.ones(n_vert_axis-1) * (-1.0)
    if np.size(L_b_beam) == 0: L_b_beam = np.ones(n_hor_axis-1) * (-1.0)

    if np.size(list_col) != n_vert_axis-1: raise WrongDimension()
    if np.size(list_beam) != n_vert_axis-1: raise WrongDimension()
    if np.size(N_G) != n_vert_axis-1: raise WrongDimension()
    if np.size(t_dp) != n_vert_axis-1: raise WrongDimension()
    if np.size(L_b_col) != n_vert_axis-1: raise WrongDimension()
    if np.size(L_b_beam) != n_hor_axis-1: raise WrongDimension()
    if geo_trans_ID < 1: raise NegativeValue()
    
    half_pz_height = np.zeros(n_vert_axis)
    if panel_zone:
        for ii, beam in enumerate(list_beam):
            half_pz_height[ii+1] = beam.d/2
    
    node_array = DefineFrameNodes(n_hor_axis, n_vert_axis, storey_width, storey_height, half_pz_height, [0, 0], 1, 1, False)
    
    beam_column_pzspring = [[], [], []]
    for xx in range(n_hor_axis):
        for yy in range(n_vert_axis):
            node_ID = node_array[xx*n_vert_axis + yy]
            if yy != 0:
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
                ele_ID = col_i_node_ID if panel_zone else -1
                tmp_col = SpringBasedElementModifiedIMKSteelIShape(col_i_node_ID, col_j_node_ID, list_col[yy-1], geo_trans_ID,
                    col_mat_i, col_mat_j, N_G[yy-1], L_b=L_b_col[yy-1], ele_ID = ele_ID)
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
                    ele_ID = beam_i_node_ID if panel_zone else -1
                    tmp_beam = SpringBasedElementModifiedIMKSteelIShape(beam_i_node_ID, beam_j_node_ID, list_beam[yy-1], geo_trans_ID,
                        beam_mat_i, beam_mat_j, L_b=L_b_beam[xx-1], ele_ID = ele_ID)
                    tmp_beam.CreateMember()
                    beam_column_pzspring[0].append(deepcopy(tmp_beam))
            else:
                if fix_support: RigidSupport(node_ID)

    
    if show_plot:
        opsplt.plot_model("nodes", "elements")
        
    return beam_column_pzspring


def DefineSubassemblageNodes(beam_left_L_cl, beam_right_L_cl, col_top_L_cl, col_bottom_L_cl, depth_col, depth_beam,
    boundary_condition = True, show_plot = True):
    """
    Function that declares and initialises the grid nodes of an interior subassemblage. The panel zone geometry is defined by the two arguments
        depth_col and depth_beam. 

    @param beam_left_L_cl (float): Centerline length of the left beam (excluding the panel zone).
    @param beam_right_L_cl (float): Centerline length of the right beam (excluding the panle zone).
    @param col_top_L_cl (float): Centerline length of the top column (excluding the panel zone).
    @param col_bottom_L_cl (float): Centerline length of the bottom column (excluding the panel zone).
    @param depth_col (float): Depth of the columns for the panel zone.
    @param depth_beam (float): Depth of the beams for the panel zone.
    @param boundary_condition (bool, optional): Option to set already the boundary condition (bottom column pinned, beams fix only y movement).
        Defaults to True.
    @param show_plot (bool, optional): Option to show the plot of the nodes declared and initialised. Defaults to True.

    @exception NegativeValue: beam_left_L_cl needs to be positive.
    @exception NegativeValue: beam_right_L_cl needs to be positive.
    @exception NegativeValue: col_top_L_cl needs to be positive.
    @exception NegativeValue: col_bottom_L_cl needs to be positive.
    @exception NegativeValue: depth_col needs to be positive.
    @exception NegativeValue: depth_beam needs to be positive.

    @returns list: List with the nodes declared. 
    """
    # origin is the bottom left corner
    if beam_left_L_cl < 0: raise NegativeValue()
    if beam_right_L_cl < 0: raise NegativeValue()
    if col_top_L_cl < 0: raise NegativeValue()
    if col_bottom_L_cl < 0: raise NegativeValue()
    if depth_col < 0: raise NegativeValue()
    if depth_beam < 0: raise NegativeValue()

    node(12, 0.0, col_bottom_L_cl+depth_beam/2)
    node(21, beam_left_L_cl+depth_col/2, 0.0)
    node(22, beam_left_L_cl+depth_col/2, col_bottom_L_cl+depth_beam)
    node(23, beam_left_L_cl+depth_col/2, col_bottom_L_cl+depth_beam+col_top_L_cl)
    node(32, beam_left_L_cl+depth_col+beam_right_L_cl, col_bottom_L_cl+depth_beam/2)
    node_array = [12, 21, 22, 23, 32]
    
    if boundary_condition:
        fix(12, 0, 1, 0)
        fix(32, 0, 1, 0)
        fix(21, 1, 1, 0)

    if show_plot:
        plot_nodes(node_array, "Subassemblage geometry template with only nodes", True)
        plt.grid()

    return node_array


# def DefineRCSSubassemblage():
#     # WIP and experimental
    

