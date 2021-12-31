# Module with useful functions
#	Carmine Schipani, 2021

import math
import numpy as np
import matplotlib.pyplot as plt
import openseespy.postprocessing.internal_plotting_functions as ipltf
from openseespy.opensees import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *


def ProgressingPercentage(max_iter, i, next_step, step = 10):
	"""Function that shows the progressing percentage of an iterative process.

	Args:
		max_iter (int): Maximal number of interations
		i (int): Current iteration
		next_step (int): Next step of the percentage (set to 0 for the first iteration and then use the return parameter)
		step (int, optional): Size of the step (should be a fraction of 100). Defaults to 10.

	Returns:
		int: The updated next step
	"""

	if i*100.0/(max_iter-1) >= next_step:
		print("The progression is {}%".format(next_step))
		return next_step + step

	return next_step


def DiscretizeLoadProtocol(SDR_LP: np.ndarray, nr_cycles_LP: np.ndarray, discr_first_cycle: int, plot = False, block = False):
	"""Discretized a load protocol maintening a similar discretisation throughout the different cycles and keeping in the output the extremes (peaks).

	Args:
		SDR_LP (np.ndarray): Array (1 dimension) that stores the peaks of the cycles (positive)
		nr_cycles_LP (np.ndarray): Array (1 dimension) that stores the number of cycles for every extreme declared before the
		discr_first_cycle (int): The number of points from peak to peak (counting the two peaks). It should be odd

	Returns:
		np.array: Array (1 dimension) that stores the new discretized load protocol
	"""
	#TODO: check that SDR_LP and nr_cycles_LP are nonnegative, that nr_cycles_LP should have only int, they are 1 dimensions and size(SDR_LP)==size(nr_cycles_LP)

	if discr_first_cycle % 2 == 0:
			discr_first_cycle = discr_first_cycle + 1
	discr_factor = discr_first_cycle / (SDR_LP[0]*2)
	discretized_LP = [0.0]
	for i in range(np.size(SDR_LP)):
		discr_i = math.ceil(discr_factor*SDR_LP[i]*2)-1;
		if discr_i % 2 == 0:
			discr_i = discr_i + 1
		length_tmp = int((discr_i+1)/2)
		tmp_up = np.linspace(0.0, SDR_LP[i], length_tmp)
		tmp_down = np.linspace(SDR_LP[i], 0.0, length_tmp)
		for j in range(nr_cycles_LP[i]):
			discretized_LP = np.append(discretized_LP, tmp_up[1:length_tmp])
			discretized_LP = np.append(discretized_LP, tmp_down[1:length_tmp])
			discretized_LP = np.append(discretized_LP, -tmp_up[1:length_tmp])
			discretized_LP = np.append(discretized_LP, -tmp_down[1:length_tmp])
	
	if plot:
		fig, ax = plt.subplots()
		ax.plot(discretized_LP)

		ax.set(xlabel='Step number [-]', ylabel='Unit of the loading protocol', 
			title='Discretized loading protocol')
		ax.grid()

		if block:
			plt.show()

	return discretized_LP


def DiscretizeLinearly(LP: np.ndarray, discr: int, plot = False, block = False):
	"""
	This function creates a discretized LP with the number of point given by discr between every point from LP (linearly).
	
	LP : np.array
		Array (1 dimension) that stores the curve that needs to be discretized
	discr : int
		The number of points to add between every two points of LP (linearly)
	"""

	#TODO: check discr nonnegative int and LP 1 dimension
	
	# Define the new discretized LP
	length = 1 + (np.size(LP)-1) * (discr+1)
	discr_LP = np.zeros(length)

	# Performa manually the first iteration
	yprev = LP[0]
	x = [0, 1]
	discr_LP[0] = yprev
	iter = 0

	# add the other points and the discretized ones
	for ynext in LP[1:]:
		y = [yprev, ynext]

		# Compute new points
		xnew = np.linspace(x[0], x[1], discr+2)
		ynew = np.interp(xnew[1:], x, y)

		# Add to the recording vector discr_LP
		index = np.array(np.arange(discr+1)+1+iter)
		discr_LP[index] = ynew

		# Prepare for next iteration
		yprev = ynext
		iter = iter + discr + 1

	if plot:
		fig, ax = plt.subplots()
		ax.plot(discr_LP)

		ax.set(xlabel='Step number [-]', ylabel='Unit of the loading protocol', 
			title='Discretized loading protocol')
		ax.grid()

		if block:
			plt.show()

	return discr_LP


def IDConvention(iNodeID: int, jNodeID: int, n_zeros_between: int = 0):
    # Convention:																							
    # 		GridNodeID: 	1xy 			with x = pier, y = floor 											o 1xy7	|	
    #	AdditionalNodeID:   1xya			with x = pier, y = floor, a:  --o 1xy  1xy2 o-----o 2xy3  2xy o--	|		o 1xy
    #																									    	|	
    #																									    	|		o 1xy
    #																									    	o 1xy6	|
    #		ElementID:		1xy(a)1xy(a)	with 1xy(a) = NodeID i and j
    #		TrussID:		1xy(a)1xy(a)	with 1xy(a) = NodeID i and j
    #		PDeltaColID:	1xy(a)1xy(a)	with 1xy(a) = NodeID i and j
    #		Spring:			1xy(a)1xy(a)	with 1xy(a) = NodeID i and j
	if n_zeros_between < 0: raise NegativeValue()
	
	return int(str(iNodeID*10**n_zeros_between) + str(jNodeID))


def OffsetNodeIDConvention(node_ID: int, orientation: str, position_i_or_j: str):
	if position_i_or_j != "i" and position_i_or_j != "j": raise WrongArgument()
	if orientation == "vertical":
		if position_i_or_j == "i":
			return IDConvention(node_ID, 6)
		else:
			return IDConvention(node_ID, 7)
	elif orientation == "horizontal":
		if position_i_or_j == "i":
			return IDConvention(node_ID, 2)
		else:
			return IDConvention(node_ID, 3)
	else: raise WrongArgument()


def GridIDConvention(n_x_axis: int, n_y_axis: int, max_n_x = -1, max_n_y = -1):
    # Convention:																							
    # 		NodeID: 		1xy 			with x = pier, y = floor

    if max_n_x != -1 and max_n_x < 0: raise NegativeValue()
    if max_n_y != -1 and max_n_y < 0: raise NegativeValue()

    max_n_x = n_x_axis if max_n_x == -1 else max_n_x
    max_n_y = n_y_axis if max_n_y == -1 else max_n_y
    if max_n_x < n_x_axis: raise WrongArgument()
    if max_n_y < n_y_axis: raise WrongArgument()

    max_x_digits = int(math.log10(max_n_x))+1
    max_y_digits = int(math.log10(max_n_y))+1
    
    return 10**(max_x_digits+max_y_digits) + n_x_axis*10**max_y_digits + n_y_axis


def NodesOrientation(iNode_ID, jNode_ID):
	iNode = np.array(nodeCoord(iNode_ID))
	jNode = np.array(nodeCoord(jNode_ID))
	if abs(iNode[0]-jNode[0]) + abs(iNode[1]-jNode[1]) == 0:
		return "zero_length"
	elif abs(iNode[0]-jNode[0]) < abs(iNode[1]-jNode[1]):
		return "vertical"
	else:
		return "horizontal"


def plot_member(element_array: list, member_name = "Member name not defined", show_element_ID = True, show_node_ID = True):
    node_style = {'color':'black', 'marker':'o', 'facecolor':'black','linewidth':0.}
    node_text_style = {'fontsize':8, 'fontweight':'regular', 'color':'green'} 
    track_node = {}

    if show_element_ID:
        show_e_ID = 'yes'
    else:
        show_e_ID = 'no'

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for ele in element_array:
        eleTag = int(ele[0])
        Nodes =ele[1:]

        if len(Nodes) == 2:
            # 2D element
            iNode = np.array(nodeCoord(Nodes[0]))
            jNode = np.array(nodeCoord(Nodes[1]))
            ipltf._plotBeam2D(iNode, jNode, ax, show_e_ID, eleTag, "solid")
            ax.scatter(*iNode, **node_style)
            ax.scatter(*jNode, **node_style)
            if show_node_ID:
                if abs(sum(iNode - jNode)) > 1e-6:
                    # beam-col
                    __plt_node(Nodes[0], track_node, iNode, ax, node_text_style)
                    __plt_node(Nodes[1], track_node, jNode, ax, node_text_style, h_align='right', v_align='bottom')
                else:
                    # zerolength
                    __plt_node(Nodes[0], track_node, iNode, ax, node_text_style, h_align='right')
                    __plt_node(Nodes[1], track_node, jNode, ax, node_text_style, v_align='bottom')
        else:
            print("Too many nodes in this elemnet (see shell elements)")
        
    ax.set_xlabel('x [{}]'.format(length_unit))
    ax.set_ylabel('y [{}]'.format(length_unit))
    plt.title("Visualisation of: {}".format(member_name))
    plt.axis('equal')
    return ax


def plot_nodes(nodes_array: list, name = "Not defined", show_node_ID = True):
    node_style = {'color':'black', 'marker':'o', 'facecolor':'black','linewidth':0.}
    node_text_style = {'fontsize':8, 'fontweight':'regular', 'color':'green'} 
    track_node = {}

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for node_ID in nodes_array:
        node_xy = np.array(nodeCoord(node_ID))
        ax.scatter(*node_xy, **node_style)
        if show_node_ID:
            __plt_node(node_ID, track_node, node_xy, ax, node_text_style)

    ax.set_xlabel('x [{}]'.format(length_unit))
    ax.set_ylabel('y [{}]'.format(length_unit))
    plt.title("Visualisation of: {}".format(name))
    plt.axis('equal')
    return ax


def __plt_node(nodeID: int, track_node: dict, NodeXY, ax, node_text_style, x_off = 0, y_off = 0, h_align = 'left', v_align='top'):
	if not nodeID in track_node:
		track_node[nodeID] = True
		ax.text(NodeXY[0]+x_off, NodeXY[1]+y_off, nodeID,**node_text_style, horizontalalignment=h_align, verticalalignment=v_align)


