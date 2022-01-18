"""Module with useful functions (discretise curves, ID conventions, etc) \n
Carmine Schipani, 2021
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import openseespy.postprocessing.internal_plotting_functions as ipltf
from openseespy.opensees import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *


def ProgressingPercentage(max_iter, i: int, next_step: int, step = 10):
	"""
	Function that shows the progressing percentage of an iterative process.

	@param max_iter (int): Maximal number of interations
	@param i (int): Current iteration
	@param next_step (int): Next step of the percentage (set to 0 for the first iteration and then use the return parameter)
	@param step (int, optional): Size of the step (should be a fraction of 100). Defaults to 10.

	@returns int: The updated next step
	"""
	if i*100.0/(max_iter-1) >= next_step:
		print("The progression is {}%".format(next_step))
		return next_step + step

	return next_step


def DiscretizeLoadProtocol(SDR_LP: np.ndarray, nr_cycles_LP: np.ndarray, discr_first_cycle: int, plot = False, block = False, show_original_peaks = True):
	"""
	Discretized a cyclic load protocol keeping a similar discretisation step throughout the different cycles and keeping in the output the extremes (peaks).

	@param SDR_LP (np.ndarray): Array (1 dimension) that stores the peaks of the cycles.
		They needs to be only the positive peaks, beacuse this function will use them as the extreme for each cycle.
	@param nr_cycles_LP (np.ndarray): Array (1 dimension) that stores the number of cycles for every extreme declared in 'SDR_LP' and its countepart negative.
		They need to be positive integers.
	@param discr_first_cycle (int): The number of points from peak to peak (counting the two peaks) in the first cycle. It should be odd.
	@param plot (bool, optional): Option to show the plot of the discretized (and also the original peaks). Defaults to False.
	@param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
	@param show_original_peaks (bool, optional): Option to show the original peaks to check if the discretized curve is correct.
		The argument plot need to be True. Defaults to True.

	@exception WrongDimension: SDR_LP and nr_cycles_LP need to be of same length.
	@exception NegativeValue: SDR_LP needs to have only positive integers.
	@exception NegativeValue: nr_cycles_LP needs to have only positive integers.
    @exception NegativeValue: discr_first_cycle needs to be a positive integer.

	@returns np.array: Array (1 dimension) that stores the new discretized load protocol curve.
	"""
	if np.size(SDR_LP) != np.size(nr_cycles_LP): raise WrongDimension()
	if any(col < 0 for col in SDR_LP): raise NegativeValue()
	if any(col < 0 for col in nr_cycles_LP): raise NegativeValue()
	if discr_first_cycle < 0: raise NegativeValue()
	
	if discr_first_cycle % 2 == 0:
			discr_first_cycle = discr_first_cycle + 1
	discr_factor = discr_first_cycle / (SDR_LP[0]*2)
	discretized_LP = [0.0]
	x_val = []
	skip_x = 0
	for i in range(np.size(SDR_LP)):
		discr_i = math.ceil(discr_factor*SDR_LP[i]*2)-1;
		if discr_i % 2 == 0:
			discr_i = discr_i + 1
		length_tmp = int((discr_i+1)/2)
		tmp_up = np.linspace(0.0, SDR_LP[i], length_tmp)
		tmp_down = np.linspace(SDR_LP[i], 0.0, length_tmp)
		for j in range(int(nr_cycles_LP[i])):
			discretized_LP = np.append(discretized_LP, tmp_up[1:length_tmp])
			discretized_LP = np.append(discretized_LP, tmp_down[1:length_tmp])
			discretized_LP = np.append(discretized_LP, -tmp_up[1:length_tmp])
			discretized_LP = np.append(discretized_LP, -tmp_down[1:length_tmp])
		# for check of original peaks
		x_val.append(length_tmp-1+skip_x)
		skip_x = (length_tmp-1)*(4*(nr_cycles_LP[i]-1)+3)+x_val[-1]


	if plot:
		fig, ax = plt.subplots()
		ax.plot(discretized_LP, '-r', label='Discretised LP')

		ax.set(xlabel='Step number [-]', ylabel='Unit of the loading protocol', 
			title='Discretized loading protocol')
		ax.grid()

		if show_original_peaks:
			ax.plot(x_val, SDR_LP, 'ob', label='Original LP')
			ax.legend()

		if block:
			plt.show()

	return discretized_LP


def DiscretizeLinearly(LP: np.ndarray, discr: int, plot = False, block = False, show_original_LP = True):
	"""
	This function discretize the curve 'LP' given adding the number of point given by 'discr' between every point (linearly).

	@param LP (np.ndarray): Array (1 dimension) that stores the curve that needs to be discretized
	@param discr (int): The number of points to add between every two points of 'LP' (linearly)
	@param plot (bool, optional): Option to show the plot of the discretized (and also the original LP). Defaults to False.
	@param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
	@param show_original_LP (bool, optional): Option to show the original LP to check if the discretized curve is correct. Defaults to True.

	@returns np.ndarray: Array (1 dimension) that stores the new discretized load protocol.
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
		ax.plot(discr_LP, '-r', label='Discretised LP')

		ax.set(xlabel='Step number [-]', ylabel='Unit of the loading protocol', 
			title='Discretized loading protocol')
		ax.grid()

		if show_original_LP:
			x_val = np.arange(0, np.size(discr_LP), discr+1)
			ax.plot(x_val, LP, 'ob', label='Original LP')
			ax.legend()

		if block:
			plt.show()

	return discr_LP


def GridIDConvention(pier_axis: int, floor_axis: int, max_pier = -1, max_floor = -1):
	"""
	Function used to construct the ID of the nodes in the grid (first nodes that define the geometry of the model).
	The conventional grid node ID is xy, with x = the pier position 'pier_axis'; y = the floor position 'floor_axis'.

	@param pier_axis (int): The pier (or x) postion of the node. 
	@param floor_axis (int): The floor (or y) position of the node.
	@param max_pier (int, optional): Maximal pier position of the model (used to identify the number of digits).
		Defaults to -1, e.g. taken equal of pier_axis.
	@param max_floor (int, optional): Maximal floor position of the model (used to identify the number of digits).
		Defaults to -1, e.g. taken equal of floor_axis.

	@exception NameError: Work In Progress: only 9 floors/bays.
	@exception NegativeValue: The argument pier_axis needs to be a positive integer.
	@exception NegativeValue: The argument floor_axis needs to be a positive integer.
	@exception NegativeValue: The argument max_pier needs to be a positive integer if different from -1.
	@exception NegativeValue: The argument max_floor needs to be a positive integer if different from -1.
	@exception WrongArgument: The argument max_pier need to be equal or bigger to pier_axis
	@exception WrongArgument: The argument max_floor need to be equal or bigger to floor_axis

	@returns int: The grid node ID
	"""
	# Convention:
    # 		GridNodeID: 	xy 			with x = pier, y = floor 
	if pier_axis > 9 or floor_axis > 9 or max_pier > 9 or max_floor > 9: raise NameError("WIP: maximal 9 floors or bays")
	max_pier = pier_axis if max_pier == -1 else max_pier
	max_floor = floor_axis if max_floor == -1 else max_floor

	if pier_axis < 0: raise NegativeValue()
	if floor_axis < 0: raise NegativeValue()
	if max_pier < 0: raise NegativeValue()
	if max_floor < 0: raise NegativeValue()
	if max_pier < pier_axis: raise WrongArgument()
	if max_floor < floor_axis: raise WrongArgument()

	max_x_digits = int(math.log10(max_pier))+1
	max_y_digits = int(math.log10(max_floor))+1

	# return 10**(max_x_digits+max_y_digits) + pier_axis*10**max_y_digits + floor_axis # with 1 as prefix (to consider more than one digit per axis, but exceed max ID)
	return pier_axis*10**max_y_digits + floor_axis


def IDConvention(prefix: int, suffix: int, n_zeros_between = 0):
	"""
	Function used to construct IDs for elements and offgrid nodes.
	It appends to a positive integer number 'prefix' a number of zeros 'n_zeros_between' and at last another positive integer 'suffix'.
	The conventional element ID is xy(a)x'y'(a') with xya = the node ID in pier x, floor y and offgrid parameter a (optional);
		x'y'a' = the node ID in pier x', floor y' and offgrid parameter a' (optional).
	For more information on x and y, see GridIDConvention; for more information on a, see OffsetNodeIDConvention.

	@param prefix (int): Prefix of the new ID. For a vertical element it should be the left node ID;
		for an horizontal one it should be the bottom node.
	@param suffix (int): Suffix of the new ID. For a vertical element it should be the right node ID;
		for an horizontal one it should be the top node.
	@param n_zeros_between (int, optional): Number of zeros to add between the two nodes. Defaults to 0.

	@exception NegativeValue: The argument prefix needs to be a positive integer.
	@exception NegativeValue: The argument suffix needs to be a positive integer.
	@exception NegativeValue: The argument n_zeros_between needs to be a positive integer.

	@returns int: The combined ID
	"""
	# Convention:
    #		ElementID:		xy(a)x'y'(a')	with xy(a) = NodeID i and x'y'(a') = NodeID j
    #		TrussID:		xy(a)x'y'(a')	with xy(a) = NodeID i and x'y'(a') = NodeID j
    #		Spring:			xy(a)x'y'(a')	with xy(a) = NodeID i and x'y'(a') = NodeID j
	if prefix < 0: raise NegativeValue()
	if suffix < 0: raise NegativeValue()
	if n_zeros_between < 0: raise NegativeValue()
	
	return int(str(prefix*10**n_zeros_between) + str(suffix))


def OffsetNodeIDConvention(node_ID: int, orientation: str, position_i_or_j: str):
	"""
	Function used to add node on top of existing ones in the extremes of memebers with springs.

	@param node_ID (int): Node that we refer to.
	@param orientation (str): Orientation of the memeber. Can be 'vertical' or 'horizontal'.
	@param position_i_or_j (str): Position at the start 'i' (left or bottom)
		or at the end 'j' (right or top) of 'node_ID' in the member.

	@exception NegativeValue: The argument node_ID needs to be a positive integer.
	@exception WrongArgument: The argument position_i_or_j needs to be 'i' or 'j'
	@exception WrongArgument: The argument orientation needs to be 'vertical' or 'horizontal' 

	@returns int: The combined ID
	"""
	# Convention:																							o xy
    # 		GridNodeID: 	 xy 			with x = pier, y = floor 										o xy7		
    #	AdditionalNodeID:    xya			with x = pier, y = floor, a:  o xy  xy2 o-------o x'y3  x'y o	|		
    #	PanelZoneNodeID: 	 xy(a)a 	see MemberModel for the panel zone ID convention					|	
    #																								   		|		
    #																								   		o xy'6	
    #																										o xy'
	if node_ID < 1: raise NegativeValue()
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


def NodesOrientation(iNode_ID: int, jNode_ID: int):
	"""
	Function that finds the orientation of the vector with direction 'jNode_ID''iNode_ID'.
	If the the nodes are on top of each other, the function returns 'zero_length'.

	@param iNode_ID (int): Node i.
	@param jNode_ID (int): Node j.

	@exception NegativeValue: The argument iNode_ID needs to be a positive integer.
	@exception NegativeValue: The argument jNode_ID needs to be a positive integer.

	@returns str: The orientation of the vector.
	"""
	if iNode_ID < 1: raise NegativeValue()
	if jNode_ID < 1: raise NegativeValue()

	iNode = np.array(nodeCoord(iNode_ID))
	jNode = np.array(nodeCoord(jNode_ID))
	if abs(iNode[0]-jNode[0]) + abs(iNode[1]-jNode[1]) == 0:
		return "zero_length"
	elif abs(iNode[0]-jNode[0]) < abs(iNode[1]-jNode[1]):
		return "vertical"
	else:
		return "horizontal"


def plot_member(element_array: list, member_name = "Member name not defined", show_element_ID = True, show_node_ID = True):
	"""
	Function that plots a set of elements. It can be used to check the correctness of a part of the model or of a member.
	If the entire model need to be plotted, use instead 'plot_model("nodes", "elements")' from openseespy.postprocessing.Get_Rendering. \n
    Inspired by plot_model written by Anurag Upadhyay and Christian Slotboom.

	@param element_array (list): An array (list of lists of one dimensions and length = 3) that store the element and nodes IDs.
		An element is stored in one list with 3 entries: the element ID, node i ID and node j ID.
	@param member_name (str, optional): The name of what is plotted. Defaults to "Member name not defined".
	@param show_element_ID (bool, optional): Option to show the element IDs. Defaults to True.
	@param show_node_ID (bool, optional): Option to show the nodes IDs. Defaults to True.

	@exception WrongDimension: element_array needs to be non-empty.
	@exception WrongDimension: The number of entries in the lists inside the argument element_array need to be 3.

	@returns matplotlib.axes._subplots.AxesSubplo: The figure's wrappr, useful to customise the plot (change axes label, etc).
	"""
	if len(element_array) == 0: raise WrongArgument()
	if len(element_array[0]) != 3: raise WrongDimension()

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
					__plt_node(Nodes[1], track_node, jNode, ax, node_text_style)
		else:
			print("Too many nodes in this elemnet (see shell elements)")
		
	ax.set_xlabel('x [{}]'.format(length_unit))
	ax.set_ylabel('y [{}]'.format(length_unit))
	plt.title("Visualisation of: {}".format(member_name))
	plt.axis('equal')
	return ax


def plot_nodes(nodes_array: list, name = "Not defined", show_node_ID = True):
	"""
	Function that plots a set of nodes. It can be used to check the correctness of the model's geometry.
	If the entire model need to be plotted, use instead 'plot_model("nodes", "elements")' from openseespy.postprocessing.Get_Rendering.

	@param nodes_array (list): List of 1 dimension with the IDs of the nodes to be displayed.
	@param name (str, optional): Name that describe what the plot will show. Defaults to "Not defined".
	@param show_node_ID (bool, optional): Option to show the node IDs. Defaults to True.

	@exception WrongArgument: nodes_array needs to be non-empty.

	@returns (matplotlib.axes._subplots.AxesSubplot): The figure's wrapper, useful to customise the plot (change axes label, etc).
	"""
	if len(nodes_array) == 0: raise WrongArgument()

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
	"""
	PRIVATE FUNCTION. Used to plot the nodes in a controlled way (no repetition, position of the IDs, text style).

	@param nodeID (int): The ID of the node.
	@param track_node (dict): A dictionary used to avoid plotting a node multiple times.
	@param NodeXY (list): List of dimension 1, length 2 with the position of the node.
	@param ax (matplotlib.axes._subplots.AxesSubplot): The figure's wrappr.
	@param node_text_style (dict): Dictionary for the text style.
	@param x_off (int, optional): Offset in x direction. Defaults to 0.
	@param y_off (int, optional): Offset in y direction. Defaults to 0.
	@param h_align (str, optional): Horizontal alignment ('left' or 'right'). Defaults to 'left'.
	@param v_align (str, optional): Vertical alignment ('center', 'top' or 'bottom'). Defaults to 'top'.
	"""
	if not nodeID in track_node:
		track_node[nodeID] = True
		ax.text(NodeXY[0]+x_off, NodeXY[1]+y_off, nodeID,**node_text_style, horizontalalignment=h_align, verticalalignment=v_align)


class IDGenerator():
	"""Class that manage the ID generation.
	USE ONLY IF EVERY NODE IS DEFINED BY THE USER (because the OpenSeesPyAssistant modules use the convention defined in the functions above).
	"""
	def __init__(self):
		"""The class constructor.
		"""
		self.current_node_ID = 0
		self.current_element_ID = 0
		self.current_mat_ID = 0
		self.current_fiber_ID = 0

	def GenerateIDNode(self):
		"""
		Method that generate a unique node ID.

		@returns int: The node ID.
		"""
		self.current_node_ID = self.current_node_ID + 1
		return self.current_node_ID

	def GenerateIDElement(self):
		"""
		Method that generate a unique element ID.

		@returns int: The element ID.
		"""
		self.current_element_ID = self.current_element_ID + 1
		return self.current_element_ID

	def GenerateIDMat(self):
		"""
		Method that generate a unique material ID.

		@returns int: The material ID.
		"""
		self.current_mat_ID = self.current_mat_ID + 1
		return self.current_mat_ID

	def GenerateIDFiber(self):
		"""
		Method that generate a unique fiber ID.

		@returns int: The fiber ID.
		"""
		self.current_fiber_ID = self.current_fiber_ID + 1
		return self.current_fiber_ID


