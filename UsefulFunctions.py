# Module with useful functions
#	Carmine Schipani, 2021

import math
import numpy as np


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

	# else
	return next_step

def DiscretizeLoadProtocol(SDR_LP: np.ndarray, nr_cycles_LP: np.ndarray, discr_first_cycle: int):
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
	return discretized_LP

def DiscretizeLinearly(LP: np.ndarray, discr: int):
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

	return discr_LP


