# Module that set the value for a set of essential constants
#   Carmine Schipani, 20

from OpenSeesPyHelper.Units import *

"""Module with the values of a set of essential constants. 
"""

TOL = 1.0e-6
ZERO = 1.0e-9 # used when defining mass that is equal to 0 (avoid convergence problem)
G_CONST = 9.810*m_unit/s_unit**2
RIGID = 100.0 			# multiply with the biggest value to have an infinitely rigid element
MAX_ITER = 100
