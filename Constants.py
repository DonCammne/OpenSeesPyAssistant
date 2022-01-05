"""Module with the values of a set of essential constants. They are consistent with the units defined in Units. \n
Carmine Schipani, 2021
"""

from OpenSeesPyAssistant.Units import *


TOL = 1.0e-6
TOL_INTEGRATION = 1.0e-12
ZERO = 1.0e-9                   # used when defining mass that is equal to 0 (avoid convergence problem)
G_CONST = 9.810*m_unit/s_unit**2
RIGID = 100.0 			        # multiply with a mechanical or geometrical value to have the corresponding infinitely rigid value  
MAX_ITER = 100
MAX_ITER_INTEGRATION = 50
