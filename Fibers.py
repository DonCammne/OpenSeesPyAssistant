"""
Module for the fibers (rectangular, circular and I shape).
Carmine Schipani, 2021
"""

from openseespy.opensees import *
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Wedge
import numpy as np
from copy import copy, deepcopy
from OpenSeesPyAssistant.Section import *
from OpenSeesPyAssistant.DataManagement import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *
from OpenSeesPyAssistant.MaterialModels import *


class Fibers(DataManagement):
    """
    Parent abstract class for the storage and manipulation of a fiber's information (mechanical
        and geometrical parameters, etc) and initialisation in the model.

    @param DataManagement: Parent abstract class.
    """
    pass


class FibersRect(Fibers):
    """
	Class that stores funcions, material properties, geometric and mechanical parameters for a rectangular RC fiber section.
    Coordinates: plotting coordinte (x, y) = fiber section coordinate (z, y) = (-x, y). For more information, see the OpenSeesPy documentation.

    @param Fibers: Parent abstract class.
    """
    def __init__(self, ID: int, b, d, Ay, D_hoops, e, unconf_mat_ID: int, conf_mat_ID: int, bars_mat_ID: int,
        bars_x: np.ndarray, ranges_y: np.ndarray, discr_core: list, discr_cover_lateral: list, discr_cover_topbottom: list, GJ = 0.0):
        """
        Constructor of the class.

        @param ID (int): Unique fiber section ID.
        @param b (float): Width of the section.
        @param d (float): Depth of the section.
        @param Ay (float): Area of one vertical reinforcing bar.
        @param D_hoops (float): Diameter of the hoops.
        @param e (float): Concrete cover.
        @param unconf_mat_ID (int): ID of material model that will be assigned to the unconfined fibers.
        @param conf_mat_ID (int): ID of material model that will be assigned to the confined fibers.
        @param bars_mat_ID (int): ID of material model that will be assigned to the reinforcing bars fibers.
        @param bars_x (np.ndarray): Array with a range of aligned vertical reinforcing bars for each row in x direction.
            Distances from border to bar centerline, bar to bar centerlines and finally bar centerline to border in the x direction (aligned).
            Starting from the left to right, from the top range to the bottom one.
            The number of bars for each range can vary; in this case, add this argument when defining the array " dtype = object"
        @param ranges_y (np.ndarray): Array of dimension 1 with the position or spacing in y of the ranges in bars_x.
            Distances from border to range centerlines, range to range centerlines and finally range centerline to border in the y direction.
            Starting from the top range to the bottom one.
        @param discr_core (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the confined core.
        @param discr_cover_lateral (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the lateral unconfined cover.
        @param discr_cover_topbottom (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the top and bottom unconfined cover.
        @param GJ (float, optional): Linear-elastic torsional stiffness assigned to the section. Defaults to 0.0, assume no torsional stiffness.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: b needs to be positive.
        @exception NegativeValue: d needs to be positive.
        @exception NegativeValue: Ay needs to be positive.
        @exception NegativeValue: D_hoops needs to be positive.
        @exception NegativeValue: e needs to be positive.
        @exception NegativeValue: unconf_mat_ID needs to be a positive integer.
        @exception NegativeValue: conf_mat_ID needs to be a positive integer.
        @exception NegativeValue: bars_mat_ID needs to be a positive integer.
        @exception WrongDimension: Number of rows in the list bars_x needs to be the same of the length of ranges_y - 1.
        @exception InconsistentGeometry: The sum of the distances for each row in bars_x should be equal to the section's width (tol = 5 mm).
        @exception InconsistentGeometry: The sum of the distances in ranges_y should be equal to the section's depth (tol = 5 mm).
        @exception InconsistentGeometry: e should be smaller than half the depth and the width of the section.
        @exception WrongDimension: discr_core has a length of 2.
        @exception WrongDimension: discr_cover_lateral has a length of 2.
        @exception WrongDimension: discr_cover_topbottom has a length of 2.
        @exception NegativeValue: GJ needs to be positive.
        """
        # Check
        if ID < 1: raise NegativeValue()
        if b < 0: raise NegativeValue()
        if d < 0: raise NegativeValue()
        if Ay < 0: raise NegativeValue()
        if D_hoops < 0: raise NegativeValue()
        if e < 0: raise NegativeValue()
        if unconf_mat_ID < 1: raise NegativeValue()
        if conf_mat_ID < 1: raise NegativeValue()
        if bars_mat_ID < 1: raise NegativeValue()
        if np.size(bars_x) != np.size(ranges_y)-1: raise WrongDimension()
        geometry_tol = 5*mm_unit
        for bars in bars_x:
            if abs(np.sum(bars) - b) > geometry_tol: raise InconsistentGeometry()
        if abs(np.sum(ranges_y) - d) > geometry_tol: raise InconsistentGeometry()
        if e > b/2 or e > d/2: raise InconsistentGeometry()
        if len(discr_core) != 2: raise WrongDimension()
        if len(discr_cover_lateral) != 2: raise WrongDimension()
        if len(discr_cover_topbottom) != 2: raise WrongDimension()
        if GJ < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.b = b
        self.d = d
        self.Ay = Ay
        self.D_hoops = D_hoops
        self.e = e
        self.unconf_mat_ID = unconf_mat_ID
        self.conf_mat_ID = conf_mat_ID
        self.bars_mat_ID = bars_mat_ID
        self.bars_x = deepcopy(bars_x)
        self.ranges_y = copy(ranges_y)
        self.discr_core = copy(discr_core)
        self.discr_cover_lateral = copy(discr_cover_lateral)
        self.discr_cover_topbottom = copy(discr_cover_topbottom)
        self.GJ = GJ

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.Initialized = False
        self.ReInit()

    def ReInit(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        # Memebers
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Parameters
        z1 = self.b/2
        y1 = self.d/2
        zc = z1-self.e-self.D_hoops/2
        yc = y1-self.e-self.D_hoops/2

        # Create the concrete core fibers
        core = [-yc, -zc, yc, zc]
        core_cmd = ['patch', 'rect', self.conf_mat_ID, *self.discr_core, *core]

        # Create the concrete cover fibers (bottom left, top right)
        cover_up = [yc, -z1, y1, z1]
        cover_down = [-y1, -z1, -yc, z1]
        cover_left = [-yc, zc, yc, z1]
        cover_right = [-yc, -z1, yc, -zc]
        cover_up_cmd = ['patch', 'rect', self.unconf_mat_ID, *self.discr_cover_topbottom, *cover_up]
        cover_down_cmd = ['patch', 'rect', self.unconf_mat_ID, *self.discr_cover_topbottom, *cover_down]
        cover_left_cmd = ['patch', 'rect', self.unconf_mat_ID, *self.discr_cover_lateral, *cover_left]
        cover_right_cmd = ['patch', 'rect', self.unconf_mat_ID, *self.discr_cover_lateral, *cover_right]
        self.fib_sec = [['section', 'Fiber', self.ID, '-GJ', self.GJ], 
            core_cmd, cover_up_cmd, cover_down_cmd, cover_left_cmd, cover_right_cmd]
        
        # Create the reinforcing fibers (top, middle, bottom)
        # NB: note the order of definition of bars_x and ranges_y
        nr_bars = 0
        for range in self.bars_x:
            nr_bars += np.size(range)-1
        rebarY = -np.cumsum(self.ranges_y[0:-1]) + y1
        self.rebarYZ = np.zeros((nr_bars, 2))

        iter = 0
        for ii, Y in enumerate(rebarY):
            rebarZ = -np.cumsum(self.bars_x[ii][0:-1]) + z1
            for Z in rebarZ:
                self.rebarYZ[iter, :] = [Y, Z]
                iter = iter + 1
        
        for YZ in self.rebarYZ:
            self.fib_sec.append(['layer', 'bar', self.bars_mat_ID, self.Ay, *YZ])

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        self.data = [["INFO_TYPE", "FibersRect"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["b", self.b],
            ["d", self.d],
            ["Ay", self.Ay],
            ["D_hoops", self.D_hoops],
            ["e", self.e],
            ["GJ", self.GJ],
            ["conf_mat_ID", self.conf_mat_ID],
            ["discr_core", self.discr_core],
            ["unconf_mat_ID", self.unconf_mat_ID],
            ["discr_cover_topbottom", self.discr_cover_topbottom],
            ["discr_cover_lateral", self.discr_cover_lateral],
            ["bars_mat_ID", self.bars_mat_ID],
            ["bars_x", self.bars_x],
            ["ranges_y", self.ranges_y],
            ["Initialized", self.Initialized]]

    def ShowInfo(self, plot = False, block = False):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the fiber. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
        """
        print("")
        print("Requested info for FibersRect, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print("Base b = {} mm and depth d = {} mm".format(self.b/mm_unit, self.d/mm_unit))
        print("Confined material model ID = {}".format(self.conf_mat_ID))
        print("Unconfined material model ID = {}".format(self.unconf_mat_ID))
        print("Bars material model ID = {}".format(self.bars_mat_ID))
        print("Discretisation in the core [IJ or x/z dir, JK or y dir] = {}".format(self.discr_core))
        print("Discretisation in the lateral covers [IJ or x/z dir, JK or y dir] = {}".format(self.discr_cover_lateral))
        print("Discretisation in the top and bottom covers [IJ or x/z dir, JK or y dir] = {}".format(self.discr_cover_topbottom))
        print("")

        if plot:
            plot_fiber_section(self.fib_sec, matcolor=['#808080', '#D3D3D3', 'k'])
            
            if block:
                plt.show()


    def CreateFibers(self):
        """
        Method that initialises the fiber by calling the OpenSeesPy commands.
        """
        create_fiber_section(self.fib_sec)
        self.Initialized = True
        self.UpdateStoredData()


class FibersRectRCRectShape(FibersRect):
    """
    Class that is the children of FibersRect and combine the class RCRectShape (section) to retrieve the information needed.  

    @param FibersRect: Parent class.
    """
    def __init__(self, ID: int, section: RCRectShape, unconf_mat_ID: int, conf_mat_ID: int, bars_mat_ID: int,
        discr_core: list, discr_cover_lateral: list, discr_cover_topbottom: list, GJ=0):
        """
        Constructor of the class.

        @param ID (int): Unique fiber section ID.
        @param section (RCRectShape): RCRectShape section object.
        @param unconf_mat_ID (int): ID of material model that will be assigned to the unconfined fibers.
        @param conf_mat_ID (int): ID of material model that will be assigned to the confined fibers.
        @param bars_mat_ID (int): ID of material model that will be assigned to the reinforcing bars fibers.
        @param discr_core (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the confined core.
        @param discr_cover_lateral (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the lateral unconfined core.
        @param discr_cover_topbottom (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the top and bottom unconfined core.
        @param GJ (float, optional): Linear-elastic torsional stiffness assigned to the section. Defaults to 0.0, assume no torsional stiffness.
        """
        self.section = deepcopy(section)
        super().__init__(ID, section.b, section.d, section.Ay, section.D_hoops, section.e, unconf_mat_ID, conf_mat_ID, bars_mat_ID,
            section.bars_position_x, section.bars_ranges_position_y, discr_core, discr_cover_lateral, discr_cover_topbottom, GJ=GJ)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()


class FibersCirc(Fibers):
    """
	Class that stores funcions, material properties, geometric and mechanical parameters for a circular RC fiber section.
    Coordinates: plotting coordinte (x, y) = fiber section coordinate (z, y) = (-x, y). For more information, see the OpenSeesPy documentation.

    @param Fibers: Parent abstract class.
    """
    def __init__(self, ID: int, b, e, D_bars, Ay, n_bars, D_hoops, unconf_mat_ID: int, conf_mat_ID: int, bars_mat_ID: int,
        discr_core: list, discr_cover: list, alpha_i = 0.0, GJ = 0.0):
        """
        Constructor of the class.

        @param ID (int): Unique fiber section ID.
        @param b (float): Width of the section.
        @param e (float): Concrete cover.
        @param D_bars (float): Diameter of vertical reinforcing bars.
        @param Ay (float): Area of one vertical reinforcing bar.
        @param n_bars (float): Number of reinforcement (allow float for computing the equivalent n_bars with different reinforcement areas).
        @param D_hoops (float): Diameter of the hoops.
        @param unconf_mat_ID (int): ID of material model that will be assigned to the unconfined fibers.
        @param conf_mat_ID (int): ID of material model that will be assigned to the confined fibers.
        @param bars_mat_ID (int): ID of material model that will be assigned to the reinforcing bars fibers.
        @param discr_core (list): List with two entries: number of subdivisions (fibers) in the circumferential direction (number of wedges),
            number of subdivisions (fibers) in the radial direction (number of rings) for the confined core.
        @param discr_cover (list): List with two entries: number of subdivisions (fibers) in the circumferential direction (number of wedges),
            number of subdivisions (fibers) in the radial direction (number of rings) for the unconfined cover.
        @param alpha_i (float, optional): Angle in deg of the first vertical rebars with respect to the y axis, counterclockwise. Defaults to 0.0.
        @param GJ (float, optional): Linear-elastic torsional stiffness assigned to the section. Defaults to 0.0, assume no torsional stiffness.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: b needs to be positive.
        @exception NegativeValue: e needs to be positive.
        @exception InconsistentGeometry: e can't be bigger than half of the width b.
        @exception NegativeValue: D_bars needs to be positive.
        @exception NegativeValue: Ay needs to be positive.
        @exception NegativeValue: n_bars needs to be positive.
        @exception NegativeValue: D_hoops needs to be positive.
        @exception NegativeValue: unconf_mat_ID needs to be a positive integer.
        @exception NegativeValue: conf_mat_ID needs to be a positive integer.
        @exception NegativeValue: bars_mat_ID needs to be a positive integer.
        @exception WrongDimension: discr_core has a length of 2.
        @exception WrongDimension: discr_cover has a length of 2.
        @exception NegativeValue: GJ needs to be positive.
        """
        # Check
        if ID < 1: raise NegativeValue()
        if b < 0: raise NegativeValue()
        if e < 0: raise NegativeValue()
        if e > b/2: raise InconsistentGeometry()
        if D_bars < 0: raise NegativeValue()
        if Ay < 0: raise NegativeValue()
        if n_bars < 0: raise NegativeValue()
        if D_hoops < 0: raise NegativeValue()
        if unconf_mat_ID < 1: raise NegativeValue()
        if conf_mat_ID < 1: raise NegativeValue()
        if bars_mat_ID < 1: raise NegativeValue()
        if len(discr_core) != 2: raise WrongDimension()
        if len(discr_cover) != 2: raise WrongDimension()
        if GJ < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.b = b
        self.e = e
        self.D_bars = D_bars
        self.Ay = Ay
        self.n_bars = n_bars
        self.D_hoops = D_hoops
        self.unconf_mat_ID = unconf_mat_ID
        self.conf_mat_ID = conf_mat_ID
        self.bars_mat_ID = bars_mat_ID
        self.discr_core = copy(discr_core)
        self.discr_cover = copy(discr_cover)
        self.alpha_i = alpha_i
        self.GJ = GJ

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.Initialized = False
        self.ReInit()

    def ReInit(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        # Memebers
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Parameters
        self.r_bars = self.b/2 - self.e - self.D_hoops - self.D_bars/2
        self.r_core = self.b/2 - self.e - self.D_hoops/2

        # Create the concrete core fibers
        core_cmd = ['patch', 'circ', self.conf_mat_ID, *self.discr_core, 0, 0, 0, self.r_core]

        # Create the concrete cover fibers
        cover_cmd = ['patch', 'circ', self.unconf_mat_ID, *self.discr_cover, 0, 0, self.r_core, self.b/2]
        self.fib_sec = [['section', 'Fiber', self.ID, '-GJ', self.GJ], 
            core_cmd, cover_cmd]
        
        # Create the reinforcing fibers
        bars_cmd = ['layer', 'circ', self.bars_mat_ID, self.n_bars, self.Ay, 0, 0, self.r_bars, self.alpha_i]
        self.fib_sec.append(bars_cmd)

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        self.data = [["INFO_TYPE", "FibersCirc"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["b", self.b],
            ["e", self.e],
            ["r_core", self.r_core],
            ["D_bars", self.D_bars],
            ["Ay", self.Ay],
            ["n_bars", self.n_bars],
            ["r_bars", self.r_bars],
            ["D_hoops", self.D_hoops],
            ["alpha_i", self.alpha_i],
            ["GJ", self.GJ],
            ["conf_mat_ID", self.conf_mat_ID],
            ["discr_core", self.discr_core],
            ["unconf_mat_ID", self.unconf_mat_ID],
            ["discr_cover", self.discr_cover],
            ["bars_mat_ID", self.bars_mat_ID],
            ["Initialized", self.Initialized]]
        

    def ShowInfo(self, plot = False, block = False):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the material model. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
        """
        print("")
        print("Requested info for FibersCirc, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print("Base b = {} mm and concrete cover e = {} mm".format(self.b/mm_unit, self.e/mm_unit))
        print("Radius of the confined core r_core = {} mm, radius of the bars range r_bars = {} mm and initial angle alpha_i = {} deg".format(self.r_core/mm_unit, self.r_bars/mm_unit, self.alpha_i))
        print("Confined material model ID = {}".format(self.conf_mat_ID))
        print("Unconfined material model ID = {}".format(self.unconf_mat_ID))
        print("Bars material model ID = {}".format(self.bars_mat_ID))
        print("Discretisation in the core [number of wedges, number of rings] = {}".format(self.discr_core))
        print("Discretisation in the lateral covers [number of wedges, number of rings] = {}".format(self.discr_cover))
        print("")

        if plot:
            plot_fiber_section(self.fib_sec, matcolor=['#808080', '#D3D3D3', 'k'])
            
            if block:
                plt.show()


    def CreateFibers(self):
        """
        Method that initialise the fiber by calling the OpenSeesPy commands.
        """
        create_fiber_section(self.fib_sec)
        self.Initialized = True
        self.UpdateStoredData()


class FibersCircRCCircShape(FibersCirc):
    """
    Class that is the children of FibersCirc and combine the class RCCircShape (section) to retrieve the information needed.  

    @param FibersCirc: Parent class.
    """
    def __init__(self, ID: int, section: RCCircShape, unconf_mat_ID: int, conf_mat_ID: int, bars_mat_ID: int,
        discr_core: list, discr_cover: list, alpha_i=0.0, GJ=0):
        """
        Constructor of the class.

        @param ID (int): Unique fiber section ID.
        @param section (RCCircShape): RCCircShape section object.
        @param unconf_mat_ID (int): ID of material model that will be assigned to the unconfined fibers.
        @param conf_mat_ID (int): ID of material model that will be assigned to the confined fibers.
        @param bars_mat_ID (int): ID of material model that will be assigned to the reinforcing bars fibers.
        @param discr_core (list): List with two entries: number of subdivisions (fibers) in the circumferential direction (number of wedges),
            number of subdivisions (fibers) in the radial direction (number of rings) for the confined core.
        @param discr_cover (list): List with two entries: number of subdivisions (fibers) in the circumferential direction (number of wedges),
            number of subdivisions (fibers) in the radial direction (number of rings) for the unconfined cover.
        @param alpha_i (float, optional): Angle in deg of the first vertical rebars with respect to the y axis, counterclockwise. Defaults to 0.0.
        @param GJ (float, optional): Linear-elastic torsional stiffness assigned to the section. Defaults to 0.0, assume no torsional stiffness.
        """
        self.section = deepcopy(section)
        super().__init__(ID, section.b, section.e, section.D_bars, section.Ay, section.n_bars, section.D_hoops, unconf_mat_ID, conf_mat_ID, bars_mat_ID,
            discr_core, discr_cover, alpha_i=alpha_i, GJ=GJ)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()


class FibersIShape(Fibers):
    """
	Class that stores funcions, material properties, geometric and mechanical parameters for a steel I shape (non double symmetric) fiber section.
    Coordinates: plotting coordinte (x, y) = fiber section coordinate (z, y) = (-x, y). For more information, see the OpenSeesPy documentation.

    @param Fibers: Parent abstract class.
    """
    def __init__(self, ID: int, d, bf_t, bf_b, tf_t, tf_b, tw, top_flange_mat_ID: int, bottom_flange_mat_ID: int, web_mat_ID: int,
        discr_top_flange: list, discr_bottom_flange: list, discr_web: list,  GJ = 0.0):
        """
        Constructor of the class.

        @param ID (int): Unique fiber section ID.
        @param d (float): Depth of the section.
        @param bf_t (float): Top flange's width of the section
        @param bf_b (float): Bottom flange's width of the section
        @param tf_t (float): Top flange's thickness of the section
        @param tf_b (float): Bottom flange's thickness of the section
        @param tw (float): Web's thickness of the section
        @param top_flange_mat_ID (int): ID of material model that will be assigned to the top flange fibers.
        @param bottom_flange_mat_ID (int): ID of material model that will be assigned to the bottom flange fibers.
        @param web_mat_ID (int): ID of material model that will be assigned to the web fibers.
        @param discr_top_flange (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the top flange.
        @param discr_bottom_flange (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the bottom flange.
        @param discr_web (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the web.
        @param GJ (float, optional): Linear-elastic torsional stiffness assigned to the section. Defaults to 0.0, assume no torsional stiffness.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: d needs to be positive.
        @exception NegativeValue: bf_t needs to be positive.
        @exception NegativeValue: bf_b needs to be positive.
        @exception NegativeValue: tf_t needs to be positive.
        @exception NegativeValue: tf_b needs to be positive.
        @exception NegativeValue: tw needs to be positive.
        @exception NegativeValue: top_flange_mat_ID needs to be a positive integer.
        @exception NegativeValue: bottom_flange_mat_ID needs to be a positive integer.
        @exception NegativeValue: web_mat_ID needs to be a positive integer.
        @exception WrongDimension: discr_top_flange has a length of 2.
        @exception WrongDimension: discr_bottom_flange has a length of 2.
        @exception WrongDimension: discr_web has a length of 2.
        @exception NegativeValue: GJ needs to be positive.
        @exception InconsistentGeometry: The sum of the flanges thickness can't be bigger than d.
        """
        # Check
        if ID < 1: raise NegativeValue()
        if d < 0: raise NegativeValue()
        if bf_t < 0: raise NegativeValue()
        if bf_b < 0: raise NegativeValue()
        if tf_b < 0: raise NegativeValue()
        if tf_t < 0: raise NegativeValue()
        if tw < 0: raise NegativeValue()
        if top_flange_mat_ID < 1: raise NegativeValue()
        if bottom_flange_mat_ID < 1: raise NegativeValue()
        if web_mat_ID < 1: raise NegativeValue()
        if len(discr_top_flange) != 2: raise WrongDimension()
        if len(discr_bottom_flange) != 2: raise WrongDimension()
        if len(discr_web) != 2: raise WrongDimension()
        if GJ < 0: raise NegativeValue()
        if tf_t+tf_b >= d: raise InconsistentGeometry()

        # Arguments
        self.ID = ID
        self.d = d
        self.bf_t = bf_t
        self.bf_b = bf_b
        self.tf_t = tf_t
        self.tf_b = tf_b
        self.tw = tw
        self.top_flange_mat_ID = top_flange_mat_ID
        self.bottom_flange_mat_ID = bottom_flange_mat_ID
        self.web_mat_ID = web_mat_ID
        self.discr_top_flange = copy(discr_top_flange)
        self.discr_bottom_flange = copy(discr_bottom_flange)
        self.discr_web = copy(discr_web)
        self.GJ = GJ

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.Initialized = False
        self.ReInit()

    def ReInit(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        # Memebers
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Parameters
        z1 = self.tw/2
        y1 = (self.d - self.tf_t - self.tf_b)/2

        # Create the flange top
        flange_top = [y1, -self.bf_t/2, y1+self.tf_t, self.bf_t/2]
        flange_top_cmd = ['patch', 'rect', self.top_flange_mat_ID, *self.discr_top_flange, *flange_top]

        # Create the flange bottom
        flange_bottom = [-y1-self.tf_b, -self.bf_b/2, -y1, self.bf_b/2]
        flange_bottom_cmd = ['patch', 'rect', self.bottom_flange_mat_ID, *self.discr_bottom_flange, *flange_bottom]

        # Create the web
        web = [-y1, -z1, y1, z1]
        web_cmd = ['patch', 'rect', self.web_mat_ID, *self.discr_web, *web]

        self.fib_sec = [['section', 'Fiber', self.ID, '-GJ', self.GJ], 
            flange_top_cmd, web_cmd, flange_bottom_cmd]

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        self.data = [["INFO_TYPE", "FibersIShape"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["d", self.d],
            ["bf_t", self.bf_t],
            ["bf_b", self.bf_b],
            ["tf_t", self.tf_t],
            ["tf_b", self.tf_b],
            ["tw", self.tw],
            ["GJ", self.GJ],
            ["top_flange_mat_ID", self.top_flange_mat_ID],
            ["bottom_flange_mat_ID", self.bottom_flange_mat_ID],
            ["web_mat_ID", self.web_mat_ID],
            ["discr_top_flange", self.discr_top_flange],
            ["discr_bottom_flange", self.discr_bottom_flange],
            ["discr_web", self.discr_web],
            ["Initialized", self.Initialized]]


    def ShowInfo(self, plot = False, block = False):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the fiber. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
        """
        print("")
        print("Requested info for FibersRect, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print("Depth d = {} mm and web thickness tw = {} mm".format(self.d/mm_unit, self.tw/mm_unit))
        print("Top flange width bf_t = {} mm and thickness tf_t = {} mm".format(self.bf_t/mm_unit, self.tf_t/mm_unit))
        print("Bottom flange width bf_b = {} mm and thickness tf_b = {} mm".format(self.bf_b/mm_unit, self.tf_b/mm_unit))
        print("Web material model ID = {}".format(self.web_mat_ID))
        print("Top flange material model ID = {}".format(self.top_flange_mat_ID))
        print("Bottom flange material model ID = {}".format(self.bottom_flange_mat_ID))
        print("Discretisation in the web [IJ or x/z dir, JK or y dir] = {}".format(self.discr_web))
        print("Discretisation in the top flange [IJ or x/z dir, JK or y dir] = {}".format(self.discr_top_flange))
        print("Discretisation in the bottom flange [IJ or x/z dir, JK or y dir] = {}".format(self.discr_bottom_flange))
        print("")

        if plot:
            plot_fiber_section(self.fib_sec, matcolor=['r', 'b', 'g', 'k'])
            
            if block:
                plt.show()


    def CreateFibers(self):
        """
        Method that initialise the fiber by calling the OpenSeesPy commands.
        """
        create_fiber_section(self.fib_sec)
        self.Initialized = True
        self.UpdateStoredData()


class FibersIShapeSteelIShape(FibersIShape):
    """
    Class that is the children of FibersIShape and combine the class SteelIShape (section) to retrieve the information needed.  

    @param FibersIShape: Parent class.
    """
    def __init__(self, ID: int, section: SteelIShape, top_flange_mat_ID: int, discr_top_flange: list, discr_bottom_flange: list, discr_web: list,
        GJ=0.0, bottom_flange_mat_ID = -1, web_mat_ID = -1):
        """
        Constructor of the class.

        @param ID (int): Unique fiber section ID.
        @param section (SteelIShape): SteelIShape section object.
        @param top_flange_mat_ID (int): ID of material model that will be assigned to the top flange fibers.
        @param discr_top_flange (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the top flange.
        @param discr_bottom_flange (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the bottom flange.
        @param discr_web (list): List with two entries: discretisation in IJ (x/z) and JK (y) for the web.
        @param GJ (float, optional): Linear-elastic torsional stiffness assigned to the section. Defaults to 0.0, assume no torsional stiffness.
        @param bottom_flange_mat_ID (int): ID of material model that will be assigned to the bottom flange fibers.
            Defaults to -1, e.g. equal to top_flange_mat_ID.
        @param web_mat_ID (int): ID of material model that will be assigned to the web fibers.
            Defaults to -1, e.g. equal to top_flange_mat_ID.
        """
        self.section = deepcopy(section)
        if bottom_flange_mat_ID == -1: bottom_flange_mat_ID = top_flange_mat_ID
        if web_mat_ID == -1: web_mat_ID = top_flange_mat_ID

        super().__init__(ID, section.d, section.bf, section.bf, section.tf, section.tf, section.tw, top_flange_mat_ID, bottom_flange_mat_ID, web_mat_ID,
            discr_top_flange, discr_bottom_flange, discr_web, GJ)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()  


def plot_fiber_section(fiber_info, fill_shapes = True, matcolor=['#808080', '#D3D3D3', 'r', 'b', 'g', 'y']):
    """
    Plot fiber cross-section. Coordinate system used: plotting coordinte = (x, y), fiber section coordinate (z, y) = (-x, y)  
    Inspired by plot_fiber_section from ops_vis written by Seweryn Kokot.

    @param fiber_info (list): List of lists (be careful with the local coordinate system!). The first list defines the fiber section: \n
        ['section', 'Fiber', ID, '-GJ', GJ] \n
        The other lists have one of the following format (coordinate input: (y, z)!): \n
        ['layer', 'bar', mat_ID, A, y, z] # one bar \n
        ['layer', 'straight', mat_ID, n_bars, A, yI, zI, yJ, zJ] # line range of bars (with I = first bar, J = last bar) \n
        ['layer', 'circ', mat_ID, n_bars, A, yC, zC, r, (a0_deg), (a1_deg)] # circular range of bars (with C = center, r = radius) \n
        ['patch', 'rect', mat_ID, *discr, -yI, zI, yK, -zK] # rectangle (with yI = yK = d/2; zI = zK = b/2) \n
        ['patch', 'quad', mat_ID, *discr, yI, zI, yJ, zJ, yK, zK, yL, zL] # quadrilateral shaped (starting from bottom left, counterclockwise: I, J, K, L) \n
        ['patch', 'circ', mat_ID, *discr, yC, zC, ri, re, (a0), (a1)] #  (with C = center, ri = internal radius, re = external radius)
    @param fill_shapes (bool, optional): Option to fill fibers with color specified in matcolor. Defaults to True.
    @param matcolor (list, optional): List of colors for various material IDs. Defaults to ['#808080', '#D3D3D3', 'r', 'b', 'g', 'y'].

    Example 1: Simple rectangle with 2 rebars (D = diameter) on top (e distance from the top and from the lateral borders).
        Rectangle with first corner =  I (bottom right) and second corner = K (top left); number of fibers = discr (list of 2)
        fib_sec = [['section', 'Fiber', ID, '-GJ', GJ], 
            ['patch', 'rect', concrete_mat_ID, *discr, yI, zI, yK, zK],
            ['layer', 'bar', bars_mat_ID, Ay, yI-e-D/2, zI-e-D/2], # left rebar
            ['layer', 'bar', bars_mat_ID, Ay, yI-e-D/2, -(zI-e-D/2)]] # right rebar

    Example 2: double symmetric I shape.
        Each rectangle (2 flanges and 1 web): first corner =  I (bottom right) and second corner = K (top left); number of fibers = discr (list of 2)
        fib_sec = [['section', 'Fiber', ID, '-GJ', GJ], 
            ['patch', 'rect', mat_ID, *discr, yI_tf, zI_tf, yK_tf, zK_tf], # top flange
            ['patch', 'rect', mat_ID, *discr, yI_bf, zI_bf, yK_bf, zK_bf], # bottom flange
            ['patch', 'rect', mat_ID, *discr, yI_w, zI_w, yK_w, zK_w]] # web
    """
    
    mat_to_col = {}
    fig, ax = plt.subplots()
    ax.grid(False)

    for item in fiber_info:
        if item[0] == 'section':
            fib_ID = item[2]
        if item[0] == 'layer':
            matID = item[2]
            mat_to_col = __assignColorToMat(matID, mat_to_col, matcolor)
            if item[1] == 'bar':
                As = item[3]
                Iy = item[4]
                Iz = item[5]
                r = np.sqrt(As / np.pi)
                bar = Circle((-Iz, Iy), r, ec='k', fc='k', zorder=10)
                ax.add_patch(bar)
            if item[1] == 'straight':
                n_bars = item[3]
                As = item[4]
                Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                r = np.sqrt(As / np.pi)
                Y = np.linspace(Iy, Jy, n_bars)
                Z = np.linspace(Iz, Jz, n_bars)
                for zi, yi in zip(Z, Y):
                    bar = Circle((-zi, yi), r, ec='k', fc=mat_to_col[matID], zorder=10)
                    ax.add_patch(bar)
            if item[1] == 'circ':
                n_bars, As = item[3], item[4]
                yC, zC, r = item[5], item[6], item[7]
                if len(item) > 8:
                    a0_deg = item[8]
                else:
                    a0_deg = 0.
                a1_deg = 360. - 360./n_bars + a0_deg
                if len(item) > 9: a1_deg = item[9]

                a0_rad, a1_rad = np.pi * a0_deg / 180., np.pi * a1_deg / 180.
                r_bar = np.sqrt(As / np.pi)
                thetas = np.linspace(a0_rad, a1_rad, n_bars)
                Y = yC + r * np.cos(thetas)
                Z = zC + r * np.sin(thetas)
                for zi, yi in zip(Z, Y):
                    bar = Circle((-zi, yi), r_bar, ec='k', fc=mat_to_col[matID], zorder=10)
                    ax.add_patch(bar)

        if (item[0] == 'patch' and (item[1] == 'quad' or item[1] == 'quadr' or
                                  item[1] == 'rect')):
            matID, nIJ, nJK = item[2], item[3], item[4]
            mat_to_col = __assignColorToMat(matID, mat_to_col, matcolor)


            if item[1] == 'quad' or item[1] == 'quadr':
                Iy, Iz, Jy, Jz = item[5], item[6], item[7], item[8]
                Ky, Kz, Ly, Lz = item[9], item[10], item[11], item[12]

            if item[1] == 'rect':
                Iy, Iz, Ky, Kz = item[5], item[6], item[7], item[8]
                Jy, Jz, Ly, Lz = Iy, Kz, Ky, Iz
                # check order of definition
                if Kz-Iz < 0 or Ky-Iy < 0: print("!!!!!!! WARNING !!!!!!! The fiber is not defined bottom right, top left")
            
            # check for convexity (vector products)
            outIJxIK = (Jy-Iy)*(Kz-Iz) - (Ky-Iy)*(Jz-Iz)
            outIKxIL = (Ky-Iy)*(Lz-Iz) - (Ly-Iy)*(Kz-Iz)
            # check if I, J, L points are colinear
            outIJxIL = (Jy-Iy)*(Lz-Iz) - (Ly-Iy)*(Jz-Iz)
            # outJKxJL = (Ky-Jy)*(Lz-Jz) - (Ly-Jy)*(Kz-Jz)

            if -outIJxIK <= 0 or -outIKxIL <= 0 or -outIJxIL <= 0:
                print('!!!!!!! WARNING !!!!!!! Patch quad is non-convex or non-counter-clockwise defined or has at least 3 colinear points in line')

            IJz, IJy = np.linspace(Iz, Jz, nIJ+1), np.linspace(Iy, Jy, nIJ+1)
            JKz, JKy = np.linspace(Jz, Kz, nJK+1), np.linspace(Jy, Ky, nJK+1)
            LKz, LKy = np.linspace(Lz, Kz, nIJ+1), np.linspace(Ly, Ky, nIJ+1)
            ILz, ILy = np.linspace(Iz, Lz, nJK+1), np.linspace(Iy, Ly, nJK+1)

            if fill_shapes:
                Z = np.zeros((nIJ+1, nJK+1))
                Y = np.zeros((nIJ+1, nJK+1))

                for j in range(nIJ+1):
                    Z[j, :] = np.linspace(IJz[j], LKz[j], nJK+1)
                    Y[j, :] = np.linspace(IJy[j], LKy[j], nJK+1)

                for j in range(nIJ):
                    for k in range(nJK):
                        zy = np.array([[-Z[j, k], Y[j, k]],
                                       [-Z[j, k+1], Y[j, k+1]],
                                       [-Z[j+1, k+1], Y[j+1, k+1]],
                                       [-Z[j+1, k], Y[j+1, k]]])
                        poly = Polygon(zy, True, ec='k', fc=mat_to_col[matID])
                        ax.add_patch(poly)

            else:
                # horizontal lines
                for az, bz, ay, by in zip(IJz, LKz, IJy, LKy):
                    plt.plot([-az, -bz], [ay, by], 'b-', zorder=1)

                # vertical lines
                for az, bz, ay, by in zip(JKz, ILz, JKy, ILy):
                    plt.plot([-az, -bz], [ay, by], 'b-', zorder=1)

        if item[0] == 'patch' and item[1] == 'circ':
            matID, nc, nr = item[2], item[3], item[4]
            mat_to_col = __assignColorToMat(matID, mat_to_col, matcolor)

            yC, zC, ri, re = item[5], item[6], item[7], item[8]
            if len(item) > 9:
                a0 = item[9]
            else:
                a0 = 0.
            a1 = 360. + a0
            if len(item) > 10: a1 = item[10]

            dr = (re - ri) / nr
            dth = (a1 - a0) / nc

            for j in range(nr):
                rj = ri + j * dr
                rj1 = rj + dr

                for i in range(nc):
                    thi = a0 + i * dth
                    thi1 = thi + dth
                    wedge = Wedge((-zC, yC), rj1, thi, thi1, width=dr, ec='k', #Seweryn Kokot: (yC, -zC), wrong??
                        lw=1, fc=mat_to_col[matID])
                    ax.add_patch(wedge)
    ax.set(xlabel='x dimension [{}]'.format(length_unit), ylabel='y dimension [{}]'.format(length_unit), 
                    title='Fiber section (ID = {})'.format(fib_ID))
    ax.axis('equal')


def __assignColorToMat(matID: int, mat_to_col: dict, matcolor: list):
    """
    PRIVATE FUNCTION. Used to assign different colors for each material model assign to the fiber section.

    @param matID (int): ID of the material model.
    @param mat_to_col (dict): Dictionary to check with material model has which color.
    @param matcolor (list): List of colors.

    @returns dict: Updated dictionary.
    """
    if not matID in mat_to_col:
        if len(mat_to_col) >= len(matcolor):
            print("Warning: not enough colors defined for fiber section plot (white used)")
            mat_to_col[matID] = 'w'
        else:
            mat_to_col[matID] = matcolor[len(mat_to_col)]
    return mat_to_col


def create_fiber_section(fiber_info):
    """
    Initialise fiber cross-section with OpenSeesPy commands.
    For examples, see plot_fiber_section.
    Inspired by fib_sec_list_to_cmds from ops_vis written by Seweryn Kokot 

    @param fiber_info (list): List of lists (be careful with the local coordinate system!). The first list defines the fiber section: \n
        ['section', 'Fiber', ID, '-GJ', GJ] \n
        The other lists have one of the following format (coordinate input: (y, z)!): \n
        ['layer', 'bar', mat_ID, A, y, z] # one bar \n
        ['layer', 'straight', mat_ID, n_bars, A, yI, zI, yJ, zJ] # line range of bars (with I = first bar, J = last bar) \n
        ['layer', 'circ', mat_ID, n_bars, A, yC, zC, r, (a0_deg), (a1_deg)] # circular range of bars (with C = center, r = radius) \n
        ['patch', 'rect', mat_ID, *discr, -yI, zI, yK, -zK] # rectangle (with yI = yK = d/2; zI = zK = b/2) \n
        ['patch', 'quad', mat_ID, *discr, yI, zI, yJ, zJ, yK, zK, yL, zL] # quadrilateral shaped (starting from bottom left, counterclockwise: I, J, K, L) \n
        ['patch', 'circ', mat_ID, *discr, yC, zC, ri, re, (a0), (a1)] #  (with C = center, ri = internal radius, re = external radius)
    """
    for dat in fiber_info:
        if dat[0] == 'section':
            fib_ID, GJ = dat[2], dat[4]
            section('Fiber', fib_ID, 'GJ', GJ)

        if dat[0] == 'layer':
            mat_ID = dat[2]
            if dat[1] == 'straight':
                n_bars = dat[3]
                As = dat[4]
                Iy, Iz, Jy, Jz = dat[5], dat[6], dat[7], dat[8]
                layer('straight', mat_ID, n_bars, As, Iy, Iz, Jy, Jz)
            if dat[1] == 'bar':
                As = dat[3]
                Iy = dat[4]
                Iz = dat[5]
                fiber(Iy, Iz, As, mat_ID)
                # layer('straight', mat_ID, 1, As, Iy, Iz, Iy, Iz)
            if dat[1] == 'circ':
                n_bars, As = dat[3], dat[4]
                yC, zC, r = dat[5], dat[6], dat[7]
                if len(dat) > 8:
                    a0_deg = dat[8]
                else:
                    a0_deg = 0.
                a1_deg = 360. - 360./n_bars + a0_deg
                if len(dat) > 9: a1_deg = dat[9]
                layer('circ', mat_ID, n_bars, As, yC, zC, r, a0_deg, a1_deg)

        if dat[0] == 'patch':
            mat_ID = dat[2]
            nIJ = dat[4] ###
            nJK = dat[3] ###

            if dat[1] == 'quad' or dat[1] == 'quadr':
                Iy, Iz, Jy, Jz = dat[5], dat[6], dat[7], dat[8]
                Ky, Kz, Ly, Lz = dat[9], dat[10], dat[11], dat[12]
                patch('quad', mat_ID, nIJ, nJK, Iy, Iz, Jy, Jz, Ky, Kz,
                        Ly, Lz)

            if dat[1] == 'rect':
                Iy, Iz, Ky, Kz = dat[5], dat[6], dat[7], dat[8]
                patch('rect', mat_ID, nIJ, nJK, Iy, Iz, Ky, Kz)
                # patch('rect', mat_ID, nIJ, nJK, Iy, Kz, Ky, Iz)

            if dat[1] == 'circ':
                mat_ID, nc, nr = dat[2], dat[3], dat[4]
                yC, zC, ri, re = dat[5], dat[6], dat[7], dat[8]
                if len(dat) > 9:
                    a0 = dat[9]
                else:
                    a0 = 0.
                a1 = 360. + a0
                if len(dat) > 10: a1 = dat[10]
                patch('circ', mat_ID, nc, nr, yC, zC, ri, re, a0, a1)
        
