# Module for the Section classes
#          Carmine Schipani, 2021

import numpy as np
import math
from DataManagement import *
from ErrorHandling import *
from Units import *

class Section(DataManagement):
    """Parent abstract class that groups every class in this module.

    Args:
        DataManagement: Parent abstract class
    """
    pass


class SteelIShape(Section):
    """Class that stores funcions, geometric and mechanical properties of a steel double symmetric I-shape profile.

    Args:
        Section: Parent abstract class
    """

    global n
    n = 10.0

    def __init__(self, Type, d, bf, tf, tw, L, r, E, Fy, Fy_web = -1, name_tag = "Not Defined"):
        """The conctructor of the class

        Args:
            Type (str): The type of element (Col or Beam)
            d (double): The depth
            bf (double): The flange's width
            tf (double): The flange's thickness
            tw (double): The web's thickness
            L (double): The element length (exclude the panel zone if present)
            r (double): The radius of the weld fillets
            E (double): The Young modulus
            Fy (double): The yield strength
            Fy_web (double, optional): The yield strength of the web. Defaults to -1 (e.g. computed in __init__()).
            name_tag (str, optional): A nametag for the section. Defaults to "Not Defined".
        """
        # Check
        if Type != "Beam" and Type != "Col": raise WrongArgument()
        if d < 0: raise NegativeValue()
        if bf < 0: raise NegativeValue()
        if tf < 0: raise NegativeValue()
        if tw < 0: raise NegativeValue()
        if L < 0: raise NegativeValue()
        if r < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if Fy < 0: raise NegativeValue()
        if Fy_web != -1 and Fy_web < 0: raise NegativeValue()
        if tw > bf: raise InconsistentGeometry()
        if tf > d/2: raise InconsistentGeometry()
        if r > bf/2: raise InconsistentGeometry()

        # Arguments
        self.Type = Type
        self.d = d
        self.bf = bf
        self.tf = tf
        self.tw = tw
        self.L = L
        self.r = r
        self.E = E
        self.Fy = Fy
        self.Fy_web = Fy_web
        if self.Fy_web == -1: self.Fy_web = Fy
        self.name_tag = name_tag

        # Initialized the parameters that are dependent from others
        self.ReInit()

    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "copy()" from the module "copy".
        """
        # Member
        self.h_1 = self.d - 2.0*self.r -2.0*self.tf
        self.A = self.ComputeA()
        self.Npl = self.A*self.Fy
        self.Iy = self.ComputeIy()
        self.Iz = self.ComputeIz()
        self.Wply = self.ComputeWply()
        self.Wplz = self.ComputeWplz()
        self.My = self.Fy*self.Wply
        self.Iy_mod = self.Iy*(n + 1.0)/n
        self.iz = self.Compute_iz()
        self.iy = self.Compute_iy()

        # Data storage for loading/saving
        self.data = ["SteelIShape", # Tag for differentiating different data
            self.name_tag, 
            self.Type, 
            self.d, 
            self.bf, 
            self.tf, 
            self.tw, 
            self.L, 
            self.r, 
            self.h_1, 
            self.E, 
            self.Fy, 
            self.Fy_web, 
            self.A, 
            self.Iy, 
            self.Iz, 
            self.Wply, 
            self.Wplz, 
            self.Iy_mod, 
            self.iy,
            self.iz,
            self.Npl, 
            self.My] 

    def ShowInfo(self):
        """Function that show the data stored in the class in the command window.
        """
        print("")
        print("Requested info for steel I shape section of type = {} and name tag = {}".format(self.Type, self.name_tag))
        print("d = {} mm".format(self.d/mm_unit))
        print("Fy = {} MPa".format(self.Fy/MPa_unit))
        print("Fy web = {} MPa".format(self.Fy_web/MPa_unit))
        print("E = {} GPa".format(self.E/GPa_unit))
        print("h_1 = {} mm".format(self.h_1/mm_unit))
        print("A = {} mm2".format(self.A/mm2_unit))
        print("Iy = {} mm4".format(self.Iy/mm4_unit))
        print("Iz = {} mm4".format(self.Iz/mm4_unit))
        print("Wply = {} mm3".format(self.Wply/mm3_unit))
        print("Wplz = {} mm3".format(self.Wplz/mm3_unit))
        print("Iy_mod = {} mm4".format(self.Iy_mod/mm4_unit))
        print("iy = {} mm".format(self.iy/mm_unit))
        print("iz = {} mm".format(self.iz/mm_unit))
        print("My = {} kNm".format(self.My/kNm_unit))
        print("Npl = {} kN".format(self.Npl/kN_unit))


    def ComputeA(self):
        """Compute the area of a double symmetric I-profile section with fillets.

        Returns:
            double: The area of the section
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        # without fillets bf*tf*2 + tw*(d-2*tf)
        return 2.0*self.bf*self.tf+self.tw*(self.d-2.0*self.tf)+0.8584*self.r**2

    def ComputeIy(self):
        """Compute the moment of inertia of a double symmetric I-profile section, with respect to its strong axis with fillets.

        Returns:
            double: The moment of inertia with respect to the strong axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        # without fillets: bf*tf/2*(d-tf)**2 + bf*tf**3/6 + (d-tf*2)**3*tf/12 
        return (self.bf*self.d**3.0-(self.bf-self.tw)*(self.d-2.0*self.tf)**3)/12.0+0.8584*self.r**2*(0.5*self.d-self.tf-0.4467*self.r/2.0)**2

    def ComputeIz(self):
        """Compute the moment of inertia of a double symmetric I-profile section, with respect to its weak axis with fillets.

        Returns:
            double: The moment of inertia with respect to the weak axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        return (self.tf*self.bf**3)/6.0+((self.d-2.0*self.tf)*self.tw**3)/12.0+0.8584*self.r**2*(0.5*self.tw+0.2234*self.r)**2

    def ComputeWply(self):
        """Compute the plastic modulus of a double symmetric I-profile section, with respect to its strong axis with fillets.

        Returns:
            double: The plastic modulus with respect to the strong axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        return self.bf*self.tf*(self.d-self.tf)+(self.d-2.0*self.tf)**2.0*(self.tw/4.0)+0.4292*self.r**2*(self.d-2.0*self.tf-0.4467*self.r)

    def ComputeWplz(self):
        """Compute the plastic modulus of a double symmetric I-profile section, with respect to its weak axis with fillets.

        Returns:
            double: The plastic modulus with respect to the weak axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius
        #TODO: add correct formula (at the moment is the formula for Wply)
        return self.bf*self.tf*(self.d-self.tf)+(self.d-2.0*self.tf)**2.0*(self.tw/4.0)+0.4292*self.r**2*(self.d-2.0*self.tf-0.4467*self.r)

    def Compute_iy(self):
        """Compute the gyration radius with respect to the strong axis.

        Returns:
            double: The gyration radius with respect to the strong axis
        """
        #  Iy :    The second moment of inertia with respect to thte strong axis
        #  A :     The area

        return math.sqrt(self.Iy/self.A)

    def Compute_iz(self):
        """Compute the gyration radius with respect to the weak axis.

        Returns:
            double: The gyration radius with respect to the weak axis
        """
        #  Iz :    The second moment of inertia with respect to thte weak axis
        #  A :     The area

        return math.sqrt(self.Iz/self.A)


class RCRectShape(Section):
    """Class that stores funcions, geometric and mechanical properties of RC rectangular shape profile.

    Args:
        Section: Parent abstract class
    """

    def __init__(self, b, d, L, e, fc, D_bars, bars_position_x: np.ndarray, bars_ranges_position_y: np.ndarray, fy, Ey, 
        D_hoops, s, fs, Es, name_tag = "Not Defined", rho_s_x = -1, rho_s_y = -1, Ec = -1):
        """The conctructor of the class.

        Args:
            b (double): Width of the section
            d (double): Depth of the section
            L (double): Length of the element
            e (double): Concrete cover
            fc (double): Unconfined concrete compressive strength (cylinder test)
            D_bars (double): Diameter of the reinforcing bars
            bars_position_x (np.ndarray): Distances from border to bar centerline, bar to bar centerlines and finally bar centerline to border in the x direction (aligned).
                Starting from the left to right, from the top range to the bottom one. The number of bars for each range can vary; in this case, add this argument when defining the array " dtype = object"
            bars_ranges_position_y (np.ndarray): Distances from border to range centerlines, range to range centerlines and finally range centerline to border in the y direction.
                Starting from the top range to the bottom one
            fy (double): Yield stress for reinforcing bars
            Ey (double): Young modulus for reinforcing bars
            D_hoops (double): Diameter of the hoops
            s (double): Centerline distance for the hoops
            fs (double): Yield stress for the hoops
            Es (double): Young modulus for the hoops
            name_tag (str, optional): A nametag for the section. Defaults to "Not Defined".
            rho_s_x (double, optional): Volumetric or not?????????????. Defaults to -1 (e.g. computed in __init__() assuming one range of hoops).
            rho_s_y (double, optional): Volumetric or not?????????????. Defaults to -1 (e.g. computed in __init__() assuming one range of hoops).
            Ec (double, optional): Young modulus for concrete. Defaults to -1 (e.g. computed in __init__()).
        """
        # Note that for the validity of the formulas, at least one bar per corner and at least one hoop closed with 135 degress!

        #TODO: ask someone good with concrete (Katrin Beyer, Diego,... ) if rho_s is volumtric or not (and how to compute it)
        
        # Check
        if b < 0: raise NegativeValue()
        if d < 0: raise NegativeValue()
        if L < 0: raise NegativeValue()
        if e < 0: raise NegativeValue()
        if fc > 0: raise PositiveValue()
        if D_bars < 0: raise NegativeValue()
        if fy < 0: raise NegativeValue()
        if Ey < 0: raise NegativeValue()
        if D_hoops < 0: raise NegativeValue()
        if s < 0: raise NegativeValue()
        if fs < 0: raise NegativeValue()
        if Es < 0: raise NegativeValue()
        if rho_s_x != -1 and rho_s_x < 0: raise NegativeValue()
        if rho_s_y != -1 and rho_s_y < 0: raise NegativeValue()
        if Ec != -1 and Ec < 0: raise NegativeValue()
        if np.size(bars_position_x) != np.size(bars_ranges_position_y)-1: raise WrongDimension()
        geometry_tol = 5*mm_unit
        for bars in bars_position_x:
            if abs(np.sum(bars) - b) > geometry_tol: raise InconsistentGeometry()
        if abs(np.sum(bars_ranges_position_y)-d) > geometry_tol: raise InconsistentGeometry()
        if e > b/2 or e > d/2: raise InconsistentGeometry()

        # Arguments
        self.b = b
        self.d = d
        self.L = L
        self.e = e
        self.fc = fc
        self.D_bars = D_bars
        self.bars_position_x = bars_position_x
        self.bars_ranges_position_y = bars_ranges_position_y
        self.fy = fy
        self.Ey = Ey
        self.D_hoops = D_hoops
        self.s = s
        self.fs = fs
        self.Es = Es
        self.name_tag = name_tag

        # Initialized the parameters that are dependent from others
        self.ReInit(rho_s_x, rho_s_y, Ec)


    def ReInit(self, rho_s_x = -1, rho_s_y = -1, Ec = -1):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "copy()" from the module "copy".
        """
        # Precompute some members
        self.cl_hoops = self.e + self.D_hoops/2.0 # centerline distance from the border of the extreme confining hoops
        self.cl_bars = self.e + self.D_bars/2.0 + self.D_hoops # centerline distance from the border of the corner bars
        self.bc = self.b - self.cl_hoops
        self.dc = self.d - self.cl_hoops
        self.As = ComputeACircle(self.D_hoops)

        # Arguments
        self.rho_s_x = rho_s_x
        if self.rho_s_x == -1: self.rho_s_x = 2.0*ComputeRho(self.As, 1, self.bc*self.s)
        self.rho_s_y = rho_s_y
        if self.rho_s_y == -1: self.rho_s_y = 2.0*ComputeRho(self.As, 1, self.dc*self.s)
        self.Ec = Ec
        if self.Ec == -1: self.Ec = self.ComputeEc()

        # Members
        self.nr_bars = self.ComputeNrBars()
        self.A = self.ComputeA()
        self.Ac = self.ComputeAc()
        self.Ay = ComputeACircle(self.D_bars)
        self.rho_bars = ComputeRho(self.Ay, self.nr_bars, self.A)
        self.Iy = self.ComputeIy()
        self.Iz = self.ComputeIz()

        # Data storage for loading/saving
        self.data = ["RCRectShape", # Tag for differentiating different data
            self.name_tag,
            self.b,
            self.d,
            self.bc,
            self.dc,
            self.L,
            self.e,
            self.A,
            self.Ac,
            self.Iy,
            self.Iz,
            self.fc,
            self.Ec,
            self.D_bars,
            self.nr_bars,
            self.Ay,
            self.bars_position_x,
            self.bars_ranges_position_y,
            self.rho_bars,
            self.cl_bars,
            self.fy,
            self.Ey,
            self.D_hoops,
            self.s,
            self.As,
            self.rho_s_x,
            self.rho_s_y,
            self.cl_hoops,
            self.fs,
            self.Es]


    def ShowInfo(self):
        """Function that show the data stored in the class in the command window.
        """
        print("")
        print("Requested info for RC rectangular section of name tag = {}".format(self.name_tag))
        print("b = {} mm".format(self.b/mm_unit))
        print("d = {} mm".format(self.d/mm_unit))
        print("e = {} mm".format(self.e/mm_unit))
        print("A = {} mm2".format(self.A/mm2_unit))
        print("Ac = {} mm2".format(self.Ac/mm2_unit))
        print("fc = {} MPa".format(self.fc/MPa_unit))
        print("Ec = {} GPa".format(self.Ec/GPa_unit))
        print("D bars = {} and Ay = {} mm2 with {} bars".format(self.D_bars/mm_unit, self.Ay/mm2_unit, self.nr_bars))
        print("D hoops = {} and As = {} mm2".format(self.D_hoops/mm_unit, self.As/mm2_unit))
        print("rho bars = {} ".format(self.rho_bars))
        print("rho hoops in x = {} ".format(self.rho_s_x))
        print("rho hoops in y = {} ".format(self.rho_s_y))
        print("Iy = {} mm4".format(self.Iy/mm4_unit))
        print("Iz = {} mm4".format(self.Iz/mm4_unit))
        print("")
    

    def ComputeNrBars(self):
        """Compute the number of vertical bars in the array bars_position_x (note that this list of lists can have different list sizes).

        Returns:
            double: Number of vertical bars
        """
        nr_bars = 0
        for range in self.bars_position_x:
            nr_bars += np.size(range)-1

        return nr_bars


    def ComputeEc(self):
        """Compute Ec using the formula from Mander et Al. 1988.

        Returns:
            double: Young modulus of concrete
        """

        return 5000.0 * math.sqrt(-self.fc/MPa_unit) * MPa_unit


    def ComputeA(self):
        """Compute the area for a rectangular section.

        Returns:
            double: Total area
        """
        return self.b * self.d


    def ComputeAc(self):
        """Compute the confined area.

        Returns:
            double: Confined area
        """
        return self.bc * self.dc


    def ComputeIy(self):
        """Compute the moment of inertia of the rectangular section with respect to the strong axis.

        Returns:
            double: Moment of inertia (strong axis)
        """
        return self.b * self.d**3 / 12.0


    def ComputeIz(self):
        """Compute the moment of inertia of the rectangular section with respect to the weak axis.

        Returns:
            double: Moment of inertia (weak axis)
        """
        return self.d * self.b**3 / 12.0


def ComputeACircle(D):
    """Compute the area of one circle (reinforcing bar or hoop).

    Args:
        D (double): Diameter of the circle (reinforcing bar of hoop)

    Returns:
        double: Area the circle (for reinforcing bars or hoops)
    """
    return D**2/4.0*math.pi


def ComputeRho(A, nr, A_tot):
    """Compute the ratio of area of a reinforcement to area of a section.

    Args:
        A (double): Area of reinforcement
        nr (int): Number of reinforcement
        A_tot (double): Area of the concrete

    Returns:
        double: Ratio
    """
    return nr * A / A_tot


