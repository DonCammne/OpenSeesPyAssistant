"""
Module for the section (steel I shape profiles, RC circular/square/rectangular sections).
Carmine Schipani, 2021
"""

import numpy as np
import math
from copy import copy, deepcopy
from OpenSeesPyAssistant.DataManagement import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *

class Section(DataManagement):
    """
    Parent abstract class for the storage and manipulation of a section's information (mechanical and geometrical parameters, etc).

    @param DataManagement: Parent abstract class.
    """
    pass


class SteelIShape(Section):
    """
    Class that stores funcions, geometric and mechanical properties of a steel double symmetric I-shape profile.
    The parameter 'n' is used as global throughout the SteelIShape sections to optimise the program (given the fact that is constant everytime).

    @param Section: Parent abstract class.
    """
    global n
    n = 10.0

    def __init__(self, Type: str, d, bf, tf, tw, L, r, E, Fy, Fy_web = -1, name_tag = "Not Defined"):
        """
        The conctructor of the class.

        @param Type (str): Type of the section. It can be 'Col' for column or 'Beam' for beams.
        @param d (float): Depth of the section.
        @param bf (float): Flange's width of the section
        @param tf (float): Flange's thickness of the section
        @param tw (float): Web's thickness of the section
        @param L (float): Effective length of the element associated with this section.
            If the panel zone is present, exclude its dimension.
        @param r (float): Radius of the weld fillets of the section.
        @param E (float): Young modulus of the section.
        @param Fy (float): Yield strength of the flange of the section. Used as the yield strength of the entire section.
        @param Fy_web (float, optional): Yield strength of the web of the section. Used for panel zone associated to this section.
            Defaults to -1, e.g. computed in __init__() as equal to Fy.
        @param name_tag (str, optional): Name TAG of the section. Defaults to "Not Defined".

        @exception WrongArgument: Type needs to be 'Col' or 'Beam'.
        @exception NegativeValue: d needs to be positive.
        @exception NegativeValue: bf needs to be positive.
        @exception NegativeValue: tf needs to be positive.
        @exception NegativeValue: tw needs to be positive.
        @exception NegativeValue: L needs to be positive.
        @exception NegativeValue: r needs to be positive.
        @exception NegativeValue: E needs to be positive.
        @exception NegativeValue: Fy needs to be positive.
        @exception NegativeValue: Fy_web needs to be positive if different from -1.
        @exception InconsistentGeometry: tw should be smaller than bf.
        @exception InconsistentGeometry: tf needs to be smaller than half of d
        @exception InconsistentGeometry: r should be less than half bf and d
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
        if r > bf/2 or r > d/2: raise InconsistentGeometry()

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
        self.Fy_web = Fy if Fy_web == -1 else Fy_web
        self.name_tag = name_tag

        # Initialized the parameters that are dependent from others
        self.ReInit()

    def ReInit(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
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
        self.UpdateStoredData()

    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        self.data = [["INFO_TYPE", "SteelIShape"], # Tag for differentiating different data
            ["name_tag", self.name_tag], 
            ["Type", self.Type], 
            ["d", self.d], 
            ["bf", self.bf], 
            ["tf", self.tf], 
            ["tw", self.tw], 
            ["L", self.L], 
            ["r", self.r], 
            ["h_1", self.h_1], 
            ["E", self.E], 
            ["Fy", self.Fy], 
            ["Fy_web", self.Fy_web], 
            ["A", self.A], 
            ["Iy", self.Iy], 
            ["Iz", self.Iz], 
            ["Wply", self.Wply], 
            ["Wplz", self.Wplz], 
            ["Iy_mod", self.Iy_mod], 
            ["iy", self.iy],
            ["iz", self.iz],
            ["Npl", self.Npl], 
            ["My", self.My]] 

    def ShowInfo(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
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
        print("")


    def ComputeA(self):
        """
        Compute the area of a double symmetric I-profile section with fillets.

        @returns float: Area of the I shape section (with fillets included)
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        # without fillets bf*tf*2 + tw*(d-2*tf)
        return 2.0*self.bf*self.tf+self.tw*(self.d-2.0*self.tf)+0.8584*self.r**2

    def ComputeIy(self):
        """
        Compute the moment of inertia of a double symmetric I-profile section, with respect to its strong axis with fillets.

        @returns float: The moment of inertia with respect to the strong axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        # without fillets: bf*tf/2*(d-tf)**2 + bf*tf**3/6 + (d-tf*2)**3*tf/12 
        return (self.bf*self.d**3.0-(self.bf-self.tw)*(self.d-2.0*self.tf)**3)/12.0+0.8584*self.r**2*(0.5*self.d-self.tf-0.4467*self.r/2.0)**2

    def ComputeIz(self):
        """
        Compute the moment of inertia of a double symmetric I-profile section, with respect to its weak axis with fillets.

        @returns float: The moment of inertia with respect to the weak axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        return (self.tf*self.bf**3)/6.0+((self.d-2.0*self.tf)*self.tw**3)/12.0+0.8584*self.r**2*(0.5*self.tw+0.2234*self.r)**2

    def ComputeWply(self):
        """
        Compute the plastic modulus of a double symmetric I-profile section, with respect to its strong axis with fillets.

        @returns float: The plastic modulus with respect to the strong axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius

        return self.bf*self.tf*(self.d-self.tf)+(self.d-2.0*self.tf)**2.0*(self.tw/4.0)+0.4292*self.r**2*(self.d-2.0*self.tf-0.4467*self.r)

    def ComputeWplz(self):
        """
        Compute the plastic modulus of a double symmetric I-profile section, with respect to its weak axis with fillets.

        @returns float: The plastic modulus with respect to the weak axis.
        """
        #  d :     The depth
        #  bf :    The flange's width
        #  tf :    The flange's thickness
        #  tw :    The web's thickness
        #  r :     The weld fillet radius
        return (self.tf*self.bf**2)/2+(self.d-2.0*self.tf)*(self.tw**2/4.0)+0.4292*self.r**2*(self.tw+0.4467*self.r)

    def Compute_iy(self):
        """
        Compute the gyration radius with respect to the strong axis.

        @returns float: The gyration radius with respect to the strong axis.
        """
        #  Iy :    The second moment of inertia with respect to thte strong axis
        #  A :     The area

        return math.sqrt(self.Iy/self.A)

    def Compute_iz(self):
        """
        Compute the gyration radius with respect to the weak axis.

        @returns float: The gyration radius with respect to the weak axis.
        """
        #  Iz :    The second moment of inertia with respect to thte weak axis
        #  A :     The area

        return math.sqrt(self.Iz/self.A)


class RCRectShape(Section):
    """
    Class that stores funcions, geometric and mechanical properties of RC rectangular shape profile.
    Note that for the validity of the formulas, at least one bar per corner and at least one hoop closed (with 135 degress possibly).

    @param Section: Parent abstract class.
    """
    def __init__(self, b, d, L, e, fc, D_bars, bars_position_x: np.ndarray, bars_ranges_position_y: np.ndarray, fy, Ey, 
        D_hoops, s, fs, Es, name_tag = "Not Defined", rho_s_x = -1, rho_s_y = -1, Ec = -1):
        """
        The conctructor of the class.

        @param b (float): Width of the section.
        @param d (float): Depth of the section.
        @param L (float): Effective length of the element associated with this section.
            If the panel zone is present, exclude its dimension.
        @param e (float): Concrete cover.
        @param fc (float): Unconfined concrete compressive strength (cylinder test).
        @param D_bars (float): Diameter of the reinforcing bars.
        @param bars_position_x (np.ndarray): Array with a range of aligned vertical reinforcing bars for each row in x direction.
            Distances from border to bar centerline, bar to bar centerlines and
            finally bar centerline to border in the x direction (aligned).
            Starting from the left to right, from the top range to the bottom one.
            The number of bars for each range can vary; in this case, add this argument when defining the array " dtype = object".
        @param bars_ranges_position_y (np.ndarray): Array of dimension 1 with the position or spacing in y of the ranges in bars_position_x.
            Distances from border to range centerlines, range to range centerlines and
            finally range centerline to border in the y direction.
            Starting from the top range to the bottom one.
        @param fy (float): Yield stress for reinforcing bars.
        @param Ey (float): Young modulus for reinforcing bars.
        @param D_hoops (float): Diameter of the hoops.
        @param s (float): Centerline distance for the hoops.
        @param fs (float): Yield stress for the hoops.
        @param Es (float): Young modulus for the hoops
        @param name_tag (str, optional): A nametag for the section. Defaults to "Not Defined".
        @param rho_s_x (float, optional): Ratio of the transversal area of the hoops to the associated concrete area in the x direction.
            Defaults to -1, e.g. computed in __init__()  and ReInit() assuming one range of hoops.
        @param rho_s_y (float, optional): Ratio of the transversal area of the hoops to the associated concrete area in the y direction.
            Defaults to -1, e.g. computed in __init__()  and ReInit() assuming one range of hoops.
        @param Ec (float, optional): Young modulus for concrete. Defaults to -1, e.g. computed in __init__() and ReInit().

        @exception NegativeValue: b needs to be positive.
        @exception NegativeValue: d needs to be positive.
        @exception NegativeValue: L needs to be positive.
        @exception NegativeValue: e needs to be positive.
        @exception PositiveValue: fc needs to be negative.
        @exception NegativeValue: D_bars needs to be positive.
        @exception NegativeValue: fy needs to be positive.
        @exception NegativeValue: Ey needs to be positive.
        @exception NegativeValue: D_hoops needs to be positive.
        @exception NegativeValue: s needs to be positive.
        @exception NegativeValue: fs needs to be positive.
        @exception NegativeValue: Es needs to be positive.
        @exception NegativeValue: rho_s_x needs to be positive if different from -1. 
        @exception NegativeValue: rho_s_y needs to be positive if different from -1.
        @exception NegativeValue: Ec needs to be positive if different from -1.
        @exception WrongDimension: Number of lists in the list bars_position_x needs to be the same of the length of bars_ranges_position_y - 1.
        @exception InconsistentGeometry: The sum of the distances for each list in bars_position_x should be equal to the section's width (tol = 5 mm).
        @exception InconsistentGeometry: The sum of the distances in bars_ranges_position_y should be equal to the section's depth (tol = 5 mm).
        @exception InconsistentGeometry: e should be smaller than half the depth and the width of the section.
        """
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
        warning_min_bars = "!!!!!!! WARNING !!!!!!! The hypothesis of one bar per corner (aligned) is not fullfilled."
        if len(bars_position_x) < 2:
            print(warning_min_bars)
        elif len(bars_position_x[0]) < 3 or len(bars_position_x[-1]) < 3:
            print(warning_min_bars)

        # Arguments
        self.b = b
        self.d = d
        self.L = L
        self.e = e
        self.fc = fc
        self.D_bars = D_bars
        self.bars_position_x = deepcopy(bars_position_x)
        self.bars_ranges_position_y = copy(bars_ranges_position_y)
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
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.

        @param rho_s_x (float, optional): Ratio of the transversal area of the hoops to the associated concrete area in the x direction.
            Defaults to -1, e.g. computed assuming one range of hoops.
        @param rho_s_y (float, optional): Ratio of the transversal area of the hoops to the associated concrete area in the y direction.
            Defaults to -1, e.g. computed assuming one range of hoops.
        @param Ec (float, optional): Young modulus for concrete. Defaults to -1, e.g. computed according to Mander et Al. 1988.
        """
        # Precompute some members
        self.cl_hoops = self.e + self.D_hoops/2.0 # centerline distance from the border of the extreme confining hoops
        self.cl_bars = self.e + self.D_bars/2.0 + self.D_hoops # centerline distance from the border of the corner bars
        self.bc = self.b - self.cl_hoops*2
        self.dc = self.d - self.cl_hoops*2
        self.As = ComputeACircle(self.D_hoops)

        # Arguments
        self.rho_s_x = 2.0*ComputeRho(self.As, 1, self.bc*self.s) if rho_s_x == -1 else rho_s_x
        self.rho_s_y = 2.0*ComputeRho(self.As, 1, self.dc*self.s) if rho_s_y == -1 else rho_s_y
        self.Ec = self.ComputeEc() if Ec == -1 else Ec

        # Members
        self.nr_bars = self.ComputeNrBars()
        self.A = self.ComputeA()
        self.Ac = self.ComputeAc()
        self.Ay = ComputeACircle(self.D_bars)
        self.rho_bars = ComputeRho(self.Ay, self.nr_bars, self.A)
        self.Iy = self.ComputeIy()
        self.Iz = self.ComputeIz()

        # Data storage for loading/saving
        self.UpdateStoredData()


    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        self.data = [["INFO_TYPE", "RCRectShape"], # Tag for differentiating different data
            ["name_tag", self.name_tag],
            ["b", self.b],
            ["d", self.d],
            ["bc", self.bc],
            ["dc", self.dc],
            ["L", self.L],
            ["e", self.e],
            ["A", self.A],
            ["Ac", self.Ac],
            ["Iy", self.Iy],
            ["Iz", self.Iz],
            ["fc", self.fc],
            ["Ec", self.Ec],
            ["D_bars", self.D_bars],
            ["nr_bars", self.nr_bars],
            ["Ay", self.Ay],
            ["bars_position_x", self.bars_position_x],
            ["bars_ranges_position_y", self.bars_ranges_position_y],
            ["rho_bars", self.rho_bars],
            ["cl_bars", self.cl_bars],
            ["fy", self.fy],
            ["Ey", self.Ey],
            ["D_hoops", self.D_hoops],
            ["s", self.s],
            ["As", self.As],
            ["rho_s_x", self.rho_s_x],
            ["rho_s_y", self.rho_s_y],
            ["cl_hoops", self.cl_hoops],
            ["fs", self.fs],
            ["Es", self.Es]]


    def ShowInfo(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        print("")
        print("Requested info for RC rectangular section of name tag = {}".format(self.name_tag))
        print("Width of the section b = {} mm".format(self.b/mm_unit))
        print("Depth of the section d = {} mm".format(self.d/mm_unit))
        print("Concrete cover e = {} mm".format(self.e/mm_unit))
        print("Concrete area A = {} mm2".format(self.A/mm2_unit))
        print("Core concrete area Ac = {} mm2".format(self.Ac/mm2_unit))
        print("Unconfined concrete compressive strength fc = {} MPa".format(self.fc/MPa_unit))
        print("Young modulus for concrete Ec = {} GPa".format(self.Ec/GPa_unit))
        print("Diameter of the reinforcing bars D_bars = {} mm and area of one bar Ay = {} mm2 with {} bars".format(self.D_bars/mm_unit, self.Ay/mm2_unit, self.nr_bars))
        print("Diameter of the hoops D_hoops = {} mm and area of one stirrup As = {} mm2".format(self.D_hoops/mm_unit, self.As/mm2_unit))
        print("Ratio of area of longitudinal reinforcement to area of concrete section rho_bars = {}".format(self.rho_bars))
        print("Ratio of area of lateral reinforcement to lateral area of concrete section in x rho_s_x = {} ".format(self.rho_s_x))
        print("Ratio of area of lateral reinforcement to lateral area of concrete section in y rho_s_y = {} ".format(self.rho_s_y))
        print("Moment of inertia of the circular section (strong axis) Iy = {} mm4".format(self.Iy/mm4_unit))
        print("Moment of inertia of the circular section (weak axis) Iz = {} mm4".format(self.Iz/mm4_unit))
        print("")
        

    def ComputeNrBars(self):
        """
        Compute the number of vertical bars in the array bars_position_x (note that this list of lists can have different list sizes).

        @returns int: Number of vertical reinforcing bars.
        """
        nr_bars = 0
        for range in self.bars_position_x:
            nr_bars += np.size(range)-1

        return nr_bars


    def ComputeEc(self):
        """
        Compute Ec using the formula from Mander et Al. 1988.

        @returns float: Young modulus of concrete.
        """

        return 5000.0 * math.sqrt(-self.fc/MPa_unit) * MPa_unit


    def ComputeA(self):
        """
        Compute the area for a rectangular section.

        @returns float: Total area.
        """
        return self.b * self.d


    def ComputeAc(self):
        """
        Compute the confined area (area inside the centerline of the hoops, according to Mander et Al. 1988).

        @returns float: Confined area.
        """
        return self.bc * self.dc


    def ComputeIy(self):
        """
        Compute the moment of inertia of the rectangular section with respect to the strong axis.

        @returns float: Moment of inertia (strong axis)
        """
        return self.b * self.d**3 / 12.0


    def ComputeIz(self):
        """
        Compute the moment of inertia of the rectangular section with respect to the weak axis.

        @returns float: Moment of inertia (weak axis)
        """
        return self.d * self.b**3 / 12.0


class RCSquareShape(RCRectShape):
    """
    Class that is the children of RCRectShape and cover the specific case of square RC sections.  

    @param RCRectShape: Parent class.
    """
    def __init__(self, b, L, e, fc, D_bars, bars_position_x: np.ndarray, bars_ranges_position_y: np.ndarray, fy, Ey, D_hoops, s, fs, Es, name_tag="Not Defined", rho_s_x=-1, rho_s_y=-1, Ec=-1):
        """
        Constructor of the class. It passes the arguments into the parent class to generate the specific case of a aquare RC section.

        @param b (float): Width/depth of the section.
        @param L (float): Effective length of the element associated with this section.
            If the panel zone is present, exclude its dimension.
        @param e (float): Concrete cover.
        @param fc (float): Unconfined concrete compressive strength (cylinder test).
        @param D_bars (float): Diameter of the reinforcing bars.
        @param bars_position_x (np.ndarray):  Distances from border to bar centerline, bar to bar centerlines and
            finally bar centerline to border in the x direction (aligned).
            Starting from the left to right, from the top range to the bottom one.
            The number of bars for each range can vary; in this case, add this argument when defining the array " dtype = object".
        @param bars_ranges_position_y (np.ndarray): Distances from border to range centerlines, range to range centerlines and
            finally range centerline to border in the y direction.
            Starting from the top range to the bottom one.
        @param fy (float): Yield stress for reinforcing bars.
        @param Ey (float): Young modulus for reinforcing bars.
        @param D_hoops (float): Diameter of the hoops.
        @param s (float): Vertical centerline spacing between hoops.
        @param fs (float): Yield stress for the hoops.
        @param Es (float): Young modulus for the hoops
        @param name_tag (str, optional): A nametag for the section. Defaults to "Not Defined".
        @param rho_s_x (float, optional): Ratio of the transversal area of the hoops to the associated concrete area in the x direction.
            Defaults to -1, e.g. computed in __init__()  and ReInit() assuming one range of hoops.
        @param rho_s_y (float, optional): Ratio of the transversal area of the hoops to the associated concrete area in the y direction.
            Defaults to -1, e.g. computed in __init__()  and ReInit() assuming one range of hoops.
        @param Ec (float, optional): Young modulus for concrete. Defaults to -1, e.g. computed in __init__() and ReInit().
        """
        super().__init__(b, b, L, e, fc, D_bars, bars_position_x, bars_ranges_position_y, fy, Ey, D_hoops, s, fs, Es, name_tag, rho_s_x, rho_s_y, Ec)


class RCCircShape(Section):
    """
    Class that stores funcions, geometric and mechanical properties of RC circular shape profile.
    Note that for the validity of the formulas, the hoops needs to be closed (with 135 degress possibly).

    @param Section: Parent abstract class.
    """
    def __init__(self, b, L, e, fc, D_bars, n_bars: int, fy, Ey, D_hoops, s, fs, Es, name_tag = "Not Defined", rho_s_vol = -1, Ec = -1):
        """
        The conctructor of the class.

        @param b (float): Width of the section.
        @param L (float): Effective length of the element associated with this section.
            If the panel zone is present, exclude its dimension.
        @param e (float): Concrete cover.
        @param fc (float): Unconfined concrete compressive strength (cylinder test).
        @param D_bars (float): Diameter of the vertical reinforcing bars.
        @param n_bars (int): Number of vertical reinforcing bars.
        @param fy (float): Yield stress for reinforcing bars.
        @param Ey (float): Young modulus for reinforcing bars.
        @param D_hoops (float): Diameter of the hoops.
        @param s (float): Vertical centerline spacing between hoops.
        @param fs (float): Yield stress for the hoops.
        @param Es (float): Young modulus for the hoops
        @param name_tag (str, optional): A nametag for the section. Defaults to "Not Defined".
        @param rho_s_vol (float, optional): Ratio of the volume of transverse confining steel to the volume of confined concrete core.
            Defaults to -1, e.g. computed according to Mander et Al. 1988.
        @param Ec (float, optional): Young modulus for concrete. Defaults to -1, e.g. computed in __init__() and ReInit().

        @exception NegativeValue: b needs to be positive.
        @exception NegativeValue: L needs to be positive.
        @exception NegativeValue: e needs to be positive.
        @exception PositiveValue: fc needs to be negative.
        @exception NegativeValue: D_bars needs to be positive.
        @exception NegativeValue: n_bars needs to be a positive integer.
        @exception NegativeValue: fy needs to be positive.
        @exception NegativeValue: Ey needs to be positive.
        @exception NegativeValue: D_hoops needs to be positive.
        @exception NegativeValue: s needs to be positive.
        @exception NegativeValue: fs needs to be positive.
        @exception NegativeValue: Es needs to be positive.
        @exception NegativeValue: Ec needs to be positive if different from -1.
        @exception InconsistentGeometry: e should be smaller than half the depth and the width of the section.
        """
        # Check
        if b < 0: raise NegativeValue()
        if L < 0: raise NegativeValue()
        if e < 0: raise NegativeValue()
        if fc > 0: raise PositiveValue()
        if D_bars < 0: raise NegativeValue()
        if n_bars < 0: raise NegativeValue()
        if fy < 0: raise NegativeValue()
        if Ey < 0: raise NegativeValue()
        if D_hoops < 0: raise NegativeValue()
        if s < 0: raise NegativeValue()
        if fs < 0: raise NegativeValue()
        if Es < 0: raise NegativeValue()
        if Ec != -1 and Ec < 0: raise NegativeValue()
        if e > b/2: raise InconsistentGeometry()

        # Arguments
        self.b = b
        self.L = L
        self.e = e
        self.fc = fc
        self.D_bars = D_bars
        self.n_bars = n_bars
        self.fy = fy
        self.Ey = Ey
        self.D_hoops = D_hoops
        self.s = s
        self.fs = fs
        self.Es = Es
        self.name_tag = name_tag

        # Initialized the parameters that are dependent from others
        self.ReInit(rho_s_vol, Ec)


    def ReInit(self, rho_s_vol = -1, Ec = -1):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.

        @param rho_s_vol (float, optional): Ratio of the volume of transverse confining steel to the volume of confined concrete core.
            Defaults to -1, e.g. computed according to Mander et Al. 1988.
        @param Ec (float): Young modulus for concrete. Defaults to -1, e.g. computed according to Mander et Al. 1988.
        """
        # Precompute some members
        self.cl_hoops = self.e + self.D_hoops/2.0 # centerline distance from the border of the extreme confining hoops
        self.cl_bars = self.e + self.D_bars/2.0 + self.D_hoops # centerline distance from the border of the corner bars
        self.bc = self.b - self.cl_hoops*2 # diameter of spiral (hoops) between bar centerline
        self.As = ComputeACircle(self.D_hoops)

        # Arguments
        self.rho_s_vol = self.ComputeRhoVol() if rho_s_vol == -1 else rho_s_vol
        self.Ec = self.ComputeEc() if Ec == -1 else Ec

        # Members
        self.A = ComputeACircle(self.b)
        self.Ac = ComputeACircle(self.bc)
        self.Ay = ComputeACircle(self.D_bars)
        self.rho_bars = ComputeRho(self.Ay, self.n_bars, self.A)
        self.I = self.ComputeI()

        # Data storage for loading/saving
        self.UpdateStoredData()


    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        self.data = [["INFO_TYPE", "RCCircShape"], # Tag for differentiating different data
            ["name_tag", self.name_tag],
            ["b", self.b],
            ["bc", self.bc],
            ["L", self.L],
            ["e", self.e],
            ["A", self.A],
            ["Ac", self.Ac],
            ["I", self.I],
            ["fc", self.fc],
            ["Ec", self.Ec],
            ["D_bars", self.D_bars],
            ["n_bars", self.n_bars],
            ["Ay", self.Ay],
            ["rho_bars", self.rho_bars],
            ["cl_bars", self.cl_bars],
            ["fy", self.fy],
            ["Ey", self.Ey],
            ["D_hoops", self.D_hoops],
            ["s", self.s],
            ["As", self.As],
            ["rho_s_vol", self.rho_s_vol],
            ["cl_hoops", self.cl_hoops],
            ["fs", self.fs],
            ["Es", self.Es]]

    def ShowInfo(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        print("")
        print("Requested info for RC circular section of name tag = {}".format(self.name_tag))
        print("Width of the section b = {} mm".format(self.b/mm_unit))
        print("Concrete cover e = {} mm".format(self.e/mm_unit))
        print("Concrete area A = {} mm2".format(self.A/mm2_unit))
        print("Core concrete area Ac = {} mm2".format(self.Ac/mm2_unit))
        print("Unconfined concrete compressive strength fc = {} MPa".format(self.fc/MPa_unit))
        print("Young modulus for concrete Ec = {} GPa".format(self.Ec/GPa_unit))
        print("Diameter of the reinforcing bars D_bars = {} mm and area of one bar Ay = {} mm2 with {} bars".format(self.D_bars/mm_unit, self.Ay/mm2_unit, self.n_bars))
        print("Diameter of the hoops D_hoops = {} mm and area of one stirrup As = {} mm2".format(self.D_hoops/mm_unit, self.As/mm2_unit))
        print("Ratio of area of longitudinal reinforcement to area of concrete section rho_bars = {} ".format(self.rho_bars))
        print("Ratio of the volume of transverse confining steel to the volume of confined concrete core rho_s = {} ".format(self.rho_s_vol)) 
        print("Moment of inertia of the circular section I = {} mm4".format(self.I/mm4_unit))
        print("")
    

    def ComputeRhoVol(self):
        """
        Compute the ratio of the volume of transverse confining steel to the volume of confined concrete core.
        (according to Mander et Al. 1988).

        @returns float: Ratio.
        """
        vol_s = self.As*math.pi*self.bc
        vol_c = math.pi/4*self.bc**2*self.s

        return vol_s/vol_c


    def ComputeEc(self):
        """
        Compute Ec using the formula from Mander et Al. 1988.

        @returns float: Young modulus of concrete.
        """

        return 5000.0 * math.sqrt(-self.fc/MPa_unit) * MPa_unit


    def ComputeI(self):
        """
        Compute the moment of inertia of the circular section.

        @returns float: Moment of inertia.
        """
        return self.b**4*math.pi/64


def ComputeACircle(D):
    """
    Function that computes the area of one circle (reinforcing bar or hoop).

    @param D (float): Diameter of the circle (reinforcing bar of hoop).

    @returns float: Area the circle (for reinforcing bars or hoops).
    """
    return D**2/4.0*math.pi


def ComputeRho(A, nr, A_tot):
    """
    Compute the ratio of area of a reinforcement to area of a section.

    @param A (float): Area of reinforcement.
    @param nr (float): Number of reinforcement (allow float for computing ratio with different area;
        just convert the other areas to one and compute the equivalent n).
    @param A_tot (float): Area of the concrete.

    @returns float: Ratio.
    """
    return nr * A / A_tot