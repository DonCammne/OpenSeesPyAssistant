# Module for the Section classes
#          Carmine Schipani, 2021

from OpenSeesPyHelper.DataManagement import DataManagement
import numpy as np
import math

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
    def __init__(self, Type, d, bf, tf, tw, L, r, E, Fy, Fy_wed = 0, NameTAG = ""):
        """The conctructor of the class

        Args:
            Type (str): The type of element (Col or Beam)
            d (double): The depth
            bf (double): The flange's width
            tf (double): The flange's thickness
            tw (double): The web's thickness
            L (double): The element length (exclude the panel zone is present)
            r (double): The radius of the weld fillets
            E (double): The Young modulus
            Fy (double): The yield strength
            Fy_wed (double, optional): The yield strength of the web. Defaults to 0 (e.g. computed in __init__()).
            NameTAG (str, optional): A nametag for the section. Defaults to "".
        """
        pass

    def ReIter(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "copy()" from the module "copy".
        """
        pass

    def ShowInfo(self):
        """Function that show the data stored in the class in the command window.
        """
        pass



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
        #TODO: add the function
        pass

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

    def __init__(self, b, d, L, e, fc, D_bars, bars_position_x: np.ndarray, bars_ranges_position_y: np.ndarray, fy, Ey, D_hoops, s, fs, Es, NameTAG = "", rho_s_x = 0, rho_s_y = 0, Ec = 0):
        """The conctructor of the class.

        Args:
            b (double): Width of the section
            d (double): Depth of the section
            L (double): Length of the element
            e (double): Concrete cover
            fc (double): Unconfined concrete compressive strength (cylinder test)
            D_bars (double): Diameter of the reinforcing bars
            bars_position_x (np.ndarray): Distances from bar to bar centerlines or border to bar centerline in the x direction (aligned).
                Starting from the left to right, from the top range to the bottom one. The number of bars for each range can vary
            bars_ranges_position_y (np.ndarray): Distances from range to range centerlines or border to range centerline in the y direction.
                Starting from the top range to the bottom one
            fy (double): Yield stress for reinforcing bars
            Ey (double): Young modulus for reinforcing bars
            D_hoops (double): Diameter of the hoops
            s (double): Centerline distance for the hoops
            fs (double): Yield stress for the hoops
            Es (double): Young modulus for the hoops
            NameTAG (str, optional): A nametag for the section. Defaults to "".
            rho_s_x (double, optional): Volumetric or not?????????????. Defaults to 0 (e.g. computed in __init__() assuming one range of hoops).
            rho_s_y (double, optional): Volumetric or not?????????????. Defaults to 0 (e.g. computed in __init__() assuming one range of hoops).
            Ec (double, optional): Young modulus for concrete. Defaults to 0 (e.g. computed in __init__()).
        """
        # Check: number row of pos_x == length position_ranges
        # Check: sum (row i of pos_x)-bf < tol (1-2 mm)
        # Note that for the validity of the formulas, at least one bar per corner and at least one hoop closed with 135 degress!

        # wx : numpy.ndarray
        #     A vector that defines the distance between bars in x direction (NOT CENTERLINE DISTANCE). One range of bars implemented
        # wy : numpy.ndarray
        #     A vector that defines the distance between bars in y direction (NOT CENTERLINE DISTANCE). One range of bars implemented
        
        # wx_top, wx_bottom and wy computed with bars_position_x (list of list) and bars_ranges_position_y (list)
        # nr_bars computed
        # if rho_s not defined, assume one range of hoops
        #TODO: ask someone good with concrete (Katrin Beyer, Diego,... ) if rho_s is volumtric or not (and how to compute it)
        pass

    def ReIter(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "copy()" from the module "copy".
        """
        pass

    def ShowInfo(self, d):
        """Function that show the data stored in the class in the command window.
        """
        print(d)


    def ComputeEc(self):
        """Compute Ec using the formula from Mander et Al. 1988.

        Returns:
            double: Young modulus of concrete
        """

        return 5.0 * math.sqrt(-self.fc*1000)


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


    def ComputeAy(self):
        """Compute the area of one reinforcing bar.

        Returns:
            double: Area of one bar
        """
        return self.D_bars**2/4.0*math.pi


    def ComputeAs(self):
        """Compute the area of one hoop.

        Returns:
            double: Area of one hoop
        """     
        return self.D_hoops**2/4.0*math.pi


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

    def ComputeRho(self, A, nr, A_tot):
        """Compute the ratio of area of a reinforcement to area of a section.

        Args:
            A (double): Area of reinforcement
            nr (int): Number of reinforcement
            A_tot (double): Area of the concrete

        Returns:
            double: Ratio
        """
        return nr * A / A_tot


