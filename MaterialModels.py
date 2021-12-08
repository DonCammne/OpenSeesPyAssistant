# Module with the material models
#   Carmine Schipani, 2021

from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import math
from abc import abstractmethod
from copy import deepcopy
from OpenSeesPyHelper.Section import *
from OpenSeesPyHelper.DataManagement import *
from OpenSeesPyHelper.ErrorHandling import *
from OpenSeesPyHelper.Units import *

# Material models

class MaterialModels(DataManagement):
    @abstractmethod
    def CheckApplicability(self):
        """Abstract function used to check the applicability of the material model.
        """
        pass

class ModifiedIMK(MaterialModels):
	# Class that stores funcions and material properties of an double symmetric I-shape profile (Bilin), the default values are valid for a simple cantelever.
    # For more information about the empirical model for the computation of the parameters, see Lignos Krawinkler 2011.
    #TODO: the units are mm and kN and Mkg but they should be m and N and kg (define and multiply with the correct unit, divide with the resulting unit to change it)
    #TODO: RBS not implemented

    global n
    n = 10.0
    
    def __init__(self, ID: int, Type, d, bf, tf, tw, h_1, Iy_mod, iz, E, Fy, Npl, My, L,
        N_G = 0, K_factor = 3, L_0 = -1, L_b = -1, Mc = -1, K = -1, theta_u = -1, safety_factors = False):
        # """
        # Parameters
        # ----------
        # ID : int
        #     The ID of the material model parameters. Can be used as ID for the material model
        # ele : SectionSteelIShape
        #     The element section
        # Lb : double
        #     The distance from the column face to the nearest lateral brace
        # N_G : double
        #     The gravity axial load (default: 0)
        # Cantilever : bool
        #     A parameter to know if the element is a cantilever or not (double fixed) (default: False)
        # Mc : double
        #     The capping moment (default: 0, thus computed here)
        # K : double
        #     The residual strength ratio (default: 0, thus computed here)
        # theta_u : double
        #     The ultimate rotation (default: 0, thus computed here)
        # """
        # Check
        if ID < 0: raise NegativeValue()
        if Type != "Beam" and Type != "Col": raise WrongArgument()
        if d < 0: raise NegativeValue()
        if bf < 0: raise NegativeValue()
        if tf < 0: raise NegativeValue()
        if tw < 0: raise NegativeValue()
        if h_1 < 0: raise NegativeValue()
        if Iy_mod < 0: raise NegativeValue()
        if iz < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if Fy < 0: raise NegativeValue()
        if Npl < 0: raise NegativeValue()
        if My < 0: raise NegativeValue()
        if L < 0: raise NegativeValue()
        if N_G < 0: raise NegativeValue()
        if L_0 != -1 and L_0 < 0: raise NegativeValue()
        if L_b != -1 and L_b < 0: raise NegativeValue()
        if Mc != -1 and Mc < 0: raise NegativeValue()
        if K != -1 and K < 0: raise NegativeValue()
        if theta_u != -1 and theta_u < 0: raise NegativeValue()
        if h_1 > d: raise InconsistentGeometry()
        if N_G > Npl: raise MemberFailure()
        if L_0 > L: raise InconsistentGeometry()

        # Arguments
        self.Type = Type
        self.ID = ID
        self.d = d
        self.bf = bf
        self.tf = tf
        self.tw = tw
        self.h_1 = h_1
        self.Iy_mod = Iy_mod
        self.iz = iz
        self.E = E
        self.Fy = Fy
        self.Npl = Npl
        self.My = My
        self.L = L
        self.N_G = N_G
        self.K_factor = K_factor
        self.L_0 = L if L_0 == -1 else L_0
        self.L_b = L if L_b == -1 else L_b

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        if safety_factors:
            self.gamma_rm = 1.25
            self.prob_factor = 1.15
        else:
            self.gamma_rm = 1.0
            self.prob_factor = 1.0
        self.ReInit(Mc, K, theta_u)


    # Methods
    def ReInit(self, Mc = -1, K = -1, theta_u = -1):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Precompute some members
        self.My_star = self.ComputeMyStar()

        # Arguments
        self.Mc = self.ComputeMc() if Mc == -1 else Mc
        self.K = self.ComputeK() if K == -1 else K
        self.theta_u = self.ComputeTheta_u() if theta_u == -1 else theta_u

        # Check applicability
        self.CheckApplicability()

        # Members
        self.Ke = self.ComputeKe()
        self.theta_y = self.ComputeTheta_y()
        self.theta_p = self.ComputeTheta_p()
        self.theta_pc = self.ComputeTheta_pc()
        self.McMy = self.Mc/self.My_star
        self.rate_det = self.ComputeRefEnergyDissipationCap()
        self.a = self.Computea()
        self.a_s = self.Computea_s()
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"

        # Data storage for loading/saving
        self.UpdateStoredData()
        

    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "ModifiedIMK"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag], 
            ["Type", self.Type],
            ["d", self.d],
            ["bf", self.bf],
            ["tf", self.tf],
            ["tw", self.tw],
            ["h_1", self.h_1],
            ["Iy_mod", self.Iy_mod],
            ["iz", self.iz],
            ["E", self.E],
            ["Fy", self.Fy],
            ["L", self.L],
            ["N_G", self.N_G],
            ["K_factor", self.K_factor],
            ["Ke", self.Ke],
            ["L_0", self.L_0],
            ["L_b", self.L_b],
            ["gamma_rm", self.gamma_rm],
            ["prob_factor", self.prob_factor],
            ["Npl", self.Npl],
            ["My", self.My],
            ["My_star", self.My_star],
            ["Mc", self.Mc],
            ["McMy", self.McMy],
            ["K", self.K],
            ["theta_y", self.theta_y],
            ["theta_p", self.theta_p],
            ["theta_pc", self.theta_pc],
            ["theta_u", self.theta_u],
            ["rate_det", self.rate_det],
            ["a", self.a],
            ["a_s", self.a_s]]


    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        Mr = self.K*self.My_star
        theta_p_plot = self.theta_p
        if self.theta_p > self.theta_u-self.theta_y:
            theta_p_plot = self.theta_u-self.theta_y
        theta_r = self.theta_y + theta_p_plot + self.theta_pc*(1.0-Mr/self.Mc)
        if theta_r > self.theta_u:
            theta_r = self.theta_u
            Mr = self.Mc*(1.0-1.0/self.theta_pc*(self.theta_u-self.theta_y-theta_p_plot))

        print("")
        print("Requested info for Modified IMK (Ibarra-Medina-Krawinkler) material model Parameters, ID = {}".format(self.ID))
        print("Section associated: {}".format(self.section_name_tag))
        print('theta y = {}'.format(self.theta_y))
        print('theta p = {}'.format(self.theta_p))
        print('theta r = {}'.format(theta_r))
        print('theta pc = {}'.format(self.theta_pc))
        print('theta u = {}'.format(self.theta_u))
        print('My star = {} kNm'.format(self.My_star/kNm_unit))
        print('Mc = {} kNm'.format(self.Mc/kNm_unit))
        print('Mr = {} kNm'.format(Mr/kNm_unit))
        print('a = {} '.format(self.a))
        print('as = {} '.format(self.a_s))
        print('lambda (deterioration rate) = {} '.format(self.rate_det))
        print("")
        
        if plot:
            # Data for plotting
            x_axis = np.array([0.0, self.theta_y, self.theta_y + theta_p_plot, theta_r, self.theta_u, self.theta_u])
            x_axis2 = np.array([self.theta_y + theta_p_plot, self.theta_y + theta_p_plot + self.theta_pc])
            y_axis = ([0.0, self.My_star/kNm_unit, self.Mc/kNm_unit, Mr/kNm_unit, Mr/kNm_unit, 0.0])
            y_axis2 = ([self.Mc/kNm_unit, 0.0])

            fig, ax = plt.subplots()
            ax.plot(x_axis, y_axis, 'k-')
            ax.plot(x_axis2, y_axis2, 'k--')

            ax.set(xlabel='Rotation [rad]', ylabel='Moment [kNm]', 
                title='Modified IMK deterioration model (ID={})'.format(self.ID))
            ax.grid()

            if block:
                plt.show()


    def CheckApplicability(self):
        Check = True
        if self.Type == "Beam":
            if self.d/self.tw < 20 or self.d/self.tw > 55:
                Check = False
                print("The d/tw check was not fullfilled")
            if self.L_b/self.iz < 20 or self.L_b/self.iz > 80:
                Check = False
                print("The Lb/iz check was not fullfilled")
            if self.bf/2/self.tf < 4 or self.bf/2/self.tf  > 8:
                Check = False
                print("The bf/2/tf check was not fullfilled")
            if self.L/self.d < 2.5 or self.L/self.d > 7:
                Check = False
                print("The check L/d was not fullfilled")
            if self.d < 0.102 or self.d > 0.914:
                Check = False
                print("The d check was not fullfilled")
            if self.Fy < 240e6 or self.Fy > 450e6:
                Check = False
                print("The Fy check was not fullfilled")
        else:
            if self.h_1/self.tw < 3.71 or self.d/self.tw > 57.5:
                Check = False
                print("The h1/tw check was not fullfilled")
            if self.L_b/self.iz < 38.4 or self.L_b/self.iz > 120:
                Check = False
                print("The Lb/iz check was not fullfilled")
            if self.N_G/self.Npl < 0 or self.N_G/self.Npl > 0.75:
                Check = False
                print("The NG/Npl check was not fullfilled")
        if not Check:
            print("The validity of the equations is not fullfilled.")
            print("!!!!!!! WARNING !!!!!!! Check material model of Modified IMK, ID=", self.ID)
            print("")

    def ComputeKe(self):
        return self.K_factor*n*self.E*self.Iy_mod/self.L

    def Computea(self):
         # strain hardening ratio of spring
        # a = (n+1)*My*(McMy-1)/(K*theta_p)
        return (n+1.0)*self.My_star*(self.McMy-1.0)/(self.Ke*self.theta_p)

    def Computea_s(self):
        # modified strain hardening ratio of spring (Ibarra & Krawinkler 2005)
        # as = a/(1+n*(1-a))
        return self.a/(1.0+n*(1.0-self.a))

    def ComputeMyStar(self):
        if self.Type == "Beam":
            return self.prob_factor*self.My*self.gamma_rm*1.1
        else:
            if self.N_G/self.Npl > 0.2:
                return 1.15*self.prob_factor*self.My*self.gamma_rm*(1-self.N_G/self.Npl)*9.0/8.0
            else:
                return 1.15*self.prob_factor*self.My*self.gamma_rm*(1-self.N_G/2.0/self.Npl)

    def ComputeMc(self):
        if self.Type == "Beam":
            return self.My_star*1.11
            # For RBS: My_star*1.09
        else:
            tmp = 12.5*(self.h_1/self.tw)**(-0.2)*(self.L_b/self.iz)**(-0.4)*(1-self.N_G/self.Npl)**0.4
            return max(min(1.3, tmp), 1.0)*self.My_star

    def ComputeK(self):
        if self.Type == "Beam":
            return 0.4
        else:
            tmp = 0.5-0.4*self.N_G/self.Npl
            return max(tmp, 0)
    
    def ComputeTheta_y(self):
        return self.My_star/self.Ke*(n+1)

    def ComputeTheta_p(self):
        if self.Type == "Beam":
            if self.d < 533.0*mm_unit:
                return 0.0865*(self.h_1/self.tw)**(-0.365)*(self.bf/2.0/self.tf)**(-0.14)*(self.L_0/self.d)**(0.34)*(self.d/533.0*mm_unit)**(-0.721)*(self.Fy/355.0*MPa_unit)**(-0.23)
            else:
                return 0.318*(self.h_1/self.tw)**(-0.550)*(self.bf/2.0/self.tf)**(-0.345)*(self.L_0/self.d)**(0.090)*(self.L_b/self.iz)**(-0.023)*(self.d/533.0*mm_unit)**(-0.330)*(self.Fy/355.0*MPa_unit)**(-0.130)
                # With RBS: ...
        else:
            tmp = 294.0*(self.h_1/self.tw)**(-1.7)*(self.L_b/self.iz)**(-0.7)*(1.0-self.N_G/self.Npl)**(1.6) # *(self.E/self.Fy/gamma_rm)**(0.2) # EC8
            if tmp > 0.2:
                tmp = 0.2
            # if tmp > self.theta_u-self.theta_y:
            #     tmp = (self.theta_u-self.theta_y)*0.799 # convergence issue
            return tmp

    def ComputeTheta_pc(self):
        if self.Type == "Beam":
            if self.d < 533.0*mm_unit:
                return 5.63*(self.h_1/self.tw)**(-0.565)*(self.bf/2.0/self.tf)**(-0.800)*(self.d/533.0*mm_unit)**(-0.280)*(self.Fy/355.0*MPa_unit)**(-0.430)
            else:
                return 7.50*(self.h_1/self.tw)**(-0.610)*(self.bf/2.0/self.tf)**(-0.710)*(self.L_b/self.iz)**(-0.110)*(self.d/533.0*mm_unit)**(-0.161)*(self.Fy/355.0*MPa_unit)**(-0.320)
                # With RBS: ...
        else:
            tmp =  90.0*(self.h_1/self.tw)**(-0.8)*(self.L_b/self.iz)**(-0.8)*(1.0-self.N_G/self.Npl)**(2.5) # *(self.E/self.Fy/gamma_rm)**(0.07) # EC8
            if tmp > 0.3:
                tmp = 0.3
            return tmp

    def ComputeTheta_u(self):
        if self.Type == "Beam":
            return 0.2
        else:
            return 0.15

    def ComputeRefEnergyDissipationCap(self):
        if self.Type == "Beam":
            if self.d < 533.0*mm_unit:
                return 495.0*(self.h_1/self.tw)**(-1.34)*(self.bf/2.0/self.tf)**(-0.595)*(self.Fy/355.0*MPa_unit)**(-0.360)
            else:
                return 536.0*(self.h_1/self.tw)**(-1.26)*(self.bf/2.0/self.tf)**(-0.525)*(self.L_b/self.iz)**(-0.130)*(self.Fy/355.0*MPa_unit)**(-0.291)
                # With RBS: ...
        else:
            if self.N_G/self.Npl > 0.35:
                tmp = 268000.0*(self.h_1/self.tw)**(-2.30)*(self.L_b/self.iz)**(-1.130)*(1.0-self.N_G/self.Npl)**(1.19)
                if tmp > 3.0:
                    tmp = 3.0
                return tmp
            else:
                tmp = 25000.0*(self.h_1/self.tw)**(-2.14)*(self.L_b/self.iz)**(-0.53)*(1.0-self.N_G/self.Npl)**(4.92)
                if tmp > 3.0:
                    tmp = 3.0
                return tmp


    def Bilin(self):
        # Generate the material model Bilin (Modified IMK) using the computed parameters
        # uniaxialMaterial("Bilin", IDMat, K, asPos, asNeg, MyPos, MyNeg, LS, LK, LA, LD, cS, cK, cA, cD, th_pP, th_pN, th_pcP, th_pcN, ResP, ResN, th_uP, th_uN, DP, DN)
        #   ID         Material Identification (integer)
        #   K          Initial stiffness after the modification for n (see Ibarra and Krawinkler, 2005)
        #   asPos      Strain hardening ratio after n modification (see Ibarra and Krawinkler, 2005)
        #   asNeg      Strain hardening ratio after n modification (see Ibarra and Krawinkler, 2005)
        #   MyPos      Positive yield moment (with sign)
        #   MyNeg      Negative yield moment (with sign)
        #   LS = 1000  Basic strength deterioration parameter (see Lignos and Krawinkler, 2009) (a very large # = no cyclic deterioration)
        #   LK = 1000  Unloading stiffness deterioration parameter (see Lignos and Krawinkler, 2009) (a very large # = no cyclic deterioration)
        #   LA = 1000  Accelerated reloading stiffness deterioration parameter (see Lignos and Krawinkler, 2009) (a very large # = no cyclic deterioration)
        #   LD = 1000  Post-capping strength deterioration parameter (see Lignos and Krawinkler, 2009) (a very large # = no cyclic deterioration)
        #   cS = 1     Exponent for basic strength deterioration (c = 1.0 for no deterioration)
        #   cK = 1     Exponent for unloading stiffness deterioration (c = 1.0 for no deterioration)
        #   cA = 1     Exponent for accelerated reloading stiffness deterioration (c = 1.0 for no deterioration)
        #   cD = 1     Exponent for post-capping strength deterioration (c = 1.0 for no deterioration)
        #   th_pP      Plastic rotation capacity for positive loading direction (exemple 0.025)
        #   th_pN      Plastic rotation capacity for negative loading direction (exemple 0.025)
        #   th_pcP     Post-capping rotation capacity for positive loading direction (exemple 0.3)
        #   th_pcN     Post-capping rotation capacity for negative loading direction (exemple 0.3)
        #   KP         Residual strength ratio for positive loading direction (exemple 0.4)
        #   KN         Residual strength ratio for negative loading direction (exemple 0.4)
        #   th_uP      Ultimate rotation capacity for positive loading direction (exemple 0.4)
        #   th_uN      Ultimate rotation capacity for negative loading direction (exemple 0.4)
        #   rateDetP   Rate of cyclic deterioration for positive loading direction (exemple 1.0)
        #   rateDetN   Rate of cyclic deterioration for negative loading direction (exemple 1.0)
        uniaxialMaterial("Bilin", self.ID, self.Ke, self.a_s, self.a_s, self.My_star, -1.0*self.My_star,
            1., 1., 1., 1., 1., 1., 1., 1., self.theta_p, self.theta_p, self.theta_pc, self.theta_pc,
            self.K, self.K, self.theta_u, self.theta_u, self.rate_det, self.rate_det)


class ModifiedIMKSteelIShape(ModifiedIMK):
    def __init__(self, ID, section: SteelIShape, N_G = 0, K_factor = 3, L_0 = -1, L_b = -1, Mc = -1, K = -1, theta_u = -1, safety_factors = False):
        super().__init__(ID, section.Type, section.d, section.bf, section.tf, section.tw, section.h_1,
            section.Iy_mod, section.iz, section.E, section.Fy, section.Npl, section.My, section.L, N_G,
            K_factor, L_0, L_b, Mc, K, theta_u, safety_factors)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()
    

class Gupta1999(MaterialModels):
    ######################################################################################
    ## PZRotSpringMaterialModel
    ######################################################################################
    ## Class that stores funcions and material properties of a panel zone rotational spring (Hysteresis). For more information, see Gupta 1999
    ## Warning: he units should be mm and kN
    ##          Carmine Schipani, 2021
    ##
    ##  ID :            Unique material model ID
    ##  col :           Object from the class SectionSteelIShape of the column
    ##  d_beam :        Beam depth
    ##  a_s :           Assumed strain hardening (default: 0.03)
    ##  Ry :            Expected value for yield strength (default: 1.2)
    ##  pinchx :        Pinching factor for strain (or deformation) during reloading (default: 1.0) 
    ##  pinchy :        Pinching factor for stress (or force) during reloading (default: 1.0) 
    ##  dmg1 :          Damage due to ductility: D1(mu-1) (default: 0.0) 
    ##  dmg2 :          Damage due to energy: D2(Eii/Eult) (default: 0.0) 
    ##  beta :          Power used to determine the degraded unloading stiffness based on ductility, mu-beta (default: 0.0) 
    ##  plot :          Bool for having the material Backbone graph (default: False)
    ##  block :         Bool for having the plots all at once (default: False)

    def __init__(self, ID: int, d_c, bf_c, tf_c, I_c, d_b, tf_b, Fy, E, t_p,
        t_dp = 0.0, a_s = 0.03, pinchx = 1.0, pinchy = 1.0, dmg1 = 0.0, dmg2 = 0.0, beta = 0.0, safety_factor = False):
        #TODO: check applicability and specify applicable for continous column

        # # Initialisation
        # E = col.E          # Young modulus
        # Fy = col.Fy_web    # Yield strength of the column web (assume continous column)
        # dc = col.d         # Column depth
        # bf_c = col.bf      # Column flange width
        # tf_c = col.tf      # Column flange thickness
        # tp = col.tw        # Panel zone thickness
        # Ic = col.Iy        # Column moment of inertia (strong axis)
        # tpz = tp+t_dp      # Thickness of the panel zone and the doubler plate

        # Check
        if ID < 1: raise WrongArgument()
        if d_c < 0: raise NegativeValue()
        if bf_c < 0: raise NegativeValue()
        if tf_c < 0: raise NegativeValue()
        if d_b < 0: raise NegativeValue()
        if tf_b < 0: raise NegativeValue()
        if Fy < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if t_p < 0: raise NegativeValue()
        if a_s < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.d_c = d_c
        self.bf_c = bf_c
        self.tf_c = tf_c
        self.I_c = I_c
        self.d_b = d_b
        self.tf_b = tf_b
        self.Fy = Fy
        self.E = E
        self.t_p = t_p
        self.t_dp = t_dp
        self.a_s = a_s
        self.pinchx = pinchx
        self.pinchy = pinchy
        self.dmg1 = dmg1
        self.dmg2 = dmg2
        self.beta = beta
        if safety_factor:
            self.Ry = 1.2
        else:
            self.Ry = 1.0

        # Initialized the parameters that are dependent from others
        self.beam_section_name_tag = "None"
        self.col_section_name_tag = "None"
        self.ReInit()


    # Methods
    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Check applicability
        self.CheckApplicability()

        # Members
        if self.beam_section_name_tag != "None": self.beam_section_name_tag = self.beam_section_name_tag + " (modified)"
        if self.col_section_name_tag != "None": self.col_section_name_tag = self.col_section_name_tag + " (modified)"

        # Trilinear Parameters
        self.t_pz = self.t_p + self.t_dp
        self.Vy = 0.55 * self.Fy * self.Ry * self.d_c * self.t_pz # Yield Shear
        self.G = self.E/(2.0 * (1.0 + 0.30)) # Shear Modulus
        self.Ke = 0.95 * self.G * self.t_pz * self.d_c # Elastic Stiffness
        self.Kp = 0.95 * self.G * self.bf_c * (self.tf_c * self.tf_c) / self.d_b # Plastic Stiffness

        # Define Trilinear Equivalent Rotational Spring
        # Yield point for Trilinear Spring at gamma1_y
        self.gamma1_y = self.Vy / self.Ke
        self.M1y = self.gamma1_y * (self.Ke * self.d_b)
        # Second Point for Trilinear Spring at 4 * gamma1_y
        self.gamma2_y = 4.0 * self.gamma1_y
        self.M2y = self.M1y + (self.Kp * self.d_b) * (self.gamma2_y - self.gamma1_y)
        # Third Point for Trilinear Spring at 100 * gamma1_y
        self.gamma3_y = 100.0 * self.gamma1_y
        self.M3y = self.M2y + (self.a_s * self.Ke * self.d_b) * (self.gamma3_y - self.gamma2_y)

        # Data storage for loading/saving
        self.UpdateStoredData()


    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "Gupta1999"], # Tag for differentiating different data
            ["ID", self.ID],
            ["beam_section_name_tag", self.beam_section_name_tag], 
            ["col_section_name_tag", self.col_section_name_tag], 
            ["d_c", self.d_c],
            ["bf_c", self.bf_c],
            ["tf_c", self.tf_c],
            ["I_c", self.I_c],
            ["d_b", self.d_b],
            ["tf_b", self.tf_b],
            ["Fy", self.Fy],
            ["E", self.E],
            ["G", self.G],
            ["t_p", self.t_p],
            ["t_dp", self.t_dp],
            ["t_pz", self.t_pz],
            ["a_s", self.a_s],
            ["pinchx", self.pinchx],
            ["pinchy", self.pinchy],
            ["dmg1", self.dmg1],
            ["dmg2", self.dmg2],
            ["beta", self.beta],
            ["Ry", self.Ry],
            ["Vy", self.Vy],
            ["Ke", self.Ke],
            ["Kp", self.Kp],
            ["gamma1_y", self.gamma1_y],
            ["M1y", self.M1y],
            ["gamma2_y", self.gamma2_y],
            ["M2y", self.M2y],
            ["gamma3_y", self.gamma3_y],
            ["M3y", self.M3y]]

    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        print("")
        print("Requested info for Gupta 1999 material model Parameters, ID = {}".format(self.ID))
        print("Sections associated, column: {} ".format(self.col_section_name_tag))
        print("Sections associated, beam: {} ".format(self.beam_section_name_tag))
        print("gamma1_y = {}".format(self.gamma1_y))
        print("gamma2_y = {}".format(self.gamma2_y))
        print("gamma3_y = {}".format(self.gamma3_y))
        print("M1y = {} kNm".format(self.M1y/kNm_unit))
        print("M2y = {} kNm".format(self.M2y/kNm_unit))
        print("M3y = {} kNm".format(self.M3y/kNm_unit))
        print("")
        
        if plot:
            # Data for plotting
            # Last point for plot
            gamma3_y_plot = 10.0 * self.gamma1_y
            M3y_plot = self.M2y + (self.a_s * self.Ke * self.d_b) * (gamma3_y_plot - self.gamma2_y)

            x_axis = np.array([0.0, self.gamma1_y, self.gamma2_y, gamma3_y_plot])
            y_axis = ([0.0, self.M1y/kNm_unit, self.M2y/kNm_unit, M3y_plot/kNm_unit])

            fig, ax = plt.subplots()
            ax.plot(x_axis, y_axis, 'k-')

            ax.set(xlabel='Rotation [rad]', ylabel='Moment [kNm]', 
                title='Gupta 1999 material model (ID={})'.format(self.ID))
            ax.grid()

            if block:
                plt.show()


    def CheckApplicability(self):
        #TODO: 
        Check = True
        # if self.d/self.tw < 20 or self.d/self.tw > 55:
        #     Check = False
        #     print("The d/tw check was not fullfilled")
        if not Check:
            print("The validity of the equations is not fullfilled.")
            print("!!!!!!! WARNING !!!!!!! Check material model of Gupta 1999, ID=", self.ID)
            print("")


    def Hysteretic(self):
        # Hysteretic Material
        # $matTag       integer tag identifying material
        # $s1p $e1p     stress and strain (or force & deformation) at first point of the envelope in the positive direction
        # $s2p $e2p     stress and strain (or force & deformation) at second point of the envelope in the positive direction
        # $s3p $e3p     stress and strain (or force & deformation) at third point of the envelope in the positive direction (optional)
        # $s1n $e1n     stress and strain (or force & deformation) at first point of the envelope in the negative direction
        # $s2n $e2n     stress and strain (or force & deformation) at second point of the envelope in the negative direction
        # $s3n $e3n     stress and strain (or force & deformation) at third point of the envelope in the negative direction (optional)
        # $pinchx       pinching factor for strain (or deformation) during reloading
        # $pinchy       pinching factor for stress (or force) during reloading
        # $damage1      damage due to ductility: D1(mu-1)
        # $damage2      damage due to energy: D2(Eii/Eult)
        # $beta         power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0) 

        uniaxialMaterial("Hysteretic", self.ID,
            self.M1y, self.gamma1_y, self.M2y, self.gamma2_y, self.M3y, self.gamma3_y,
            -self.M1y, -self.gamma1_y, -self.M2y, -self.gamma2_y, -self.M3y, -self.gamma3_y,
            self.pinchx, self.pinchy, self.dmg1, self.dmg2, self.beta)


class Gupta1999SteelIShape(Gupta1999):
    def __init__(self, ID: int, col: SteelIShape, beam: SteelIShape,
        t_dp = 0.0, a_s = 0.03, pinchx = 1.0, pinchy = 1.0, dmg1 = 0.0, dmg2 = 0.0, beta = 0.0, safety_factor = False):
        super().__init__(ID, col.d, col.bf, col.tf, col.Iy, beam.d, beam.tf, col.Fy_web, col.E, col.tw,
            t_dp, a_s, pinchx, pinchy, dmg1, dmg2, beta, safety_factor)
        self.beam_section_name_tag = beam.name_tag
        self.col_section_name_tag = col.name_tag
        self.UpdateStoredData()


class Skiadopoulos2021(MaterialModels):
    ######################################################################################
    ## PZRotSpringMaterialModel
    ######################################################################################
    ## Class that stores funcions and material properties of a panel zone rotational spring (Hysteresis). For more information, see Gupta 1999
    ## Warning: he units should be mm and kN
    ##          Carmine Schipani, 2021
    ##
    ##  ID :            Unique material model ID
    ##  col :           Object from the class SectionSteelIShape of the column
    ##  d_beam :        Beam depth
    ##  a_s :           Assumed strain hardening (default: 0.03)
    ##  Ry :            Expected value for yield strength (default: 1.2)
    ##  pinchx :        Pinching factor for strain (or deformation) during reloading (default: 1.0) 
    ##  pinchy :        Pinching factor for stress (or force) during reloading (default: 1.0) 
    ##  dmg1 :          Damage due to ductility: D1(mu-1) (default: 0.0) 
    ##  dmg2 :          Damage due to energy: D2(Eii/Eult) (default: 0.0) 
    ##  beta :          Power used to determine the degraded unloading stiffness based on ductility, mu-beta (default: 0.0) 
    ##  plot :          Bool for having the material Backbone graph (default: False)
    ##  block :         Bool for having the plots all at once (default: False)

    global Kf_Ke_tests, Cw1_tests, Cf1_tests, Cw4_tests, Cf4_tests, Cw6_tests, Cf6_tests

    Kf_Ke_tests = [1.000, 0.153, 0.120, 0.090, 0.059, 0.031, 0.019, 0.009, 0.005, 0.004, 0.000]
    Kf_Ke_tests.reverse()
    Cw1_tests = [0.96, 0.96, 0.955, 0.94, 0.93, 0.90, 0.89, 0.89, 0.88, 0.88, 0.88]
    Cw1_tests.reverse()
    Cf1_tests = [0.035, 0.035, 0.033, 0.031, 0.018, 0.015, 0.013, 0.009, 0.009, 0.010, 0.010]
    Cf1_tests.reverse()
    Cw4_tests = [1.145, 1.145, 1.140, 1.133, 1.120, 1.115, 1.115, 1.11, 1.10, 1.10, 1.10]
    Cw4_tests.reverse()
    Cf4_tests = [0.145, 0.145, 0.123, 0.111, 0.069, 0.040, 0.040, 0.018, 0.010, 0.012, 0.012]
    Cf4_tests.reverse()
    Cw6_tests = [1.205, 1.2050, 1.2000, 1.1925, 1.1740, 1.1730, 1.1720, 1.1690, 1.1670, 1.1650, 1.1650]
    Cw6_tests.reverse()
    Cf6_tests = [0.165, 0.1650, 0.1400, 0.1275, 0.0800, 0.0500, 0.0500, 0.0180, 0.0140, 0.0120, 0.0120]
    Cf6_tests.reverse()

    def __init__(self, ID: int, d_c, bf_c, tf_c, I_c, d_b, tf_b, Fy, E, t_p,
        t_dp = 0.0, a_s = 0.03, pinchx = 1.0, pinchy = 1.0, dmg1 = 0.0, dmg2 = 0.0, beta = 0.0, safety_factor = False):
        #TODO: check applicability and specify applicable for continous column

        # # Initialisation
        # E = col.E          # Young modulus
        # Fy = col.Fy_web    # Yield strength of the column web (assume continous column)
        # dc = col.d         # Column depth
        # bf_c = col.bf      # Column flange width
        # tf_c = col.tf      # Column flange thickness
        # tp = col.tw        # Panel zone thickness
        # Ic = col.Iy        # Column moment of inertia (strong axis)
        # tpz = tp+t_dp      # Thickness of the panel zone and the doubler plate

        # Check
        if ID < 1: raise WrongArgument()
        if d_c < 0: raise NegativeValue()
        if bf_c < 0: raise NegativeValue()
        if tf_c < 0: raise NegativeValue()
        if d_b < 0: raise NegativeValue()
        if tf_b < 0: raise NegativeValue()
        if Fy < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if t_p < 0: raise NegativeValue()
        if a_s < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.d_c = d_c
        self.bf_c = bf_c
        self.tf_c = tf_c
        self.I_c = I_c
        self.d_b = d_b
        self.tf_b = tf_b
        self.Fy = Fy
        self.E = E
        self.t_p = t_p
        self.t_dp = t_dp
        self.a_s = a_s
        self.pinchx = pinchx
        self.pinchy = pinchy
        self.dmg1 = dmg1
        self.dmg2 = dmg2
        self.beta = beta
        if safety_factor:
            self.Ry = 1.2
        else:
            self.Ry = 1.0

        # Initialized the parameters that are dependent from others
        self.beam_section_name_tag = "None"
        self.col_section_name_tag = "None"
        self.ReInit()


    # Methods
    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Check applicability
        self.CheckApplicability()

        # Memebers
        if self.beam_section_name_tag != "None": self.beam_section_name_tag = self.beam_section_name_tag + " (modified)"
        if self.col_section_name_tag != "None": self.col_section_name_tag = self.col_section_name_tag + " (modified)"
        self.t_pz = self.t_p + self.t_dp
        self.G = self.E/(2.0 * (1.0 + 0.30)) # Shear Modulus

        # Refined computation of the parameters for the backbone curve for the panel zone spring (Skiadopoulos et al. (2021))
        # Panel Zone Elastic Stiffness
        self.Ks = self.t_pz*(self.d_c-self.tf_c)*self.G
        self.Kb = 12.0*self.E*(self.I_c+self.t_dp*(self.d_c-2.0*self.tf_c)**3/12.0)/(self.d_b-0)**2
        self.Ke = self.Ks*self.Kb/(self.Ks+self.Kb)

        # Column Flange Stiffness
        self.Ksf = 2.0*(self.tf_c*self.bf_c*self.G)
        self.Kbf = 2.0*(12.0*self.E*self.bf_c*self.tf_c**3/12.0/(self.d_b-0)**2)
        self.Kf = self.Ksf*self.Kbf/(self.Ksf+self.Kbf)

        # Kf/Ke Calculation for Panel Zone Categorization
        self.Kf_Ke = self.Kf/self.Ke

        # Panel Zone Strength Coefficients (results from tests for a_w_eff and a_f_eff)
        self.Cw1 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cw1_tests)
        self.Cf1 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cf1_tests)
        self.Cw4 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cw4_tests)
        self.Cf4 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cf4_tests)
        self.Cw6 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cw6_tests)
        self.Cf6 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cf6_tests)

        # Panel Zone Model
        self.V1 = self.Fy*self.Ry/math.sqrt(3)*(self.Cw1*(self.d_c-self.tf_c)*self.t_pz + self.Cf1*2*(self.bf_c-self.t_p)*self.tf_c)
        self.V4 = self.Fy*self.Ry/math.sqrt(3)*(self.Cw4*(self.d_c-self.tf_c)*self.t_pz + self.Cf4*2*(self.bf_c-self.t_p)*self.tf_c)
        self.V6 = self.Fy*self.Ry/math.sqrt(3)*(self.Cw6*(self.d_c-self.tf_c)*self.t_pz + self.Cf6*2*(self.bf_c-self.t_p)*self.tf_c)

        self.M1 = self.V1*(self.d_b-self.tf_b)
        self.M4 = self.V4*(self.d_b-self.tf_b)
        self.M6 = self.V6*(self.d_b-self.tf_b)

        self.Gamma_1 = self.V1/self.Ke
        self.Gamma_4 = 4*self.Gamma_1
        self.Gamma_6 = 6*self.Gamma_1

        # Data storage for loading/saving
        self.UpdateStoredData()


    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "Skiadopoulos2021"], # Tag for differentiating different data
            ["ID", self.ID],
            ["beam_section_name_tag", self.beam_section_name_tag], 
            ["col_section_name_tag", self.col_section_name_tag], 
            ["d_c", self.d_c],
            ["bf_c", self.bf_c],
            ["tf_c", self.tf_c],
            ["I_c", self.I_c],
            ["d_b", self.d_b],
            ["tf_b", self.tf_b],
            ["Fy", self.Fy],
            ["E", self.E],
            ["G", self.G],
            ["t_p", self.t_p],
            ["t_dp", self.t_dp],
            ["t_pz", self.t_pz],
            ["a_s", self.a_s],
            ["pinchx", self.pinchx],
            ["pinchy", self.pinchy],
            ["dmg1", self.dmg1],
            ["dmg2", self.dmg2],
            ["beta", self.beta],
            ["Ry", self.Ry],
            ["Ks", self.Ks],
            ["Kb", self.Kb],
            ["Ke", self.Ke],
            ["Ksf", self.Ksf],
            ["Kbf", self.Kbf],
            ["Kf", self.Kf],
            ["Kf_Ke", self.Kf_Ke],
            ["V1", self.V1],
            ["V4", self.V4],
            ["V6", self.V6],
            ["M1", self.M1],
            ["M4", self.M4],
            ["M6", self.M6],
            ["Gamma_1", self.Gamma_1],
            ["Gamma_4", self.Gamma_4],
            ["Gamma_6", self.Gamma_6]]


    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        print("")
        print("Requested info for Skiadopoulos 2021 material model Parameters, ID = {}".format(self.ID))
        print("Sections associated, column: {} ".format(self.col_section_name_tag))
        print("Sections associated, beam: {} ".format(self.beam_section_name_tag))
        print("Gamma_1 = {}".format(self.Gamma_1))
        print("Gamma_4 = {}".format(self.Gamma_4))
        print("Gamma_6 = {}".format(self.Gamma_6))
        print("M1 = {} kNm".format(self.M1/kNm_unit))
        print("M4 = {} kNm".format(self.M4/kNm_unit))
        print("M6 = {} kNm".format(self.M6/kNm_unit))
        print("")
        
        if plot:
            # Data for plotting
            x_axis = np.array([0.0, self.Gamma_1, self.Gamma_4, self.Gamma_6])
            y_axis = ([0.0, self.M1/kNm_unit, self.M4/kNm_unit, self.M6/kNm_unit])

            fig, ax = plt.subplots()
            ax.plot(x_axis, y_axis, 'k-')

            ax.set(xlabel='Rotation [rad]', ylabel='Moment [kNm]', 
                title='Skiadopoulos 2021 material model (ID={})'.format(self.ID))
            ax.grid()

            if block:
                plt.show()


    def CheckApplicability(self):
        #TODO: 
        Check = True
        # if self.d/self.tw < 20 or self.d/self.tw > 55:
        #     Check = False
        #     print("The d/tw check was not fullfilled")
        if not Check:
            print("The validity of the equations is not fullfilled.")
            print("!!!!!!! WARNING !!!!!!! Check material model of Skiadopoulos 2021, ID=", self.ID)
            print("")


    def Hysteretic(self):
        # Refined backbone curve for the panel zone spring (Skiadopoulos et al. (2021))

        # Hysteretic Material
        # $matTag       integer tag identifying material
        # $s1p $e1p     stress and strain (or force & deformation) at first point of the envelope in the positive direction
        # $s2p $e2p     stress and strain (or force & deformation) at second point of the envelope in the positive direction
        # $s3p $e3p     stress and strain (or force & deformation) at third point of the envelope in the positive direction (optional)
        # $s1n $e1n     stress and strain (or force & deformation) at first point of the envelope in the negative direction
        # $s2n $e2n     stress and strain (or force & deformation) at second point of the envelope in the negative direction
        # $s3n $e3n     stress and strain (or force & deformation) at third point of the envelope in the negative direction (optional)
        # $pinchx       pinching factor for strain (or deformation) during reloading
        # $pinchy       pinching factor for stress (or force) during reloading
        # $damage1      damage due to ductility: D1(mu-1)
        # $damage2      damage due to energy: D2(Eii/Eult)
        # $beta         power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0) 

        uniaxialMaterial("Hysteretic", self.ID,
            self.M1, self.Gamma_1, self.M4, self.Gamma_4, self.M6, self.Gamma_6,
            -self.M1, -self.Gamma_1, -self.M4, -self.Gamma_4, -self.M6, -self.Gamma_6,
            self.pinchx, self.pinchy, self.dmg1, self.dmg2, self.beta)


class Skiadopoulos2021SteelIShape(Skiadopoulos2021):
    def __init__(self, ID: int, col: SteelIShape, beam: SteelIShape,
        t_dp=0, a_s=0.03, pinchx=1, pinchy=1, dmg1=0, dmg2=0, beta=0, safety_factor=False):
        super().__init__(ID, col.d, col.bf, col.tf, col.Iy, beam.d, beam.tf, col.Fy_web, col.E, col.tw,
            t_dp=t_dp, a_s=a_s, pinchx=pinchx, pinchy=pinchy, dmg1=dmg1, dmg2=dmg2, beta=beta, safety_factor=safety_factor)
        self.beam_section_name_tag = beam.name_tag
        self.col_section_name_tag = col.name_tag
        self.UpdateStoredData()
        

class UnconfMander1988(MaterialModels):
    # Class that stores funcions and material properties of a rectangular shape RC  profile (Concrete04). For more information see Mander et Al. 1988
    # Warning: the units should be m and N
      #TODO: validity check: warning if concrete fc is bigger than something (see article Lee)
    
    def __init__(self, ID: int, fc, Ec, ec = 1, ecp = 1, fct = -1, et = -1, beta = 0.1, safety_factors = False):
        #TODO: security factor
        # self.ec = -0.002
        # self.ecp = 2.0*self.ec                                  # concrete spalling
        # self.esu = 0.05                                         # strain capacity of stirrups 0.012-0.05 experimental results (BOOK BEYER)

        # Check
        if ID < 0: raise NegativeValue()
        if fc > 0: raise PositiveValue()
        if Ec < 0: raise NegativeValue()
        if ec != 1 and ec > 0: raise PositiveValue()
        if ecp != 1 and ecp > 0: raise PositiveValue()
        if fct != -1 and fct < 0: raise NegativeValue()
        if et != -1 and et < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.fc = fc
        self.Ec = Ec
        self.ec = -0.002 if ec == 1 else ec
        self.beta = beta

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        if safety_factors:
            #TODO: insert the correct value and see where to put it
            self.Ry = 1.5
        else:
            self.Ry = 1.0
        self.ReInit(ecp, fct, et)

    # Methods
    def ReInit(self, ecp = 1, fct = -1, et = -1):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Check applicability
        self.CheckApplicability()

        # Arguments
        self.ecp = self.Compute_ecp() if ecp == 1 else ecp
        self.fct = self.Compute_fct() if fct == -1 else fct
        self.et = self.Compute_et() if et == -1 else et

        # Members
        self.ecu = self.Compute_ecu()
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"

        # Data storage for loading/saving
        self.UpdateStoredData()


    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "UnconfMander1988"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["fc", self.fc],
            ["Ec", self.Ec],
            ["ec", self.ec],
            ["ecp", self.ecp],
            ["ecu", self.ecu],
            ["fct", self.fct],
            ["et", self.et],
            ["Ry", self.Ry],
            ["beta", self.beta]]


    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        print("")
        print("Requested info for Unconfined Mander 1988 material model Parameters, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print('Concrete strength fc = {} MPa'.format(self.fc/MPa_unit))
        print('Strain at maximal strength ec = {}'.format(self.ec))
        print('Maximal strain ecu = {}'.format(self.ecu))
        print("")

        if plot:
            fig, ax = plt.subplots()
            PlotConcrete04(self.fc, self.Ec, self.ec, self.ecu, "U", ax, self.ID)

            if block:
                plt.show()


    def CheckApplicability(self):
        #TODO: 
        Check = True
        # if self.d/self.tw < 20 or self.d/self.tw > 55:
        #     Check = False
        #     print("The d/tw check was not fullfilled")
        if not Check:
            print("The validity of the equations is not fullfilled.")
            print("!!!!!!! WARNING !!!!!!! Check material model of Unconfined Mander 1988, ID=", self.ID)
            print("")


    def Compute_ecp(self):
        return 2.0*self.ec


    def Compute_fct(self):
        #TODO: check units
        return 0.30 * math.pow(-self.fc/MPa_unit, 2/3) * MPa_unit               # SIA 262 2013


    def Compute_et(self):
        return self.fct/self.Ec                               # Mander Eq 45


    def Compute_ecu(self):
        #TODO: add refernce Katrin book
        return -0.004     # FROM BRIDGE BOOK OF KATRIN BEYER!


    def Concrete04(self):
        # Generate the material model using the given parameters

        # Define Concrete04 Popovics Concrete material model for unconfined concrete
        # uniaxialMaterial("Concrete04", matTag, fc, ec, ecu, Ec, <fct et> <beta>)
        # matTag     integer tag identifying material
        # fc     floating point values defining concrete compressive strength at 28 days (compression is negative)*
        # ec     floating point values defining concrete strain at maximum strength*
        # ecu    floating point values defining concrete strain at crushing strength*
        # Ec     floating point values defining initial stiffness**
        # fct    floating point value defining the maximum tensile strength of concrete
        # et     floating point value defining ultimate tensile strain of concrete
        # beta   loating point value defining the exponential curve parameter to define the residual stress (as a factor of ft) at etu 
        
        uniaxialMaterial("Concrete04", self.ID, self.fc, self.ec, self.ecu, self.Ec, self.fct, self.et, self.beta)
        #TODO: ft, et for unconf?


def Concrete04Funct(fc, discretized_eps, ec, Ec):
    x = discretized_eps/ec
    r = Ec / (Ec - fc/ec)
    return fc*x*r / (r-1+x**r)


def PlotConcrete04(fc, Ec, ec, ecu, Type, ax, ID = 0):
    # For ax, just use:
    # fig, ax = plt.subplots()
    # to generate the figure
    # Type is "C" or "U" for confined or unconfined

    if Type == "C":
        name = "Confined"
    elif Type == "U":
        name = "Unconfined"
    else:
        raise NameError("Type should be C or U (ID={})".format(ID))

    # Data for plotting
    N = int(-ecu*1e5)
    x_axis = np.zeros(N)
    y_axis = np.zeros(N)
    for i in range(N):
        x_axis[i] = -i/1e5
        y_axis[i] = Concrete04Funct(fc, x_axis[i], ec, Ec)


    ax.plot(-x_axis*100.0, -y_axis/MPa_unit, 'k-', label = name)
    ax.set(xlabel='Strain [%]', ylabel='Stress MPa]', 
                    title='Mander 1988 material model (ID={})'.format(ID))
    plt.legend()
    plt.grid()


class UnconfMander1988RCRectShape(UnconfMander1988):
    def __init__(self, ID: int, ele: RCRectShape, ec=1, ecp=1, fct=-1, et=-1, beta=0.1, safety_factors=False):
        super().__init__(ID, ele.fc, ele.Ec, ec=ec, ecp=ecp, fct=fct, et=et, beta=beta, safety_factors=safety_factors)
        self.section_name_tag = ele.name_tag
        self.UpdateStoredData()


class ConfMander1988(MaterialModels):
    # Class that stores funcions and material properties of a rectangular shape RC  profile (Concrete04). For more information see Mander et Al. 1988
    # Warning: the units should be m and N
      #TODO: validity check: warning if concrete fc is bigger than something (see article Lee)
    
    def __init__(self, ID: int, bc, dc, Ac, fc, Ec, nr_bars, D_bars, wx_top: np.ndarray, wx_bottom: np.ndarray, wy: np.ndarray, s, D_hoops, rho_s_x, rho_s_y, fs, 
        ec = 1, ecp = 1, fct = -1, et = -1, esu = -1, beta = 0.1, k1 = 4.1, safety_factors = False):
        #TODO: security factor
        # self.ec = -0.002
        # self.ecp = 2.0*self.ec                                  # concrete spalling
        # self.esu = 0.05                                         # strain capacity of stirrups 0.012-0.05 experimental results (BOOK BEYER)
        # wx_top : numpy.ndarray
            # A vector that defines the distance between bars in x direction (NOT CENTERLINE DISTANCE). One range of bars implemented
        # wx_bottom : numpy.ndarray
            # A vector that defines the distance between bars in x direction (NOT CENTERLINE DISTANCE). One range of bars implemented
        # wy : numpy.ndarray
            # A vector that defines the distance between bars in y direction (NOT CENTERLINE DISTANCE). One range of bars implemented

        # Check
        if ID < 0: raise NegativeValue()
        if bc < 0: raise NegativeValue()
        if dc < 0: raise NegativeValue()
        if Ac < 0: raise NegativeValue()
        if fc > 0: raise PositiveValue()
        if Ec < 0: raise NegativeValue()
        if nr_bars < 0: raise NegativeValue()
        if D_bars < 0: raise NegativeValue()
        if s < 0: raise NegativeValue()
        if D_hoops < 0: raise NegativeValue()
        if rho_s_x < 0: raise NegativeValue()
        if rho_s_y < 0: raise NegativeValue()
        if fs < 0: raise NegativeValue()
        if ec != 1 and ec > 0: raise PositiveValue()
        if ecp != 1 and ecp > 0: raise PositiveValue()
        if fct != -1 and fct < 0: raise NegativeValue()
        if et != -1 and et < 0: raise NegativeValue()
        if esu != -1 and esu < 0: raise NegativeValue()
        if k1 < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.bc = bc
        self.dc = dc
        self.Ac = Ac
        self.fc = fc
        self.Ec = Ec
        self.nr_bars = nr_bars
        self.D_bars = D_bars
        self.wx_top = deepcopy(wx_top)
        self.wx_bottom = deepcopy(wx_bottom)
        self.wy = deepcopy(wy)
        self.s = s
        self.D_hoops = D_hoops
        self.rho_s_x = rho_s_x
        self.rho_s_y = rho_s_y
        self.fs = fs
        self.ec = -0.002 if ec == 1 else ec
        self.esu = 0.05 if esu == -1 else esu
        self.beta = beta
        self.k1 = k1 # 4.1 from Richart et al. 1928 or 5.6 from Balmer 1949 

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        if safety_factors:
            #TODO: insert the correct value
            self.Ry = 1.5
        else:
            self.Ry = 1.0
        self.ReInit(ecp, fct, et)

    def ReInit(self, ecp = 1, fct = -1, et = -1):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Check applicability
        self.CheckApplicability()

        # Arguments
        self.ecp = self.Compute_ecp() if ecp == 1 else ecp
        self.fct = self.Compute_fct() if fct == -1 else fct
        self.et = self.Compute_et() if et == -1 else et

        # Members
        #TODO: check formulas and units one time more
        self.ecu = self.Compute_ecu()
        self.k2 = 5.0*self.k1
        self.Ai = self.ComputeAi()
        self.Ae = (self.Ac - self.Ai) * (1.0 - (self.s-self.D_hoops)/2.0/self.bc)*(1.0 - (self.s-self.D_hoops)/2.0/self.dc)
        self.rho_cc = self.nr_bars*self.D_bars**2/4.0*math.pi / self.Ac
        self.Acc = self.Ac*(1.0-self.rho_cc)
        self.ke = self.Ae/self.Acc
        self.fl_x = -self.rho_s_x * self.fs
        self.fl_y = -self.rho_s_y * self.fs
        #TODO: K_combo = combination of flx and fly factor (interpolation graph)
        self.fl_prime = (self.fl_x + self.fl_y)/2.0 * self.ke # derive the interpolating curve
        self.K_combo = -1.254 + 2.254 * math.sqrt(1.0+7.94*self.fl_prime/self.fc) - 2.0*self.fl_prime/self.fc # in Mander, it has a prime
        self.fcc = self.fc * self.K_combo
        self.ecc = (1.0 + 5.0 * (self.K_combo-1.0)) * self.ec
        self.eccu = -0.004 + (1.4*(self.rho_s_x+self.rho_s_y)*self.esu*self.fs) / self.fcc # FROM BRIDGE BOOK OF KATRIN BEYER!
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "ConfMander1988"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["bc", self.bc],
            ["dc", self.dc],
            ["Ac", self.Ac],
            ["fc", self.fc],
            ["Ec", self.Ec],
            ["ec", self.ec],
            ["ecp", self.ecp],
            ["ecu", self.ecu],
            ["fct", self.fct],
            ["et", self.et],
            ["k1", self.k1],
            ["fcc", self.fcc],
            ["ecc", self.ecc],
            ["eccu", self.eccu],
            ["Ry", self.Ry],
            ["beta", self.beta],
            ["nr_bars", self.nr_bars],
            ["D_bars", self.D_bars],
            ["wx_top", self.wx_top],
            ["wx_bottom", self.wx_bottom],
            ["wy", self.wy],
            ["s", self.s],
            ["D_hoops", self.D_hoops],
            ["rho_s_x", self.rho_s_x],
            ["rho_s_y", self.rho_s_y],
            ["fs", self.fs],
            ["esu", self.esu]]


    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        print("")
        print("Requested info for Confined Mander 1988 material model Parameters, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print('Concrete strength fc = {} MPa'.format(self.fc/MPa_unit))
        print('Concrete strength confined fcc = {} MPa'.format(self.fcc/MPa_unit))
        print('Strain at maximal strength ec = {}'.format(self.ec))
        print('Strain at maximal strength confined ecc = {}'.format(self.ecc))
        print('Maximal strain ecu = {}'.format(self.ecu))
        print('Maximal strain confined eccu = {}'.format(self.eccu))
        print("")

        if plot:
            fig, ax = plt.subplots()
            PlotConcrete04(self.fcc, self.Ec, self.ecc, self.eccu, "C", ax, self.ID)

            if block:
                plt.show()


    def CheckApplicability(self):
        #TODO: 
        Check = True
        if len(self.wy) == 0 or len(self.wx_top) == 0 or len(self.wx_bottom) == 0: 
            Check = False
            print("Hypothesis of one bar per corner not fullfilled.")
        if not Check:
            print("The validity of the equations is not fullfilled.")
            print("!!!!!!! WARNING !!!!!!! Check material model of Confined Mander 1988, ID=", self.ID)
            print("")


    def Compute_ecp(self):
        return 2.0*self.ec


    def Compute_fct(self):
        #TODO: check units. Also, use fcc?
        return 0.30 * math.pow(-self.fc/MPa_unit, 2/3) * MPa_unit               # SIA 262 2013


    def Compute_et(self):
        return self.fct/self.Ec                               # Mander Eq 45


    def Compute_ecu(self):
        #TODO: add refernce Katrin book
        return -0.004     # FROM BRIDGE BOOK OF KATRIN BEYER!


    def ComputeAi(self):
        return ( np.sum(np.multiply(self.wy, self.wy))*2.0 +
            np.sum(np.multiply(self.wx_top, self.wx_top)) +
            np.sum(np.multiply(self.wx_bottom, self.wx_bottom)) ) / 6.0 # ineffectual area


    def Concrete04(self):
        # Generate the material model using the given parameters

        # Define Concrete04 Popovics Concrete material model for confined concrete
        # uniaxialMaterial("Concrete04", matTag, fcc, ecc, eccu, Ec, <fct et> <beta>)
        # matTag     integer tag identifying material
        # fcc    floating point values defining concrete compressive strength at 28 days (compression is negative)*
        # ecc    floating point values defining concrete strain at maximum strength*
        # eccu   floating point values defining concrete strain at crushing strength*
        # Ec     floating point values defining initial stiffness**
        # fct    floating point value defining the maximum tensile strength of concrete
        # et     floating point value defining ultimate tensile strain of concrete
        # beta   loating point value defining the exponential curve parameter to define the residual stress (as a factor of ft) at etu 

        uniaxialMaterial("Concrete04", self.ID, self.fcc, self.ecc, self.eccu, self.Ec, self.fct, self.et, self.beta)
        # TODO: ft, et for conf (change way it is compute because it untilize fcc not fc??


class ConfMander1988RCRectShape(ConfMander1988):
    def __init__(self, ID: int, ele: RCRectShape, ec=1, ecp=1, fct=-1, et=-1, esu=-1, beta=0.1, k1=4.1, safety_factors=False):
        ranges = ele.bars_ranges_position_y
        bars = ele.bars_position_x
        wy = self.__Compute_w(ranges, ele.D_bars)
        wx_top = self.__Compute_w(bars[0], ele.D_bars)
        wx_bottom = self.__Compute_w(bars[-1], ele.D_bars)

        super().__init__(ID, ele.bc, ele.dc, ele.Ac, ele.fc, ele.Ec, ele.nr_bars, ele.D_bars, wx_top, wx_bottom, wy, ele.s, ele.D_hoops, ele.rho_s_x, ele.rho_s_y, ele.fs,
            ec=ec, ecp=ecp, fct=fct, et=et, esu=esu, beta=beta, k1=k1, safety_factors=safety_factors)
        self.section_name_tag = ele.name_tag
        self.UpdateStoredData()

    def __Compute_w(self, vector, D_bars):
        l = len(vector)
        w = np.zeros(l-2)
        for i, elem in enumerate(vector[1:l-1]):
            w[i] = elem - D_bars
        return w


class UniaxialBilinear(MaterialModels):
    # Class that stores funcions and material properties of a simple uniaxial bilinear model (Steel01). 
    # Warning: the units should be m and N
    #     #     Strain hardening factor (default: 0.01)


    def __init__(self, ID: int, fy, Ey, b = 0.01, safety_factors = False):
        #TODO: security factor

        # Check
        if ID < 0: raise NegativeValue()
        if fy < 0: raise NegativeValue()
        if Ey < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.fy = fy
        self.Ey = Ey
        self.b = b

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        if safety_factors:
            #TODO: insert the correct value
            self.Ry = 1.5
        else:
            self.Ry = 1.0
        self.ReInit()

    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Check applicability
        self.CheckApplicability()

        # Arguments

        # Members
        self.ey = self.fy / self.Ey
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "UniaxialBilinear"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["fy", self.fy],
            ["Ey", self.Ey],
            ["ey", self.ey],
            ["b", self.b],
            ["Ry", self.Ry]]


    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        print("")
        print("Requested info for Uniaxial Bilinear material model Parameters, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print('Yielding stress fy = {} MPa'.format(self.fy/MPa_unit))
        print('Young modulus Ey = {} MPa'.format(self.Ey/MPa_unit))
        print('Maximal elastic strain epsilon y = {}'.format(self.ey))
        print('Hardening factor b = {}'.format(self.b))
        print("")

        if plot:
            # Data for plotting
            e_pl = 10.0 * self.ey # to show that if continues with this slope
            sigma_pl = self.b * self.Ey * e_pl

            x_axis = ([0.0, self.ey*100, (self.ey+e_pl)*100])
            y_axis = ([0.0, self.fy/MPa_unit, (self.fy+sigma_pl)/MPa_unit])

            fig, ax = plt.subplots()
            ax.plot(x_axis, y_axis, 'k-')

            ax.set(xlabel='Strain [%]', ylabel='Stress [MPa]', 
                title='Uniaxial Bilinear model for material ID={}'.format(self.ID))
            ax.grid()

            if block:
                plt.show()


    def CheckApplicability(self):
        #TODO: 
        Check = True
        # if len(self.wy) == 0 or len(self.wx_top) == 0 or len(self.wx_bottom) == 0: 
        #     Check = False
        #     print("Hypothesis of one bar per corner not fullfilled.")
        if not Check:
            print("The validity of the equations is not fullfilled.")
            print("!!!!!!! WARNING !!!!!!! Check material model of Uniaxial Bilinear, ID=", self.ID)
            print("")


    def Steel01(self):
        # Generate the material model using the given parameters

        # Define a uniaxial bilinear steel material object with kinematic hardening and optional isotropic hardening described by a non-linear evolution equation (REF: Fedeas).
        # uniaxialMaterial Steel01 $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
        # $matTag 	integer tag identifying material
        # $Fy 	yield strength
        # $E0 	initial elastic tangent
        # $b 	strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
        # $a1 	isotropic hardening parameter, increase of compression yield envelope as proportion of yield strength after a plastic strain of $a2*($Fy/E0). (optional)
        # $a2 	isotropic hardening parameter (see explanation under $a1). (optional).
        # $a3 	isotropic hardening parameter, increase of tension yield envelope as proportion of yield strength after a plastic strain of $a4*($Fy/E0). (optional)
        # $a4 	isotropic hardening parameter (see explanation under $a3). (optional) 

        uniaxialMaterial("Steel01", self.ID, self.fy, self.Ey, self.b)


class UniaxialBilinearSteelIShape(UniaxialBilinear):
    def __init__(self, ID: int, ele: SteelIShape, b=0.01, safety_factors=False):
        super().__init__(ID, ele.Fy, ele.E, b=b, safety_factors=safety_factors)
        self.section_name_tag = ele.name_tag
        self.UpdateStoredData()


class GMP1970(MaterialModels):
    # Class that stores funcions and material properties of GMP1970. For more information see Giuffré, Menegotto and Pinto 1970 and Carreno et Al. 2020
    # Warning: the units should be m and N
    
    def __init__(self, ID: int, fy, Ey, b = 0.02, R0 = 20, cR1 = 0.9, cR2 = 0.08, a1 = 0.039, a2 = 1.0, a3 = 0.029, a4 = 1.0, safety_factors = False):
        #TODO: security factor
        # OpenSees Docu suggest b = 0.015, R0 = 10, cR1 = 0.925, cR2 = 0.15 and a-parameters as default

        # Check
        if ID < 0: raise NegativeValue()

        # Arguments
        self.ID = ID
        self.fy = fy
        self.Ey = Ey
        self.b = b
        self.R0 = R0
        self.cR1 = cR1
        self.cR2 = cR2
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        if safety_factors:
            #TODO: insert the correct value
            self.Ry = 1.5
        else:
            self.Ry = 1.0
        self.ReInit()

    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Check applicability
        self.CheckApplicability()

        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "GMP1970"], # Tag for differentiating different data
            ["ID", self.ID],
            ["section_name_tag", self.section_name_tag],
            ["fy", self.fy],
            ["Ey", self.Ey],
            ["b", self.b],
            ["R0", self.R0],
            ["cR1", self.cR1],
            ["cR2", self.cR2],
            ["a1", self.a1],
            ["a2", self.a2],
            ["a3", self.a3],
            ["a4", self.a4],
            ["Ry", self.Ry]]


    def ShowInfo(self):
        """Function that show the data stored in the class in the command window and plots the material model (optional).
        """
        print("")
        print("Requested info for GMP1970 (Giuffré-Menegotto-Pinto) material model Parameters, ID = {}".format(self.ID))
        print("Section associated: {} ".format(self.section_name_tag))
        print("Yield stress fy = {} MPa".format(self.fy/MPa_unit))
        print("Young modulus Ey = {} MPa".format(self.Ey/MPa_unit))
        print("Strain hardening ratio b = {}".format(self.b))
        print("Bauschinger effect factors R0 = {}, cR1 = {} and cR2 = {}".format(self.R0, self.cR1, self.cR2))
        print("Isotropic hardening factors a1 = {}, a2 = {}, a3 = {} and a4 = {}".format(self.a1, self.a2, self.a3, self.a4))
        print("")

        #TODO: add plot option (difficult to implement)
        # if plot:
        #     # Data for plotting
        #     e_pl = 10.0 * self.ey # to show that if continues with this slope
        #     sigma_pl = self.b * self.Ey * e_pl

        #     x_axis = ([0.0, self.ey*100, (self.ey+e_pl)*100])
        #     y_axis = ([0.0, self.fy/MPa_unit, (self.fy+sigma_pl)/MPa_unit])

        #     fig, ax = plt.subplots()
        #     ax.plot(x_axis, y_axis, 'k-')

        #     ax.set(xlabel='Strain [%]', ylabel='Stress [MPa]', 
        #         title='GMP1970 model for material ID={}'.format(self.ID))
        #     ax.grid()
            

        #     if block:
        #         plt.show()


    def CheckApplicability(self):
        #TODO: 
        Check = True
        # if len(self.wy) == 0 or len(self.wx_top) == 0 or len(self.wx_bottom) == 0: 
        #     Check = False
        #     print("Hypothesis of one bar per corner not fullfilled.")
        if not Check:
            print("The validity of the equations is not fullfilled.")
            print("!!!!!!! WARNING !!!!!!! Check material model of GMP1970, ID=", self.ID)
            print("")


    def Steel02(self):
        # Generate the material model using the given parameters

        # Define uniaxial Giuffre-Menegotto-Pinto steel material object with isotropic strain hardening.
        # uniaxialMaterial Steel02 $matTag $Fy $E $b $R0 $cR1 $cR2 <$a1 $a2 $a3 $a4 $sigInit>
        # $matTag 	    integer tag identifying material
        # $Fy 	        yield strength
        # $E0 	        initial elastic tangent
        # $b 	        strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
        # $R0 $CR1 $CR2 parameters to control the transition from elastic to plastic branches.
        #                       Recommended values: $R0=between 10 and 20, $cR1=0.925, $cR2=0.15
        # $a1 	        isotropic hardening parameter, increase of compression yield envelope as proportion of yield strength after a plastic strain of $a2*($Fy/E0). (optional)
        # $a2 	        isotropic hardening parameter (see explanation under $a1). (optional default = 1.0).
        # $a3 	        isotropic hardening parameter, increase of tension yield envelope as proportion of yield strength after a plastic strain of $a4*($Fy/E0). (optional default = 0.0)
        # $a4 	        isotropic hardening parameter (see explanation under $a3). (optional default = 1.0)
        # $sigInit 	    Initial Stress Value (optional, default: 0.0) the strain is calculated from epsP=$sigInit/$E
        #                   if (sigInit!= 0.0) { double epsInit = sigInit/E; eps = trialStrain+epsInit; } else eps = trialStrain; 

        uniaxialMaterial('Steel02', self.ID, self.fy, self.Ey, self.b, *[self.R0, self.cR1, self.cR2], self.a1, self.a2, self.a3, self.a4)


class GMP1970RCRectShape(GMP1970):
    def __init__(self, ID: int, ele: RCRectShape, b=0.02, R0=20.0, cR1=0.9, cR2=0.08, a1=0.039, a2=1.0, a3=0.029, a4=1.0, safety_factors=False):
        super().__init__(ID, ele.fy, ele.Ey, b=b, R0=R0, cR1=cR1, cR2=cR2, a1=a1, a2=a2, a3=a3, a4=a4, safety_factors=safety_factors)
        self.section_name_tag = ele.name_tag
        self.UpdateStoredData()




#TODO: CalibratedUVC SteelIShape?

# class NAME(MaterialModels):
#     # Class that stores funcions and material properties of (...). For more information see (...)
#     # Warning: the units should be m and N
    
#     def __init__(self, ID: int, ...):
#         #TODO: security factor

#         # Check
#         if ID < 0: raise NegativeValue()

#         # Arguments
#         self.ID = ID

#         # Initialized the parameters that are dependent from others
#         self.section_name_tag = "None"
#         if safety_factors:
#             #TODO: insert the correct value
#             self.Ry = 1.5
#         else:
#             self.Ry = 1.0
#         self.ReInit()

#     def ReInit(self):
#         """Function that computes the value of the parameters that are computed with respect of the arguments.
#         Use after changing the value of argument inside the class (to update the values accordingly). 
#         This function can be very useful in combination with the function "deepcopy()" from the module "copy".
#         """
#         # Check applicability
#         self.CheckApplicability()

#         # Arguments

#         # Members
#         if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"


#         # Data storage for loading/saving
#         self.UpdateStoredData()


#     # Methods
#     def UpdateStoredData(self):
#         self.data = [["INFO_TYPE", "NAME"], # Tag for differentiating different data
#             ["ID", self.ID],
#             ["section_name_tag", self.section_name_tag]
#             ]


#     def ShowInfo(self, plot = False, block = False):
#         """Function that show the data stored in the class in the command window and plots the material model (optional).
#         """
#         print("")
#         print("Requested info for NAME material model Parameters, ID = {}".format(self.ID))
#         print("Section associated: {} ".format(self.section_name_tag))
#         print("")

#         if plot:
#             # Data for plotting
#             e_pl = 10.0 * self.ey # to show that if continues with this slope
#             sigma_pl = self.b * self.Ey * e_pl

#             x_axis = ([0.0, self.ey*100, (self.ey+e_pl)*100])
#             y_axis = ([0.0, self.fy/MPa_unit, (self.fy+sigma_pl)/MPa_unit])

#             fig, ax = plt.subplots()
#             ax.plot(x_axis, y_axis, 'k-')

#             ax.set(xlabel='Strain [%]', ylabel='Stress [MPa]', 
#                 title='Uniaxial Bilinear model for material ID={}'.format(self.ID))
#             ax.grid()
            

#             if block:
#                 plt.show()


#     def CheckApplicability(self):
#         #TODO: 
#         Check = True
#         # if len(self.wy) == 0 or len(self.wx_top) == 0 or len(self.wx_bottom) == 0: 
#         #     Check = False
#         #     print("Hypothesis of one bar per corner not fullfilled.")
#         if not Check:
#             print("The validity of the equations is not fullfilled.")
#             print("!!!!!!! WARNING !!!!!!! Check material model of NAME, ID=", self.ID)
#             print("")


#     def UNIAXIALMATFUNCT(self):
#         # Generate the material model using the given parameters

#         # Define Concrete04 Popovics Concrete material model for confined concrete
#         # uniaxialMaterial("Concrete04", matTag, fcc, ecc, eccu, Ec, <fct et> <beta>)
#         # matTag     integer tag identifying material
#         # fcc    floating point values defining concrete compressive strength at 28 days (compression is negative)*
#         # ecc    floating point values defining concrete strain at maximum strength*
#         # eccu   floating point values defining concrete strain at crushing strength*
#         # Ec     floating point values defining initial stiffness**
#         # fct    floating point value defining the maximum tensile strength of concrete
#         # et     floating point value defining ultimate tensile strain of concrete
#         # beta   loating point value defining the exponential curve parameter to define the residual stress (as a factor of ft) at etu 

#         uniaxialMaterial("Concrete04", self.ID, self.fcc, self.ecc, self.eccu, self.Ec, self.fct, self.et, self.beta)


# class CHILD(NAME):
#     def __init__(self, ID: int, ele: RCRectShape, ec=1, ecp=1, fct=-1, et=-1, esu=-1, beta=0.1, k1=4.1, safety_factors=False):
#         ranges = ele.bars_ranges_position_y
#         bars = ele.bars_position_x
#         wy = self.__Compute_w(ranges, ele.D_bars)
#         wx_top = self.__Compute_w(bars[0], ele.D_bars)
#         wx_bottom = self.__Compute_w(bars[-1], ele.D_bars)

#         super().__init__(ID, ele.bc, ele.dc, ele.Ac, ele.fc, ele.Ec, ele.nr_bars, ele.D_bars, wx_top, wx_bottom, wy, ele.s, ele.D_hoops, ele.rho_s_x, ele.rho_s_y, ele.fs,
#             ec=ec, ecp=ecp, fct=fct, et=et, esu=esu, beta=beta, k1=k1, safety_factors=safety_factors)
#         self.section_name_tag = ele.name_tag
#         self.UpdateStoredData()