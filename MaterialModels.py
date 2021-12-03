# Module with the material models
#   Carmine Schipani, 2021

from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import math
from Section import *
from DataManagement import *
from ErrorHandling import *
from Units import *

# Material models

class MaterialModels(DataManagement):
    pass

class ModifiedIMK(MaterialModels):
	# Class that stores funcions and material properties of an double symmetric I-shape profile, the default values are valid for a simple cantelever.
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
        self.L_0 = L_0
        if self.L_0 == -1: self.L_0 = L
        self.L_b = L_b
        if self.L_b == -1: self.L_b = L
        # self.Mc = Mc
        # self.K = K
        # self.theta_u = theta_u

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        if safety_factors:
            self.gamma_rm = 1.25
            self.prob_factor = 1.15
        else:
            self.gamma_rm = 1.0
            self.prob_factor = 1.0
        # self.ReInit(self.Mc, self.K, self.theta_u)
        self.ReInit(Mc, K, theta_u)


    # Methods
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
            print("The empirical equation validity are not fullfilled.")
            print("Check material model of Modified IMK, ID=", self.ID)
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

    def ReInit(self, Mc = -1, K = -1, theta_u = -1):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "copy()" from the module "copy".
        """
        # Precompute some members
        self.My_star = self.ComputeMyStar()

        # Arguments
        self.Mc = Mc
        self.K = K
        self.theta_u = theta_u
        if self.Mc == -1: self.Mc = self.ComputeMc()
        if self.K == -1: self.K = self.ComputeK()
        if self.theta_u == -1: self.theta_u = self.ComputeTheta_u()

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

        # Data storage for loading/saving
        self.data = ["ModifiedIMK", # Tag for differentiating different datas
            self.ID,
            self.section_name_tag, 
            self.Type,
            self.d,
            self.bf,
            self.tf,
            self.tw,
            self.h_1,
            self.Iy_mod,
            self.iz,
            self.E,
            self.Fy,
            self.L,
            self.N_G,
            self.K_factor,
            self.Ke,
            self.L_0,
            self.L_b,
            self.gamma_rm,
            self.prob_factor,
            self.Npl,
            self.My,
            self.My_star,
            self.Mc,
            self.McMy,
            self.K,
            self.theta_y,
            self.theta_p,
            self.theta_pc,
            self.theta_u,
            self.rate_det,
            self.a,
            self.a_s]

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
        print("IMK Material Model Parameters, ID = {}".format(self.ID))
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
                title='Moodified IMK deterioration model (ID={})'.format(self.ID))
            ax.grid()

            if block:
                plt.show()

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
    

# class PZRotSpringMaterialModel(DataManagement):
#     ######################################################################################
#     ## PZRotSpringMaterialModel
#     ######################################################################################
#     ## Class that stores funcions and material properties of a panel zone rotational spring. For more information, see Gupta 1999
#     ## Warning: he units should be mm and kN
#     ##          Carmine Schipani, 2021
#     ##
#     ##  ID :            Unique material model ID
#     ##  col :           Object from the class SectionSteelIShape of the column
#     ##  d_beam :        Beam depth
#     ##  a_s :           Assumed strain hardening (default: 0.03)
#     ##  Ry :            Expected value for yield strength (default: 1.2)
#     ##  pinchx :        Pinching factor for strain (or deformation) during reloading (default: 1.0) 
#     ##  pinchy :        Pinching factor for stress (or force) during reloading (default: 1.0) 
#     ##  dmg1 :          Damage due to ductility: D1(mu-1) (default: 0.0) 
#     ##  dmg2 :          Damage due to energy: D2(Eii/Eult) (default: 0.0) 
#     ##  beta :          Power used to determine the degraded unloading stiffness based on ductility, mu-beta (default: 0.0) 
#     ##  plot :          Bool for having the material Backbone graph (default: False)
#     ##  block :         Bool for having the plots all at once (default: False)

#     global Kf_Ke_tests, Cw1_tests, Cf1_tests, Cw4_tests, Cf4_tests, Cw6_tests, Cf6_tests

#     Kf_Ke_tests = [1.000, 0.153, 0.120, 0.090, 0.059, 0.031, 0.019, 0.009, 0.005, 0.004, 0.000]
#     Kf_Ke_tests.reverse()
#     Cw1_tests = [0.96, 0.96, 0.955, 0.94, 0.93, 0.90, 0.89, 0.89, 0.88, 0.88, 0.88]
#     Cw1_tests.reverse()
#     Cf1_tests = [0.035, 0.035, 0.033, 0.031, 0.018, 0.015, 0.013, 0.009, 0.009, 0.010, 0.010]
#     Cf1_tests.reverse()
#     Cw4_tests = [1.145, 1.145, 1.140, 1.133, 1.120, 1.115, 1.115, 1.11, 1.10, 1.10, 1.10]
#     Cw4_tests.reverse()
#     Cf4_tests = [0.145, 0.145, 0.123, 0.111, 0.069, 0.040, 0.040, 0.018, 0.010, 0.012, 0.012]
#     Cf4_tests.reverse()
#     Cw6_tests = [1.205, 1.2050, 1.2000, 1.1925, 1.1740, 1.1730, 1.1720, 1.1690, 1.1670, 1.1650, 1.1650]
#     Cw6_tests.reverse()
#     Cf6_tests = [0.165, 0.1650, 0.1400, 0.1275, 0.0800, 0.0500, 0.0500, 0.0180, 0.0140, 0.0120, 0.0120]
#     Cf6_tests.reverse()

#     def __init__(self, ID, col: SectionSteelIShape, d_beam, tf_b, t_dp = 0.0, a_s = 0.03, Ry = 1.2, pinchx = 1.0, pinchy = 1.0, dmg1 = 0.0, dmg2 = 0.0, beta = 0.0):
#         self.ID = ID
#         self.col = col
#         self.d_beam = d_beam
#         self.tf_b = tf_b
#         self.t_dp = t_dp
#         self.a_s = a_s
#         self.Ry = Ry
#         self.pinchx = pinchx
#         self.pinchy = pinchy
#         self.dmg1 = dmg1
#         self.dmg2 = dmg2
#         self.beta = beta

#         # Initialisation
#         E = col.E          # Young modulus
#         Fy = col.Fy_web    # Yield strength of the column web (assume continous column)
#         dc = col.d         # Column depth
#         bf_c = col.bf      # Column flange width
#         tf_c = col.tf      # Column flange thickness
#         tp = col.tw        # Panel zone thickness
#         Ic = col.Iy        # Column moment of inertia (strong axis)
#         tpz = tp+t_dp      # Thickness of the panel zone and the doubler plate

#         # Trilinear Spring (Simple)
#         # Yield Shear
#         self.Vy = 0.55 * Fy * Ry * dc * tpz
#         # Shear Modulus
#         self.G = E/(2.0 * (1.0 + 0.30))
#         # Elastic Stiffness
#         self.Ke = 0.95 * self.G * tpz * dc
#         # Plastic Stiffness
#         self.Kp = 0.95 * self.G * bf_c * (tf_c * tf_c) / self.d_beam

#         # Define Trilinear Equivalent Rotational Spring (Simple)
#         # Yield point for Trilinear Spring at gamma1_y
#         self.gamma1_y = self.Vy / self.Ke
#         self.M1y = self.gamma1_y * (self.Ke * self.d_beam)
#         # Second Point for Trilinear Spring at 4 * gamma1_y
#         self.gamma2_y = 4.0 * self.gamma1_y
#         self.M2y = self.M1y + (self.Kp * self.d_beam) * (self.gamma2_y - self.gamma1_y)
#         # Third Point for Trilinear Spring at 100 * gamma1_y
#         self.gamma3_y = 100.0 * self.gamma1_y
#         self.M3y = self.M2y + (self.a_s * self.Ke * self.d_beam) * (self.gamma3_y - self.gamma2_y)


#         # Refined computation of the parameters for the backbone curve for the panel zone spring (Skiadopoulos et al. (2020))
#         # Panel Zone Elastic Stiffness
#         self.Ks_ref = tpz*(dc-tf_c)*self.G
#         self.Kb_ref = 12.0*E*(Ic+t_dp*(dc-2.0*tf_c)**3/12.0)/(self.d_beam-0)**2
#         self.Ke_ref = self.Ks_ref*self.Kb_ref/(self.Ks_ref+self.Kb_ref)

#         # Column Flange Stiffness
#         self.Ksf = 2.0*(tf_c*bf_c*self.G)
#         self.Kbf = 2.0*(12.0*E*bf_c*tf_c**3/12.0/(self.d_beam-0)**2)
#         self.Kf = self.Ksf*self.Kbf/(self.Ksf+self.Kbf)

#         # Kf/Ke Calculation for Panel Zone Categorization
#         self.Kf_Ke = self.Kf/self.Ke_ref

#         # Panel Zone Strength Coefficients (results from tests for a_w_eff and a_f_eff)
#         self.Cw1 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cw1_tests)
#         self.Cf1 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cf1_tests)
#         self.Cw4 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cw4_tests)
#         self.Cf4 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cf4_tests)
#         self.Cw6 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cw6_tests)
#         self.Cf6 = np.interp(self.Kf_Ke, Kf_Ke_tests, Cf6_tests)

#         # Panel Zone Model
#         self.V1 = Fy*Ry/math.sqrt(3)*(self.Cw1*(dc-tf_c)*tpz + self.Cf1*2*(bf_c-tp)*tf_c)
#         self.V4 = Fy*Ry/math.sqrt(3)*(self.Cw4*(dc-tf_c)*tpz + self.Cf4*2*(bf_c-tp)*tf_c)
#         self.V6 = Fy*Ry/math.sqrt(3)*(self.Cw6*(dc-tf_c)*tpz + self.Cf6*2*(bf_c-tp)*tf_c)

#         self.M1 = self.V1*(self.d_beam-tf_b)
#         self.M4 = self.V4*(self.d_beam-tf_b)
#         self.M6 = self.V6*(self.d_beam-tf_b)

#         self.Gamma_1 = self.V1/self.Ke_ref
#         self.Gamma_4 = 4*self.Gamma_1
#         self.Gamma_6 = 6*self.Gamma_1


#         # List of all indormation of the class (to be used for example in the SaveData/LoadData functions)
#         self.data = ["Data for PZ Hysteretic material model", 
#             self.ID,
#             self.col.NameTAG,
#             self.d_beam,
#             self.tf_b,
#             self.t_dp,
#             self.a_s,
#             self.Ry,
#             self.pinchx,
#             self.pinchy,
#             self.dmg1,
#             self.dmg2,
#             self.beta,
#             self.Vy,
#             self.G,
#             self.Ke,
#             self.Kp,
#             self.gamma1_y,
#             self.gamma2_y,
#             self.gamma3_y,
#             self.M1y,
#             self.M2y,
#             self.M3y,
#             self.Ks_ref,             # refined parameters
#             self.Kb_ref,
#             self.Ke_ref,
#             self.Ksf,
#             self.Kbf,
#             self.Kf,
#             self.Kf_Ke,
#             self.Cw1,
#             self.Cf1,
#             self.Cw4,
#             self.Cf4,
#             self.Cw6,
#             self.Cf6,
#             self.V1,
#             self.V4,
#             self.V6,
#             self.M1,
#             self.M4,
#             self.M6,
#             self.Gamma_1,
#             self.Gamma_4,
#             self.Gamma_6]

#     # Methods
#     def Hysteretic(self,  plot = False, block = False):
#         # Hysteretic Material
#         # $matTag       integer tag identifying material
#         # $s1p $e1p     stress and strain (or force & deformation) at first point of the envelope in the positive direction
#         # $s2p $e2p     stress and strain (or force & deformation) at second point of the envelope in the positive direction
#         # $s3p $e3p     stress and strain (or force & deformation) at third point of the envelope in the positive direction (optional)
#         # $s1n $e1n     stress and strain (or force & deformation) at first point of the envelope in the negative direction
#         # $s2n $e2n     stress and strain (or force & deformation) at second point of the envelope in the negative direction
#         # $s3n $e3n     stress and strain (or force & deformation) at third point of the envelope in the negative direction (optional)
#         # $pinchx       pinching factor for strain (or deformation) during reloading
#         # $pinchy       pinching factor for stress (or force) during reloading
#         # $damage1      damage due to ductility: D1(mu-1)
#         # $damage2      damage due to energy: D2(Eii/Eult)
#         # $beta         power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0) 

#         uniaxialMaterial("Hysteretic", self.ID,
#             self.M1y, self.gamma1_y, self.M2y, self.gamma2_y, self.M3y, self.gamma3_y,
#             -self.M1y, -self.gamma1_y, -self.M2y, -self.gamma2_y, -self.M3y, -self.gamma3_y,
#             self.pinchx, self.pinchy, self.dmg1, self.dmg2, self.beta)


#         if plot:
#             # Data for plotting
#             # Last point for plot
#             gamma3_y_plot = 10.0 * self.gamma1_y
#             M3y_plot = self.M2y + (self.a_s * self.Ke * self.d_beam) * (gamma3_y_plot - self.gamma2_y)

#             x_axis = np.array([0.0, self.gamma1_y, self.gamma2_y, gamma3_y_plot])
#             y_axis = ([0.0, self.M1y, self.M2y, M3y_plot])

#             fig, ax = plt.subplots()
#             ax.plot(x_axis, y_axis, 'k-')

#             ax.set(xlabel='Rotation [rad]', ylabel='Moment [kNmm]', 
#                 title='Backbone curve for trilinear PZ spring model for material ID={}'.format(self.ID))
#             ax.grid()

#             print("")
#             print("Trilinear PZ Spring Material Model, ID = {}".format(self.ID))
#             print("gamma1_y = {}".format(self.gamma1_y))
#             print("gamma2_y = {}".format(self.gamma2_y))
#             print("gamma3_y = {}".format(self.gamma3_y))
#             print("M1y = {} kNm".format(self.M1y/1000))
#             print("M2y = {} kNm".format(self.M2y/1000))
#             print("M3y = {} kNm".format(self.M3y/1000))
#             print("")

#             if block:
#                 plt.show()

#     def RefinedHysteretic(self,  plot = False, block = False):
#         # Refined backbone curve for the panel zone spring (Skiadopoulos et al. (2020))

#         # Hysteretic Material
#         # $matTag       integer tag identifying material
#         # $s1p $e1p     stress and strain (or force & deformation) at first point of the envelope in the positive direction
#         # $s2p $e2p     stress and strain (or force & deformation) at second point of the envelope in the positive direction
#         # $s3p $e3p     stress and strain (or force & deformation) at third point of the envelope in the positive direction (optional)
#         # $s1n $e1n     stress and strain (or force & deformation) at first point of the envelope in the negative direction
#         # $s2n $e2n     stress and strain (or force & deformation) at second point of the envelope in the negative direction
#         # $s3n $e3n     stress and strain (or force & deformation) at third point of the envelope in the negative direction (optional)
#         # $pinchx       pinching factor for strain (or deformation) during reloading
#         # $pinchy       pinching factor for stress (or force) during reloading
#         # $damage1      damage due to ductility: D1(mu-1)
#         # $damage2      damage due to energy: D2(Eii/Eult)
#         # $beta         power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0) 

#         uniaxialMaterial("Hysteretic", self.ID,
#             self.M1, self.Gamma_1, self.M4, self.Gamma_4, self.M6, self.Gamma_6,
#             -self.M1, -self.Gamma_1, -self.M4, -self.Gamma_4, -self.M6, -self.Gamma_6,
#             self.pinchx, self.pinchy, self.dmg1, self.dmg2, self.beta)


#         if plot:
#             # Data for plotting
#             x_axis = np.array([0.0, self.Gamma_1, self.Gamma_4, self.Gamma_6])
#             y_axis = ([0.0, self.M1, self.M4, self.M6])

#             fig, ax = plt.subplots()
#             ax.plot(x_axis, y_axis, 'k-')

#             ax.set(xlabel='Rotation [rad]', ylabel='Moment [kNmm]', 
#                 title='Backbone curve for refined trilinear PZ spring model for material ID={}'.format(self.ID))
#             ax.grid()

#             print("")
#             print("Refined trilinear PZ Spring Material Model, ID = {}".format(self.ID))
#             print("gamma1 = {}".format(self.Gamma_1))
#             print("gamma4 = {}".format(self.Gamma_4))
#             print("gamma6 = {}".format(self.Gamma_6))
#             print("M1 = {} kNm".format(self.M1/1000))
#             print("M4 = {} kNm".format(self.M4/1000))
#             print("M6 = {} kNm".format(self.M6/1000))
#             print("")

#             if block:
#                 plt.show()



# class Concrete04MaterialModel(DataManagement):
#     # Class that stores funcions and material properties of a rectangular shape RC  profile. For more information see Mander et Al. 1988
#     # Warning: the units should be mm and kN
#       #TODO: warning if concrete fc is bigger than something (see article Lee)
    
#     def __init__(self, ConfinedID, UnconfinedID, ele: SectionRCRectShape, N_G = 0, Cantilever = False):
#         """
#         Parameters
#         ----------
#         ConfinedID : int
#             The ID of the material model for confined
#         UnconfinedID : int
#             The ID of the material model for unconfined
#         ele : SectionRCRectShape
#             The element section
#         N_G : double
#             The gravity axial load (default: 0)
#         Cantilever : bool
#             A parameter to know if the element is a cantilever or not (double fixed) (default: False)
#         """
#         self.ConfinedID = ConfinedID 
#         self.UnconfinedID = UnconfinedID
#         self.ele = ele
#         self.N_G = N_G
#         if Cantilever:
#             L_o = ele.L
#             K_factor = 3
#         else:
#             L_o = ele.L/2
#             K_factor = 6
#         self.L_o = L_o
#         self.K_factor = K_factor

#         # Initialization
#         self.ec = -0.002
#         self.ecp = 2.0*self.ec                                  # concrete spalling
#         self.esu = 0.05                                         # strain capacity of stirrups 0.012-0.05 experimental results (BOOK BEYER)
#         self.ecu = -0.004 + (1.4*0*self.esu*ele.fs) / ele.fc    # FROM BRIDGE BOOK OF KATRIN BEYER!
#         self.fct = 0.30 * math.pow(-ele.fc, 2/3)                # SIA 262 2013
#         self.et = self.fct/ele.Ec                               # Mander Eq 45
#         wx_top : numpy.ndarray
#             A vector that defines the distance between bars in x direction (NOT CENTERLINE DISTANCE). One range of bars implemented
#         wx_bottom : numpy.ndarray
#             A vector that defines the distance between bars in x direction (NOT CENTERLINE DISTANCE). One range of bars implemented
#         wy : numpy.ndarray
#             A vector that defines the distance between bars in y direction (NOT CENTERLINE DISTANCE). One range of bars implemented


#         # Confined
#         self.k1 = 4.1                                           # 4.1 from Richart et al. 1928 or 5.6 from Balmer 1949 
#         self.k2 = 5.0*self.k1
#         self.ineffectual_area = (np.sum(np.multiply(ele.wy, ele.wy)) + np.sum(np.multiply(ele.wx, ele.wx)))*2.0/6.0
#         self.Ae = (ele.Ac - self.ineffectual_area) * (1.0 - (ele.s-ele.D_hoops)/2.0/ele.bc)*(1.0 - (ele.s-ele.D_hoops)/2.0/ele.dc)
#         self.rho_cc = ele.nr_bars*ele.D_bars**2/4.0*math.pi / ele.Ac
#         self.Acc = ele.Ac*(1.0-self.rho_cc)
#         self.ke = self.Ae/self.Acc
#         self.fl_x = -ele.rho_s_x * ele.fs
#         self.fl_y = -ele.rho_s_y * ele.fs
#         self.fl_prime = (self.fl_x + self.fl_y)/2.0 * self.ke # derive the interpolating curve
#         self.confinement_factor = -1.254 + 2.254 * math.sqrt(1.0+7.94*self.fl_prime/ele.fc) - 2.0*self.fl_prime/ele.fc # in Mander, it has a prime
#         self.fcc = ele.fc * self.confinement_factor
#         self.ecc = (1.0 + 5.0 * (self.confinement_factor-1.0)) * self.ec
#         self.eccu = -0.004 + (1.4*(ele.rho_s_x+ele.rho_s_y)*self.esu*ele.fs) / self.fcc # FROM BRIDGE BOOK OF KATRIN BEYER!
        

#         # List of all indormation of the class (to be used for example in the SaveData/LoadData functions)
#         self.data = ["Data for Concrete04 material model", 
#             self.ConfinedID,
#             self.UnconfinedID,
#             self.ele.NameTAG,
#             self.ec,
#             self.ecp,
#             self.esu,
#             self.ecu,
#             self.k1,
#             self.k2,
#             self.ineffectual_area,
#             self.Ae,
#             self.rho_cc,
#             self.Acc,
#             self.ke,
#             self.fl_x,
#             self.fl_y,
#             self.fl_prime,
#             self.confinement_factor,
#             self.fcc,
#             self.ecc,
#             self.eccu,
#             self.N_G, 
#             self.L_o,
#             self.K_factor]



#     # Methods
#     def Confined(self, plot = False, block = False):
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

#         uniaxialMaterial("Concrete04", self.ConfinedID, self.fcc, self.ecc, self.eccu, self.ele.Ec, self.fct, self.et, 0.1)

#         if plot:
#             fig, ax = plt.subplots()
#             self.PlotConcrete04("C", ax)

#             print("")
#             print("Concrete04 Material Model (Confined), ID = {}".format(self.ConfinedID))
#             print('Stress max fcc = {} kN/mm2'.format(self.fcc))
#             print('Strain with max stress ecc = {}'.format(self.ecc))
#             print('Strain max eccu = {}'.format(self.eccu))
#             print("")

#             if block:
#                 plt.show()

#     def Unconfined(self, plot = False, block = False):
#         # Generate the material model using the given parameters

#         # Define Concrete04 Popovics Concrete material model for unconfined concrete
#         # uniaxialMaterial("Concrete04", matTag, fc, ec, ecu, Ec, <fct et> <beta>)
#         # matTag     integer tag identifying material
#         # fc     floating point values defining concrete compressive strength at 28 days (compression is negative)*
#         # ec     floating point values defining concrete strain at maximum strength*
#         # ecu    floating point values defining concrete strain at crushing strength*
#         # Ec     floating point values defining initial stiffness**
#         # fct    floating point value defining the maximum tensile strength of concrete
#         # et     floating point value defining ultimate tensile strain of concrete
#         # beta   loating point value defining the exponential curve parameter to define the residual stress (as a factor of ft) at etu 
        
#         uniaxialMaterial("Concrete04", self.UnconfinedID, self.ele.fc, self.ec, self.ecu, self.ele.Ec)

#         if plot:
#             fig, ax = plt.subplots()
#             self.PlotConcrete04("U", ax)

#             print("")
#             print("Concrete04 Material Model (Unconfined), ID = {}".format(self.UnconfinedID))
#             print('Stress max fc = {} kN/mm2'.format(self.ele.fc))
#             print('Strain with max stress ec = {}'.format(self.ec))
#             print('Strain max ecu = {}'.format(self.ecu))
#             print("")

#             if block:
#                 plt.show()



#     def Concrete04Funct(self, fcc, ec, ecc, Ec):
#         x = ec/ecc
#         r = Ec / (Ec - fcc/ecc)
#         return fcc*x*r / (r-1+x**r)


#     def PlotConcrete04(self, Type, ax):
#         # For ax, just use:
#         # fig, ax = plt.subplots()
#         # to generate the figure
#         # Type is "C" or "U" for confined or unconfined

#         if Type == "C":
#             fcc = self.fcc
#             ecc = self.ecc
#             eccu = self.eccu
#             name = "Confined"
#         elif Type == "U":
#             fcc = self.ele.fc
#             ecc = self.ec
#             eccu = self.ecu
#             name = "Unconfined"
#         else:
#             raise NameError("Type should be C or U for Confined and Unconfined in material model (ConfinedID={})".format(self.ConfinedID))
#         # Data for plotting
#         N = int(-eccu*1e5)
#         x_axis = np.zeros(N)
#         y_axis = np.zeros(N)
#         for i in range(N):
#             x_axis[i] = -i/1e5
#             y_axis[i] = self.Concrete04Funct(fcc, x_axis[i], ecc, self.ele.Ec)


#         ax.plot(-x_axis*100.0, -y_axis*1000.0, 'k-', label = name)
#         ax.set(xlabel='Strain [%]', ylabel='Stress MPa]', 
#                         title='Backbone curve for Concrete04 material model')
#         plt.grid()
        




# class Steel01forRCMaterialModel(DataManagement):
#     # Class that stores funcions and material properties longitudinal bars for a rectangular shape RC section.
#     # Warning: the units should be mm and kN
    
#     def __init__(self, ID, ele: SectionRCRectShape, b = 0.01):
#         """
#         Parameters
#         ----------
#         ID : int
#             The ID of the material model
#         ele : Section
#             The element section
#         b : double
#             Strain hardening factor (default: 0.01)
#         """
#         self.ID = ID
#         self.ele = ele
#         self.b = b

#         # List of all indormation of the class (to be used for example in the SaveData/LoadData functions)
#         self.data = ["Data for Steel01 material model", 
#             self.ID, 
#             self.ele.NameTAG,
#             self.b]



#     # Methods
#     def Steel01(self, plot = False, block = False):
#         # Generate the material model using the given parameters

#         # Define Steel01 material model
#         # uniaxialMaterial("Steel01", ID, fy, Ey, b)
        
#         uniaxialMaterial("Steel01", self.ID, self.ele.fy, self.ele.Ey, self.b)

#         if plot:
#             # Data for plotting
#             ey = self.ele.fy/self.ele.Ey
#             e_pl = 10.0 * ey # to show that if continues with this slope
#             sigma_pl = self.b * self.ele.Ey * e_pl

#             x_axis = ([0.0, ey, ey+e_pl])
#             y_axis = ([0.0, self.ele.fy, self.ele.fy+sigma_pl])

#             fig, ax = plt.subplots()
#             ax.plot(x_axis, y_axis, 'k-')

#             ax.set(xlabel='Strain [%]', ylabel='Stress [kN]', 
#                 title='Backbone curve for Steel01 model for material ID={}'.format(self.ID))
#             ax.grid()

#             print("")
#             print("Steel01 Material Model, ID = {}".format(self.ID))
#             print('Max elastic strain eps y = {}'.format(ey))
#             print('Yielding stress fy = {} kN/mm2'.format(self.ele.fy))
#             print('Young modulus = {} kN/mm2'.format(self.ele.Ey))
#             print('Hardening factor = {}'.format(self.b))
#             print("")

#             if block:
#                 plt.show()

