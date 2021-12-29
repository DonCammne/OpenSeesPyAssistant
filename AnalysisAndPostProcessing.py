# Module with premade analysis and postprocessing functions
#   Carmine Schipani, 2021

from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import openseespy.postprocessing.Get_Rendering as opsplt
from OpenSeesPyHelper.ErrorHandling import *
from OpenSeesPyHelper.Units import *
from OpenSeesPyHelper.Constants import *
from OpenSeesPyHelper.FunctionalFeatures import *

"""Module with pre-made analysis and postprocessing functions.
"""

class Analysis():
    def __init__(self, data_dir: str,  name_ODB: str, algo = "KrylovNewton", test_type = "NormDispIncr", test_opt = 0, max_iter = MAX_ITER, tol = TOL, allow_smaller_step = False):
        if max_iter < 0: raise NegativeValue()
        if tol < 0: raise NegativeValue()

        self.data_dir = data_dir
        self.name_ODB = name_ODB
        self.algo = algo
        self.test_type = test_type
        self.tol = tol
        self.test_opt = test_opt
        self.max_iter = max_iter
        self.allow_smaller_step = allow_smaller_step


    def Gravity(self, loaded_nodes: list, Fy: list, timeSeries_ID: int, pattern_ID: int, n_step = 10, timeSeries_type = "Linear", pattern_type = "Plain",
        constraints_type = "Plain", numberer_type = "RCM", system_type = "BandGeneral", analysis_type = "Static", show_plot = False):
        # To save the analysis output for deformed shape, use createODB command before running the analysis
        # The following command saves the model data, and output for pushover analysis and the first nEigen modes 
        # in a folder "NAME_ODB"
        # first node of loaded_nodes is recorded
        if len(loaded_nodes) != len(Fy): raise WrongDimension()
        if timeSeries_ID < 1: raise NegativeValue()
        if pattern_ID < 1: raise NegativeValue()

        # for mass defined: opsplt.createODB(self.name_ODB, "Gravity", Nmodes = nEigen); 
        # opsplt.createODB(self.name_ODB, "Gravity");

        # Create load pattern
        timeSeries(timeSeries_type, timeSeries_ID)
        pattern(pattern_type, timeSeries_ID, pattern_ID)
        for ii, node_ID in enumerate(loaded_nodes):
            load(node_ID, 0.0, Fy[ii], 0.0)    # load(IDNode, Fx, Fy, Mz)
        DGravity = 1.0/n_step               # load increment                  

        # Set up analysis options
        constraints(constraints_type)   # how it handles boundary conditions
        numberer(numberer_type)         # renumber dof's to minimize band-width (optimization)
        system(system_type)             # how to store and solve the system of equations in the analysis
                                        # For static model, BandGeneral, for transient and/or big model, UmfPack
        integrator("LoadControl", DGravity) # LoadControl and DisplacementControl only with static model, linear TimeSeries w/ factor of 1
                                            # Newmark used for transient model  
        algorithm("Newton")             # placeholder
        analysis(analysis_type)         # define type of analysis: static for pushover

        # Analysis
        dataG = np.zeros((n_step+1,2))
        next_step = 0
        print("Gravity analysis starts")
        for iteration in range(n_step):
            # next_step = ProgressingPercentage(n_step, iteration, next_step)
            convergence = self.__VertLoadCtrlLoop(DGravity, iteration,
                self.algo, self.test_type, self.tol, self.test_opt, self.max_iter)
            if convergence != 0: break
            dataG[iteration+1,0] = nodeDisp(loaded_nodes[0], 2)/mm_unit
            dataG[iteration+1,1] = getLoadFactor(pattern_ID)*Fy[0]/kN_unit

        if show_plot:
            plt.plot(dataG[:,0], dataG[:,1])
            plt.xlabel('Vertical Displacement [mm]')
            plt.ylabel('Vertical Load [kN]')
            plt.title('Gravity curve')
            plt.show()

        loadConst("-time", 0.0)

        print("Gravity complete")


    def Pushover(self, CtrlNode: int, Dmax, Dincr, timeSeries_ID: int, pattern_ID: int, Fx = 1*kN_unit, fiber_ID_analysed = -1, fiber_section = 1,
        timeSeries_type = "Linear", pattern_type = "Plain", constraints_type = "Plain", numberer_type = "RCM", system_type = "UmfPack", analysis_type = "Static",
        show_plot = True):
        # To save the analysis output for deformed shape, use createODB command before running the analysis
        # The following command saves the model data, and output for pushover analysis and the first nEigen modes 
        # in a folder "NAME_ODB"
        # fiber:section: section number, i.e. the Gauss integratio number
        # CtrlNode = 12;                  # node where disp is read for disp control (w/ PZ use middle left xy02)
        # CtrlDOF = 1;                    # degree of freedom read for disp control (1 = x displacement)
        # Dmax = 0.01*L                   # maximum displacement of pushover: 1% roof drift
        # Dincr = 0.05                    # displacement increment (0.01-0.05)
        if CtrlNode < 1: raise NegativeValue()
        if timeSeries_ID < 1: raise NegativeValue()
        if pattern_ID < 1: raise NegativeValue()
        if fiber_ID_analysed != -1 and fiber_ID_analysed < 1: raise NegativeValue()

        # for mass defined: opsplt.createODB(self.name_ODB, "Pushover", Nmodes = nEigen); 
        opsplt.createODB(self.name_ODB, "Pushover");
        if fiber_ID_analysed != -1: opsplt.saveFiberData2D(self.name_ODB, "Pushover", fiber_ID_analysed, fiber_section)


        # Create load pattern
        timeSeries(timeSeries_type, timeSeries_ID)
        pattern(pattern_type, timeSeries_ID, pattern_ID)
        load(CtrlNode, Fx, 0.0, 0.0)    # load(IDNode, Fx, Fy, Mz)
        Nsteps = int(abs(Dmax/Dincr))   # number of pushover analysis steps
        
        # Set up analysis options
        constraints(constraints_type)   # how it handles boundary conditions
        numberer(numberer_type)         # renumber dof's to minimize band-width (optimization)
        system(system_type)             # how to store and solve the system of equations in the analysis
                                        # For static model, BandGeneral, for transient and/or big model, UmfPack
        analysis(analysis_type)         # define type of analysis: static for pushover

        # Analysis
        dataPO = np.zeros((Nsteps+1,2))
        next_step = 0
        print("Pushover analysis starts")
        for iteration in range(Nsteps):
            next_step = ProgressingPercentage(Nsteps, iteration, next_step)
            convergence = self.__LatDispCtrlLoop(CtrlNode, Dincr, iteration,
                self.algo, self.test_type, self.tol, self.test_opt, self.max_iter, self.allow_smaller_step)
            if convergence != 0: break
            dataPO[iteration+1,0] = nodeDisp(CtrlNode, 1)/mm_unit
            dataPO[iteration+1,1] = getLoadFactor(pattern_ID)*Fx/kN_unit
        
        if show_plot:
            plt.plot(dataPO[:,0], dataPO[:,1])
            plt.xlabel('Horizontal Displacement [mm]')
            plt.ylabel('Horizontal Load [kN]')
            plt.title('Pushover curve')
            plt.show()
        
        print("Pushover complete")

        wipe()


    def LoadingProtocol(self, CtrlNode: int, discr_LP: np.ndarray, timeSeries_ID: int, pattern_ID: int, Fx = 1*kN_unit, fiber_ID_analysed = -1, fiber_section = 1,
        timeSeries_type = "Linear", pattern_type = "Plain", constraints_type = "Plain", numberer_type = "RCM", system_type = "UmfPack", analysis_type = "Static",
        show_plot = True):
        # To save the analysis output for deformed shape, use createODB command before running the analysis
        # The following command saves the model data, and output for pushover analysis and the first nEigen modes 
        # in a folder "NAME_ODB"
        # fiber:section: section number, i.e. the Gauss integratio number
        # CtrlNode = 12;                  # node where disp is read for disp control (w/ PZ use middle left xy02)
        # CtrlDOF = 1;                    # degree of freedom read for disp control (1 = x displacement)
        # Dmax = 0.01*L                   # maximum displacement of pushover: 1% roof drift
        # Dincr = 0.05                    # displacement increment (0.01-0.05)
        if CtrlNode < 1: raise NegativeValue()
        if timeSeries_ID < 1: raise NegativeValue()
        if pattern_ID < 1: raise NegativeValue()
        if fiber_ID_analysed != -1 and fiber_ID_analysed < 1: raise NegativeValue()

        # for mass defined: opsplt.createODB(self.name_ODB, "LoadingProtocol", Nmodes = nEigen); 
        opsplt.createODB(self.name_ODB, "LoadingProtocol");
        if fiber_ID_analysed != -1: opsplt.saveFiberData2D(self.name_ODB, "LoadingProtocol", fiber_ID_analysed, fiber_section)

        # Create load pattern
        timeSeries(timeSeries_type, timeSeries_ID)
        pattern(pattern_type, timeSeries_ID, pattern_ID)
        load(CtrlNode, Fx, 0.0, 0.0)    # load(IDNode, Fx, Fy, Mz)
        dU_prev = 0
        Nsteps = np.size(discr_LP)      # number of pushover analysis steps

        # Set up analysis options
        constraints(constraints_type)   # how it handles boundary conditions
        numberer(numberer_type)         # renumber dof's to minimize band-width (optimization)
        system(system_type)             # how to store and solve the system of equations in the analysis
                                        # For static model, BandGeneral, for transient and/or big model, UmfPack
        analysis(analysis_type)         # define type of analysis: static for LoadingProtocol

        # Analysis
        dataLP = np.zeros((Nsteps+1,2))
        next_step = 0
        print("Loading Protocol analysis starts")
        for iteration in range(Nsteps):
            # Compute displacement usinf the given loading protocol (discretized)
            dU_next = discr_LP[iteration]
            dU = dU_next - dU_prev
            dU_prev = dU_next
            
            next_step = ProgressingPercentage(Nsteps, iteration, next_step)
            convergence = self.__LatDispCtrlLoop(CtrlNode, dU, iteration, 
                self.algo, self.test_type, self.tol, self.test_opt, self.max_iter, self.allow_smaller_step)
            if convergence != 0: break
            dataLP[iteration+1,0] = nodeDisp(CtrlNode, 1)/mm_unit
            dataLP[iteration+1,1] = getLoadFactor(pattern_ID)*Fx/kN_unit
        
        if show_plot:
            plt.plot(dataLP[:,0], dataLP[:,1])
            plt.xlabel('Horizontal Displacement [mm]')
            plt.ylabel('Horizontal Load [kN]')
            plt.title('Loading Protocol curve')
            plt.show()
        
        print("Loading Protocol complete")
        wipe()


    def __VertLoadCtrlLoop(self, DGravity, iteration, algo = "KrylovNewton", test_type = "NormDispIncr", tol = TOL, test_opt = 0, max_iter = MAX_ITER):
        # Default analysis
        integrator("LoadControl", DGravity)         # LoadControl and DisplacementControl only with static model, linear TimeSeries w/ factor of 1
                                                    # Newmark used for transient model 
        test(test_type, tol, max_iter, test_opt)    # type of convergence criteria with tolerance, max iterations;
                                                    # Normally use EnergyIncr, if conv issues, try NormDispIncr; optional: test_opt = 2 for debugging  
        algorithm(algo)                             # use Newton's solution algorithm: updates tangent stiffness at every iteration
        convergence = analyze(1)                    # this will return zero if no convergence problems were encountered
        
        # Convergence analysis
        if convergence != 0:
            print("---------------------------------------------------------------------------------------")
            print("Vertical Load-control analysis failed at iteration {}".format(iteration))
            print("Convergence analysis starts")
            print("---------------------------------------------------------------------------------------")
            
            if convergence != 0:
                print("Try 1")
                convergence = _ConvergenceTest("KrylovNewton", "NormDispIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 2")
                convergence = _ConvergenceTest("KrylovNewton", "NormDispIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 3")
                convergence = _ConvergenceTest("KrylovNewton", "EnergyIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 4")
                convergence = _ConvergenceTest("KrylovNewton", "EnergyIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 5")
                convergence = _ConvergenceTest("Newton", "NormDispIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 6")
                convergence = _ConvergenceTest("Newton", "NormDispIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 7")
                convergence = _ConvergenceTest("Newton", "EnergyIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 8")
                convergence = _ConvergenceTest("Newton", "EnergyIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 9")
                convergence = _ConvergenceTest("ModifiedNewton", "NormDispIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 10")
                convergence = _ConvergenceTest("ModifiedNewton", "NormDispIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 11")
                convergence = _ConvergenceTest("ModifiedNewton", "EnergyIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 12")
                convergence = _ConvergenceTest("ModifiedNewton", "EnergyIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("")
                print("#############################################################################")
                print("NO CONVERGENCE! Vertical Load-control analysis stops at iteration {} ".format(iteration))
                print("#############################################################################")
                print("")
            else:
                print("Convergence reached, convergence analysis ends.")

        return convergence


    def __LatDispCtrlLoop(self, CtrlNode, dU, iteration, algo = "KrylovNewton", test_type = "NormDispIncr", tol = TOL, test_opt = 0, max_iter = MAX_ITER, allow_smaller_step = False):
        # Default analysis
        CtrlDOF = 1
        integrator("DisplacementControl", CtrlNode, CtrlDOF, dU, 1, dU, dU)  # use displacement-controlled analysis
        test(test_type, tol, max_iter, test_opt)    # type of convergence criteria with tolerance, max iterations;
                                                    # Normally use EnergyIncr, if conv issues, try NormDispIncr; optional: test_opt = 2 for debugging  
        algorithm(algo)                             # use Newton's solution algorithm: updates tangent stiffness at every iteration
        convergence = analyze(1)                    # this will return zero if no convergence problems were encountered
        
        # Convergence analysis
        if convergence != 0:
            print("---------------------------------------------------------------------------------------")
            print("Lateral Displacement-control analysis failed at iteration {} and control node lateral displacement = {} mm.".format(iteration, nodeDisp(CtrlNode, CtrlDOF)/mm_unit))
            print("Convergence analysis starts")
            print("---------------------------------------------------------------------------------------")
            
            if convergence != 0:
                print("Try 1")
                convergence = _ConvergenceTest("KrylovNewton", "NormDispIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 2")
                convergence = _ConvergenceTest("KrylovNewton", "NormDispIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 3")
                convergence = _ConvergenceTest("KrylovNewton", "EnergyIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 4")
                convergence = _ConvergenceTest("KrylovNewton", "EnergyIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 5")
                convergence = _ConvergenceTest("Newton", "NormDispIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 6")
                convergence = _ConvergenceTest("Newton", "NormDispIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 7")
                convergence = _ConvergenceTest("Newton", "EnergyIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 8")
                convergence = _ConvergenceTest("Newton", "EnergyIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 9")
                convergence = _ConvergenceTest("ModifiedNewton", "NormDispIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 10")
                convergence = _ConvergenceTest("ModifiedNewton", "NormDispIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0:
                print("Try 11")
                convergence = _ConvergenceTest("ModifiedNewton", "EnergyIncr", tol, max_iter, test_opt)
            
            if convergence != 0:
                print("Try 12")
                convergence = _ConvergenceTest("ModifiedNewton", "EnergyIncr", tol*10, max_iter*10, test_opt)
            
            if convergence != 0 and allow_smaller_step:
                print("Use smaller steps")
                print("10 times more intergator iteration, min_dU = dU/10")
                integrator("DisplacementControl", CtrlNode, CtrlDOF, dU, 10, dU/10, dU)

                if convergence != 0:
                    print("Try 13")
                    convergence = _ConvergenceTest("KrylovNewton", "NormDispIncr", tol, max_iter, test_opt)

                if convergence != 0:
                    print("Try 14")
                    convergence = _ConvergenceTest("KrylovNewton", "EnergyIncr", tol, max_iter, test_opt)

                if convergence != 0:
                    print("Try 15")
                    convergence = _ConvergenceTest("Newton", "NormDispIncr", tol, max_iter, test_opt)
                
                if convergence != 0:
                    print("Try 16")
                    convergence = _ConvergenceTest("Newton", "EnergyIncr", tol, max_iter, test_opt)

                if convergence != 0:
                    print("Try 17")
                    convergence = _ConvergenceTest("ModifiedNewton", "NormDispIncr", tol, max_iter, test_opt)
                
                if convergence != 0:
                    print("Try 19")
                    convergence = _ConvergenceTest("ModifiedNewton", "EnergyIncr", tol, max_iter, test_opt)
                
            if convergence != 0:
                print("")
                print("#############################################################################")
                print("NO CONVERGENCE! Lateral Displacement-control analysis stops at iteration {} and control node lateral displacement = {} mm.".format(iteration, nodeDisp(CtrlNode, CtrlDOF)/mm_unit))
                print("#############################################################################")
                print("")
            else:
                print("Convergence reached, convergence analysis ends.")

        return convergence


    def DeformedShape(self, load_case, scale = 1, animate = False, dt = 0.01):
        # Display deformed shape, the scaling factor needs to be adjusted for each model
        opsplt.plot_deformedshape(Model = self.name_ODB, LoadCase=load_case, scale = scale)
        if animate:
            opsplt.animate_deformedshape(Model = self.name_ODB, LoadCase=load_case, dt = dt, scale = scale)


    def FiberResponse(self, load_case, fiber_ID_analysed, fiber_section = 1, animate_stress = False, animate_strain = False, fps = 25):
        opsplt.plot_fiberResponse2D(self.name_ODB, load_case, fiber_ID_analysed, fiber_section, InputType = 'stress')
        if animate_stress:
            ani1 = opsplt.animate_fiberResponse2D(self.name_ODB, load_case, fiber_ID_analysed, fiber_section, InputType = 'stress', fps = fps)
        if animate_strain:
            ani1 = opsplt.animate_fiberResponse2D(self.name_ODB, load_case, fiber_ID_analysed, fiber_section, InputType = 'strain', fps = fps)


def _ConvergenceTest(algo_type, test_type, tol, max_iter, test_opt):
    # private module function
    print("algorithm: {}".format(algo_type))
    print("test: {}, tol = {}, max iter = {}".format(test_type, tol, max_iter))
    test(test_type, tol, max_iter, test_opt)
    algorithm(algo_type)
    return analyze(1)



