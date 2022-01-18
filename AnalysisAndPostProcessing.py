"""Module with pre-made analysis and postprocessing functions. \n
Carmine Schipani, 2021
"""

from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import os
import openseespy.postprocessing.Get_Rendering as opsplt
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *
from OpenSeesPyAssistant.Constants import *
from OpenSeesPyAssistant.FunctionalFeatures import *


class Analysis():
    """Class dedicated to the analysis of the OpenSeesPy model. The Gravity method should be run first to perform the Load-control analysis (apply the vertical load). If no vertical load, this method can be omitted. \n
    Then only one of the Displacement-control (Pushover or LoadingProtocol) or Load-control (LateralForce) analysis can ran. \n
    After the analysis reach convergence in the last step, for the postprocessing, the DeformedShape method can be used to see the final deformed shape and the animation of the entire loading protocol;
    the FiberResponse method can be used to see the animation of the same fiber section recorded during the analysis (strain and/or stress).
    """
    def __init__(self, data_dir: str,  name_ODB: str, algo = "KrylovNewton", test_type = "NormDispIncr", test_opt = 0, max_iter = MAX_ITER, tol = TOL, allow_smaller_step = False):
        """
        The constructor of the class.

        @param data_dir (str): Directory in which the results from the analysis will be stored. Use the recorders (from OpenSeesPy) or the Record method from MemberModel.
        @param name_ODB (str): Name for the folder in which the data for the animations and the fibers are stored.
        @param algo (str, optional): Type of alghoritm chosen for the analysis. It detemines how to construct a SolutionAlgorithm object, which determines the sequence of steps taken to solve the non-linear equation.
            For more information on the available types, see the OpenSeesPy documentation. Defaults to "KrylovNewton".
        @param test_type (str, optional): Type of test chosen for the analysis. It determines how to construct a ConvergenceTest object.
            Certain SolutionAlgorithm objects require a ConvergenceTest object to determine if convergence has been achieved at the end of an iteration step.
            For more information on the available types, see the OpenSeesPy documentation. Defaults to "NormDispIncr".
        @param test_opt (int, optional): Print-flag from 0 to 5 used to receive more info during the iteration
            (for example: 0 print nothing and 2 print information on norms and number of iterations at end of successful test).
            For more information, see the OpenSeesPy documentation. Defaults to 0.
        @param max_iter (float, optional): Maximal number of iterations to check. Defaults to MAX_ITER (from Constants Module).
        @param tol (float, optional): Tolerance criteria used to check for convergence. Defaults to TOL (from Constants Module).
        @param allow_smaller_step (bool, optional): Allow smaller steps in the displacement-control analysis. Defaults to False.

        @exception NegativeValue: The argument max_iter should be positive.
        @exception NegativeValue: The argument tol should be positive.
        """
        if max_iter < 0: raise NegativeValue()
        if tol < 0: raise NegativeValue()
        if not os.path.exists(data_dir):
            print("Folder {} not found in this directory; creating one".format(data_dir))
            os.makedirs(data_dir)

        self.data_dir = data_dir
        self.name_ODB = name_ODB
        self.algo = algo
        self.test_type = test_type
        self.tol = tol
        self.test_opt = test_opt
        self.max_iter = max_iter
        self.allow_smaller_step = allow_smaller_step
        self.load_case = "None"


    def Gravity(self, loaded_nodes: list, Fy: list, timeSeries_ID: int, pattern_ID: int, n_step = 10, timeSeries_type = "Linear", pattern_type = "Plain",
        constraints_type = "Plain", numberer_type = "RCM", system_type = "BandGeneral", analysis_type = "Static", show_plot = False):
        """
        Method to perform the gravity analyisis with vertical loadings (load-control).
        It can be used before calling the Pushover or LoadingProtocol methods that perform the actual anlysis. If no vertical loadings present, this method can be avoided.

        @param loaded_nodes (list): List of nodes that are loaded by the the forces in Fy. The first node will be recorded (thus usually should be in the roof).
        @param Fy (list): List of vertical loadings (negative is toward the ground, thus compression; see global coordinate system).
        @param timeSeries_ID (int): ID of the timeseries.
        @param pattern_ID (int): ID of the pattern.
        @param n_step (int, optional): Number of steps used to during the analysis to reach the objective state (with 100% vertical loadings imposed). Defaults to 10.
        @param timeSeries_type (str, optional): Type of timeseries chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Linear".
        @param pattern_type (str, optional): Type of pattern chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param constraints_type (str, optional): Type of contraints chosen. It detemines how the constraint equations are enforced in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param numberer_type (str, optional): Type of numberer chosen. It determines the mapping between equation numbers and degrees-of-freedom.
            For more information, see the OpenSeesPy documentation. Defaults to "RCM".
        @param system_type (str, optional): Type of system of equations chosen. It determines how to construct the LinearSOE and LinearSolver objects to store and solve the system of equations in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "BandGeneral".
        @param analysis_type (str, optional): Type of analysis chosen. It determines how to construct the Analysis object, which defines what type of analysis is to be performed.
            For more information, see the OpenSeesPy documentation. Defaults to "Static".
        @param show_plot (bool, optional): Option to show the 'vertical displacement vs. vertical loading' curve after the analysis. Defaults to False.

        @exception WrongDimension: The dimension of the loaded_nodes and Fy arguments needs to be the same.
        @exception NegativeValue: The ID of timeSeries_ID needs to be a positive integer.
        @exception NegativeValue: The ID of pattern_ID needs to be a positive integer.
        """
        if len(loaded_nodes) != len(Fy): raise WrongDimension()
        if timeSeries_ID < 1: raise NegativeValue()
        if pattern_ID < 1: raise NegativeValue()

        # for mass defined: opsplt.createODB(self.name_ODB, "Gravity", Nmodes = nEigen); 
        # for tracking gravity with ODB: opsplt.createODB(self.name_ODB, "Gravity");

        # Create load pattern
        timeSeries(timeSeries_type, timeSeries_ID)
        pattern(pattern_type, timeSeries_ID, pattern_ID)
        for ii, node_ID in enumerate(loaded_nodes):
            load(node_ID, 0.0, Fy[ii], 0.0)     # load(IDNode, Fx, Fy, Mz)
        DGravity = 1.0/n_step                   # load increment                  

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
        print("")
        print("Gravity analysis starts")
        for iteration in range(n_step):
            convergence = self.__LoadCtrlLoop(DGravity, iteration,
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

        print("")
        print("Gravity complete")


    def LateralForce(self, loaded_nodes: list, Fx: list, timeSeries_ID: int, pattern_ID: int, n_step = 1000, fiber_ID_analysed = -1, fiber_section = 1,
        timeSeries_type = "Linear", pattern_type = "Plain", constraints_type = "Plain", numberer_type = "RCM", system_type = "BandGeneral", analysis_type = "Static",
        show_plot = True, block = False):
        """
        Method to perform the lateral force analyisis with lateral loading (load-control).
        If this method is called, the LoadingProtocol and Pushover methods should be avoided.

        @param loaded_nodes (list): List of nodes that are loaded by the the forces in Fx. The first node will be recorded (thus usually should be in the roof).
        @param Fx (list): List of horizontal loadings (negative is toward left; see global coordinate system).
        @param timeSeries_ID (int): ID of the timeseries.
        @param pattern_ID (int): ID of the pattern.
        @param n_step (int, optional): Number of steps used to during the analysis to reach the objective state (with 100% horizontal loadings imposed). Defaults to 1000.
        @param fiber_ID_analysed (int, optional): The ID of the analysed fiber. If fibers are present in the model and the user wants to save ODB data
            (to use in the post-processing with for example FiberResponse), assign to this argument the ID of the fiber chosen.
            -1 will ignore the storage of data for fibers. Defaults to -1.
        @param fiber_section (int, optional): The section number, i.e. the Gauss integratio number.
            If the fiber_ID_analysed is equal to -1, this argument is not used. Defaults to 1.
        @param timeSeries_type (str, optional): Type of timeseries chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Linear".
        @param pattern_type (str, optional): Type of pattern chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param constraints_type (str, optional): Type of contraints chosen. It detemines how the constraint equations are enforced in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param numberer_type (str, optional): Type of numberer chosen. It determines the mapping between equation numbers and degrees-of-freedom.
            For more information, see the OpenSeesPy documentation. Defaults to "RCM".
        @param system_type (str, optional): Type of system of equations chosen. It determines how to construct the LinearSOE and LinearSolver objects to store and solve the system of equations in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "BandGeneral".
        @param analysis_type (str, optional): Type of analysis chosen. It determines how to construct the Analysis object, which defines what type of analysis is to be performed.
            For more information, see the OpenSeesPy documentation. Defaults to "Static".
        @param show_plot (bool, optional): Option to show the 'Horizontal displacement vs. Horizontal loading' curve after the analysis. Defaults to True.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.

        @exception WrongDimension: The dimension of the loaded_nodes and Fx arguments needs to be the same.
        @exception NegativeValue: The ID of timeSeries_ID needs to be a positive integer.
        @exception NegativeValue: The ID of pattern_ID needs to be a positive integer.
        @exception NegativeValue: The ID of fiber_ID_analysed needs to be a positive integer.
        """
        if len(loaded_nodes) != len(Fx): raise WrongDimension()
        if timeSeries_ID < 1: raise NegativeValue()
        if pattern_ID < 1: raise NegativeValue()
        if fiber_ID_analysed != -1 and fiber_ID_analysed < 1: raise NegativeValue()

        # for mass defined: opsplt.createODB(self.name_ODB, "LateralForce", Nmodes = nEigen); 
        opsplt.createODB(self.name_ODB, "LateralForce");
        if fiber_ID_analysed != -1: opsplt.saveFiberData2D(self.name_ODB, "LateralForce", fiber_ID_analysed, fiber_section)

        # Create load pattern
        timeSeries(timeSeries_type, timeSeries_ID)
        pattern(pattern_type, timeSeries_ID, pattern_ID)
        for ii, node_ID in enumerate(loaded_nodes):
            load(node_ID, Fx[ii], 0.0, 0.0)     # load(IDNode, Fx, Fy, Mz)
        force = 1.0/n_step                   # load increment                  

        # Set up analysis options
        constraints(constraints_type)   # how it handles boundary conditions
        numberer(numberer_type)         # renumber dof's to minimize band-width (optimization)
        system(system_type)             # how to store and solve the system of equations in the analysis
                                        # For static model, BandGeneral, for transient and/or big model, UmfPack
        integrator("LoadControl", force)# LoadControl and DisplacementControl only with static model, linear TimeSeries w/ factor of 1
                                        # Newmark used for transient model  
        algorithm("Newton")             # placeholder
        analysis(analysis_type)         # define type of analysis: static for pushover

        # Analysis
        dataLF = np.zeros((n_step+1,2))
        print("")
        print("Lateral Force analysis starts")
        for iteration in range(n_step):
            convergence = self.__LoadCtrlLoop(force, iteration,
                self.algo, self.test_type, self.tol, self.test_opt, self.max_iter)
            if convergence != 0: break
            dataLF[iteration+1,0] = nodeDisp(loaded_nodes[0], 1)/mm_unit
            dataLF[iteration+1,1] = getLoadFactor(pattern_ID)*Fx[0]/kN_unit

        if show_plot:
            plt.plot(dataLF[:,0], dataLF[:,1])
            plt.xlabel('Lateral Displacement [mm]')
            plt.ylabel('Lateral Load [kN]')
            plt.title('Lateral force curve')
            if block:
                plt.show()

        loadConst("-time", 0.0)

        print("")
        print("Lateral force complete")
        self.load_case = "LateralForce"

        wipe()


    def Pushover(self, CtrlNode: int, Dmax, Dincr, timeSeries_ID: int, pattern_ID: int, Fx = 1*kN_unit, ele_fiber_ID_analysed = -1, fiber_section = 1,
        timeSeries_type = "Linear", pattern_type = "Plain", constraints_type = "Plain", numberer_type = "RCM", system_type = "UmfPack", analysis_type = "Static",
        show_plot = True, block = False):
        """
        Method to perform a pushover analysis (displacement-control). If this method is called, the LoadingProtocol and LateralForce methods should be avoided.

        @param CtrlNode (int): The node that will be used to impose the displacement Dmax of the pushover analysis.
            If the show_plot option is True, the curve displayed follows this node.
        @param Dmax (float): The imposed displacement.
        @param Dincr (float): The incremental displacement to reach Dmax. To converge, it should be small enough (1000 times smaller of Dmax).
        @param timeSeries_ID (int): ID of the timeseries.
        @param pattern_ID (int): ID of the pattern.
        @param Fx (float, optional): The force imposed at the control node CtrlNode. It is used for convergence reasons and it can be arbitrarly small.
            Defaults to 1*kN_unit.
        @param ele_fiber_ID_analysed (int, optional): The ID of the analysed element with fibers. If fibers are present in the model and the user wants to save ODB data
            (to use in the post-processing with for example FiberResponse), assign to this argument the ID of the element with fibers chosen.
            -1 will ignore the storage of data for fibers. Defaults to -1.
        @param fiber_section (int, optional): The section number, i.e. the Gauss integratio number.
            If the fiber_ID_analysed is equal to -1, this argument is not used. Defaults to 1.
        @param timeSeries_type (str, optional): Type of timeseries chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Linear".
        @param pattern_type (str, optional): Type of pattern chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param constraints_type (str, optional): Type of contraints chosen. It detemines how the constraint equations are enforced in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param numberer_type (str, optional): Type of numberer chosen. It determines the mapping between equation numbers and degrees-of-freedom.
            For more information, see the OpenSeesPy documentation. Defaults to "RCM".
        @param system_type (str, optional): Type of system of equations chosen. It determines how to construct the LinearSOE and LinearSolver objects to store and solve the system of equations in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "UmfPack".
        @param analysis_type (str, optional): Type of analysis chosen. It determines how to construct the Analysis object, which defines what type of analysis is to be performed.
            For more information, see the OpenSeesPy documentation. Defaults to "Static".
        @param show_plot (bool, optional): Option to show the 'lateral displacement vs. lateral loading' curve after the analysis. Defaults to True.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.

        @exception NegativeValue: The ID of CtrlNode needs to be a positive integer.
        @exception NegativeValue: The ID of timeSeries_ID needs to be a positive integer.
        @exception NegativeValue: The ID of pattern_ID needs to be a positive integer.
        @exception NegativeValue: The ID of ele_fiber_ID_analysed needs to be a positive integer if is different from -1.
        """
        if CtrlNode < 1: raise NegativeValue()
        if timeSeries_ID < 1: raise NegativeValue()
        if pattern_ID < 1: raise NegativeValue()
        if ele_fiber_ID_analysed != -1 and ele_fiber_ID_analysed < 1: raise NegativeValue()

        # for mass defined: opsplt.createODB(self.name_ODB, "Pushover", Nmodes = nEigen); 
        opsplt.createODB(self.name_ODB, "Pushover");
        if ele_fiber_ID_analysed != -1: opsplt.saveFiberData2D(self.name_ODB, "Pushover", ele_fiber_ID_analysed, fiber_section)

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
        integrator("LoadControl", 1)    # placeholder
        algorithm("Newton")             # placeholder
        analysis(analysis_type)         # define type of analysis: static for pushover

        # Analysis
        dataPO = np.zeros((Nsteps+1,2))
        next_step = 0
        print("")
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
            if block:
                plt.show()
        
        print("")
        print("Pushover complete")
        self.load_case = "Pushover"

        wipe()


    def LoadingProtocol(self, CtrlNode: int, discr_LP: np.ndarray, timeSeries_ID: int, pattern_ID: int, Fx = 1*kN_unit, ele_fiber_ID_analysed = -1, fiber_section = 1,
        timeSeries_type = "Linear", pattern_type = "Plain", constraints_type = "Plain", numberer_type = "RCM", system_type = "UmfPack", analysis_type = "Static",
        show_plot = True, block = False):
        """
        Method to perform a loading protocol analysis (displacement-control). If this method is called, the Pushover and LateralForce methods should be avoided.

        @param CtrlNode (int): The node that will be used to impose the displacement from the discr_LP to perform the analysis.
        @param discr_LP (np.ndarray): The loading protocol array (1 dimension) discretised. It needs to be filled with imposed displacement, not SDR.
            Use the functions DiscretizeLoadProtocol and DiscretizeLinearly in FunctionalFeatures module to help create and/or discretise one.
        @param timeSeries_ID (int): ID of the timeseries.
        @param pattern_ID (int): ID of the pattern.
        @param Fx (float, optional): The force imposed at the control node CtrlNode. It is used for convergence reasons and it can be arbitrarly small.
            Defaults to 1*kN_unit.
        @param ele_fiber_ID_analysed (int, optional): The ID of the analysed element with fibers. If fibers are present in the model and the user wants to save ODB data
            (to use in the post-processing with for example FiberResponse), assign to this argument the ID of the element with fibers chosen.
            -1 will ignore the storage of data for fibers. Defaults to -1.
        @param fiber_section (int, optional): The section number, i.e. the Gauss integratio number.
            If the fiber_ID_analysed is equal to -1, this argument is not used. Defaults to 1.
        @param timeSeries_type (str, optional): Type of timeseries chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Linear".
        @param pattern_type (str, optional): Type of pattern chosen.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param constraints_type (str, optional): Type of contraints chosen. It detemines how the constraint equations are enforced in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "Plain".
        @param numberer_type (str, optional): Type of numberer chosen. It determines the mapping between equation numbers and degrees-of-freedom.
            For more information, see the OpenSeesPy documentation. Defaults to "RCM".
        @param system_type (str, optional): Type of system of equations chosen. It determines how to construct the LinearSOE and LinearSolver objects to store and solve the system of equations in the analysis.
            For more information, see the OpenSeesPy documentation. Defaults to "UmfPack".
        @param analysis_type (str, optional): Type of analysis chosen. It determines how to construct the Analysis object, which defines what type of analysis is to be performed.
            For more information, see the OpenSeesPy documentation. Defaults to "Static".
        @param show_plot (bool, optional): Option to show the 'lateral displacement vs. lateral loading' curve after the analysis. Defaults to True.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.

        @exception NegativeValue: The ID of CtrlNode needs to be a positive integer.
        @exception NegativeValue: The ID of timeSeries_ID needs to be a positive integer.
        @exception NegativeValue: The ID of pattern_ID needs to be a positive integer.
        @exception NegativeValue: The ID of fiber_ID_analysed needs to be a positive integer if is different from -1.
        """
        if CtrlNode < 1: raise NegativeValue()
        if timeSeries_ID < 1: raise NegativeValue()
        if pattern_ID < 1: raise NegativeValue()
        if ele_fiber_ID_analysed != -1 and ele_fiber_ID_analysed < 1: raise NegativeValue()

        # for mass defined: opsplt.createODB(self.name_ODB, "LoadingProtocol", Nmodes = nEigen); 
        opsplt.createODB(self.name_ODB, "LoadingProtocol");
        if ele_fiber_ID_analysed != -1: opsplt.saveFiberData2D(self.name_ODB, "LoadingProtocol", ele_fiber_ID_analysed, fiber_section)

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
        integrator("LoadControl", 1)    # placeholder
        algorithm("Newton")             # placeholder
        analysis(analysis_type)         # define type of analysis: static for LoadingProtocol

        # Analysis
        dataLP = np.zeros((Nsteps+1,2))
        next_step = 0
        print("")
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
            if block:
                plt.show()
        
        print("")
        print("Loading Protocol complete")
        self.load_case = "LoadingProtocol"
        
        wipe()


    def __LoadCtrlLoop(self, force, iteration: int, algo = "KrylovNewton", test_type = "NormDispIncr", tol = TOL, test_opt = 0, max_iter = MAX_ITER):
        """
        PRIVATE METHOD. It is used perform one load increment 'force' load-control analysis step using 'algo' and 'test_type' as algorithm and test.
        The integrator is LoadControl. If convergence issues are encountered, the method performa a convergence analysis trying different ways to converge.

        @param force (dougle): The load increment performed.
        @param iteration (int): The current iteration.
        @param algo (str, optional): Type of alghoritm chosen for the analysis. It detemines how to construct a SolutionAlgorithm object, which determines the sequence of steps taken to solve the non-linear equation.
            For more information on the available types, see the OpenSeesPy documentation. Defaults to "KrylovNewton".
        @param test_type (str, optional): Type of test chosen for the analysis. It determines how to construct a ConvergenceTest object.
            Certain SolutionAlgorithm objects require a ConvergenceTest object to determine if convergence has been achieved at the end of an iteration step.
            For more information on the available types, see the OpenSeesPy documentation. Defaults to "NormDispIncr".
        @param tol (float, optional): Tolerance criteria used to check for convergence. Defaults to TOL (from Constants Module).
        @param test_opt (int, optional): Print-flag from 0 to 5 used to receive more info during the iteration
            (for example: 0 print nothing and 2 print information on norms and number of iterations at end of successful test).
            For more information, see the OpenSeesPy documentation. Defaults to 0.
        @param max_iter (float, optional): Maximal number of iterations to check. Defaults to MAX_ITER (from Constants Module).

        @exception NegativeValue: iteration needs to be a positive integer.
        @exception NegativeValue: tol needs to be positive.
        @exception NegativeValue: max_iter needs to be positive.
        
        @returns int: 0 if the interation converged. 
        """
        if iteration < 0: raise NegativeValue()
        if tol < 0: raise NegativeValue()
        if max_iter < 0: raise NegativeValue()

        # Default analysis
        integrator("LoadControl", force)            # LoadControl and DisplacementControl only with static model, linear TimeSeries w/ factor of 1
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
                print("NO CONVERGENCE! Load-control analysis stops at iteration {} ".format(iteration))
                print("#############################################################################")
                print("")
            else:
                print("Convergence reached, convergence analysis ends.")

        return convergence


    def __LatDispCtrlLoop(self, CtrlNode, dU, iteration, algo = "KrylovNewton", test_type = "NormDispIncr", tol = TOL, test_opt = 0, max_iter = MAX_ITER, allow_smaller_step = False):
        """
        PRIVATE METHOD. It is used perform one imposed displacement increment 'dU' displacement-control analysis step using 'algo' and 'test_type' as algorithm and test.
        The integrator is DisplacementControl. If convergence issues are encountered, the method performa a convergence analysis trying different ways to converge.

        @param CtrlNode ([type]): [description]
        @param dU (float): The imposed displacement increment for the current iteration
        @param iteration (int): The current iteration.
        @param algo (str, optional): Type of alghoritm chosen for the analysis. It detemines how to construct a SolutionAlgorithm object, which determines the sequence of steps taken to solve the non-linear equation.
            For more information on the available types, see the OpenSeesPy documentation. Defaults to "KrylovNewton".
        @param test_type (str, optional): Type of test chosen for the analysis. It determines how to construct a ConvergenceTest object.
            Certain SolutionAlgorithm objects require a ConvergenceTest object to determine if convergence has been achieved at the end of an iteration step.
            For more information on the available types, see the OpenSeesPy documentation. Defaults to "NormDispIncr".
        @param tol (float, optional): Tolerance criteria used to check for convergence. Defaults to TOL (from Constants Module).
        @param test_opt (int, optional): Print-flag from 0 to 5 used to receive more info during the iteration
            (for example: 0 print nothing and 2 print information on norms and number of iterations at end of successful test).
            For more information, see the OpenSeesPy documentation. Defaults to 0.
        @param max_iter (float, optional): Maximal number of iterations to check. Defaults to MAX_ITER (from Constants Module).
        @param allow_smaller_step (bool, optional): Allow smaller steps in the displacement-control analysis. Defaults to False.

        @exception NegativeValue: iteration needs to be a positive integer.
        @exception NegativeValue: tol needs to be positive.
        @exception NegativeValue: max_iter needs to be positive.
        
        @returns int: 0 if the interation converged. 
        """
        if iteration < 0: raise NegativeValue()
        if tol < 0: raise NegativeValue()
        if max_iter < 0: raise NegativeValue()

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


    def DeformedShape(self, scale = 1, animate = False, dt = 0.01):
        """
        Method that shows the final deformed shape of the model. It can also show the animation that shows how the model behaved during the analysis.

        @param scale (int, optional): The scaling factor to magnify the deformation. The value should be adjusted for each model. Defaults to 1.
        @param animate (bool, optional): Option to show the animation of the model during the analysis. Defaults to False.
        @param dt (float, optional): The time step between every iteration. Defaults to 0.01.
        
        @exception NameError: The methods for the analysis were not called.
        """
        if self.load_case == "None": raise NameError("The analysis is not complete.")
        
        # Display deformed shape, the scaling factor needs to be adjusted for each model
        opsplt.plot_deformedshape(Model = self.name_ODB, LoadCase=self.load_case, scale = scale)
        if animate:
            opsplt.animate_deformedshape(Model = self.name_ODB, LoadCase=self.load_case, dt = dt, scale = scale)


    def FiberResponse(self, ele_fiber_ID_analysed: int, fiber_section = 1, animate_stress = False, animate_strain = False, fps = 25):
        """
        Method that shows the final stress response of the fiber section chosen.
        It can also show the animation that shows how the fiber section behaved during the analysis. The fiber ID and section needs to be recorded during the analysis,
        thus if the method LateralForce, Pushover or LoadingProtocol was used, the same fiber ID and section need to be used. 

        @param ele_fiber_ID_analysed (int): The ID of the analysed fiber. If fibers are present in the model and the user wants to save ODB data
            (to use in the post-processing with for example FiberResponse), assign to this argument the ID of the fiber chosen.
            -1 will ignore the storage of data for fibers.
        @param fiber_section (int, optional): The section number, i.e. the Gauss integratio number.
            If the fiber_ID_analysed is equal to -1, this argument is not used. Defaults to 1.
        @param animate_stress (bool, optional): Option to show the animation of the fiber stress during the analysis. Defaults to False.
        @param animate_strain (bool, optional): Option to show the animation of the fiber strain during the analysis. Defaults to False.
        @param fps (int, optional): Number of frame per seconds for the animations. Defaults to 25.
        @exception NameError: The methods for the analysis were not called.
        """
        if self.load_case == "None": raise NameError("The analysis is not complete.")

        opsplt.plot_fiberResponse2D(self.name_ODB, self.load_case, ele_fiber_ID_analysed, fiber_section, InputType = 'stress')
        if animate_stress:
            ani1 = opsplt.animate_fiberResponse2D(self.name_ODB, self.load_case, ele_fiber_ID_analysed, fiber_section, InputType = 'stress', fps = fps)
        if animate_strain:
            ani1 = opsplt.animate_fiberResponse2D(self.name_ODB, self.load_case, ele_fiber_ID_analysed, fiber_section, InputType = 'strain', fps = fps)


def _ConvergenceTest(algo_type: str, test_type: str, tol, max_iter, test_opt = 0):
    """
    PRIVATE FUNCTION. It is used during the convergence analysis to test different ways to reach convergence.

    @param algo_type (str): Type of alghoritm chosen for the analysis. It detemines how to construct a SolutionAlgorithm object, which determines the sequence of steps taken to solve the non-linear equation.
        For more information on the available types, see the OpenSeesPy documentation.
    @param test_type (str): Type of test chosen for the analysis. It determines how to construct a ConvergenceTest object.
        Certain SolutionAlgorithm objects require a ConvergenceTest object to determine if convergence has been achieved at the end of an iteration step.
        For more information on the available types, see the OpenSeesPy documentation.
    @param tol (float): Tolerance criteria used to check for convergence.
    @param max_iter (float): Maximal number of iterations to check.
    @param test_opt (int, optional): Print-flag from 0 to 5 used to receive more info during the iteration
        (for example: 0 print nothing and 2 print information on norms and number of iterations at end of successful test).
        For more information, see the OpenSeesPy documentation. Defaults to 0.

    @returns int: 0 if the interation converged.
    """
    print("algorithm: {}".format(algo_type))
    print("test: {}, tol = {}, max iter = {}".format(test_type, tol, max_iter))
    test(test_type, tol, max_iter, test_opt)
    algorithm(algo_type)

    return analyze(1)



