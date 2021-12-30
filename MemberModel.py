# Module with the member model
#   Carmine Schipani, 2021

from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from abc import abstractmethod
from copy import copy, deepcopy
import openseespy.postprocessing.internal_plotting_functions as ipltf
import openseespy.postprocessing.Get_Rendering as opsplt
from OpenSeesPyAssistant.Section import *
from OpenSeesPyAssistant.DataManagement import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *
from OpenSeesPyAssistant.Constants import *
from OpenSeesPyAssistant.Fibers import *
from OpenSeesPyAssistant.Connections import *
from OpenSeesPyAssistant.FunctionalFeatures import *

# Member model
class MemberModel(DataManagement):
    pass

class PanelZone(MemberModel):
    # Class that stores funcions and material properties of a panel zone.
    # Warning: the units should be m and N
    def __init__(self, master_node_ID: int, mid_panel_zone_width, mid_panel_zone_height, E, A_rigid, I_rigid, geo_transf_ID: int, mat_ID):

        # Check
        if master_node_ID < 1: raise NegativeValue()
        if master_node_ID > 99: raise WrongNodeIDConvention(master_node_ID)
        if mid_panel_zone_width < 0: raise NegativeValue()
        if mid_panel_zone_height < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if A_rigid < 0: raise NegativeValue()
        if I_rigid < 0: raise NegativeValue()
        if geo_transf_ID > 1: raise NegativeValue()
        if mat_ID < 0: raise NegativeValue()

        # Arguments
        self.master_node_ID = master_node_ID
        self.mid_panel_zone_width = mid_panel_zone_width
        self.mid_panel_zone_height = mid_panel_zone_height
        self.E = E
        self.A_rigid = A_rigid
        self.I_rigid = I_rigid
        self.geo_transf_ID = geo_transf_ID
        self.mat_ID = mat_ID 

        # Initialized the parameters that are dependent from others
        self.col_section_name_tag = "None"
        self.beam_section_name_tag = "None"
        self.Initialized = False
        self.ReInit()

    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Arguments

        # Members
        if self.col_section_name_tag != "None": self.col_section_name_tag = self.col_section_name_tag + " (modified)"
        if self.beam_section_name_tag != "None": self.beam_section_name_tag = self.beam_section_name_tag + " (modified)"

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "PanelZone"], # Tag for differentiating different data
            ["master_node_ID", self.master_node_ID],
            ["col_section_name_tag", self.col_section_name_tag],
            ["beam_section_name_tag", self.beam_section_name_tag],
            ["mat_ID", self.mat_ID],
            ["mid_panel_zone_width", self.mid_panel_zone_width],
            ["mid_panel_zone_height", self.mid_panel_zone_height],
            ["E", self.E],
            ["A_rigid", self.A_rigid],
            ["I_rigid", self.I_rigid],
            ["tranf_ID", self.geo_transf_ID],
            ["Initialized", self.Initialized]]

    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the member model (optional).
        """
        print("")
        print("Requested info for Panel Zone member model, master node ID = {}".format(self.master_node_ID))
        print("Section associated, column: {} ".format(self.col_section_name_tag))
        print("Section associated, beam: {} ".format(self.beam_section_name_tag))
        print("Material model of the panel zone ID = {}".format(self.mat_ID))
        print("Mid panel zone width = {} mm".format(self.mid_panel_zone_width/mm_unit))
        print("Mid panel zone height = {} mm".format(self.mid_panel_zone_height/mm_unit))
        print("Young modulus E = {} GPa".format(self.E/GPa_unit))
        print("Area of the elements (rigid) = {} mm2".format(self.A_rigid/mm2_unit))
        print("Moment of inetia of the elements (strong axis, rigid) = {} mm4".format(self.I_rigid/mm4_unit))
        print("Geometric transformation = {}".format(self.geo_transf_ID))
        print("")

        if plot:
            if self.Initialized:
                plot_member(np.array(self.element_array))
                if block:
                    plt.show()
            else:
                print("The panel zone is not initialized (node and elements not created) for master node ID = {}".format(self.master_node_ID))


    def CreateMember(self):
        # Define nodes
        DefinePanelZoneNodes(self.master_node_ID, self.mid_panel_zone_width, self.mid_panel_zone_height)
        xy1 = IDConvention(self.master_node_ID, 1)
        xy01 = IDConvention(self.master_node_ID, 1, 1)
        xy03 = IDConvention(self.master_node_ID, 3, 1)
        xy04 = IDConvention(self.master_node_ID, 4, 1)
        xy06 = IDConvention(self.master_node_ID, 6, 1)
        xy07 = IDConvention(self.master_node_ID, 7, 1)
        xy09 = IDConvention(self.master_node_ID, 9, 1)
        xy10 = IDConvention(self.master_node_ID, 10)

        # Define rigid elements
        self.element_array = DefinePanelZoneElements(self.master_node_ID, self.E, self.A_rigid, self.I_rigid, self.geo_transf_ID)
        
        # Define zero length element
        spring_ID = IDConvention(xy1, xy01)
        RotationalSpring(spring_ID, xy1, xy01, self.mat_ID)
        self.element_array.append([spring_ID, xy1, xy01])

        # Pin connections
        Pin(xy03, xy04)
        Pin(xy06, xy07)
        Pin(xy09, xy10)

        # Update class
        self.Initialized = True
        self.UpdateStoredData()


class PanelZoneSteelIShape(PanelZone):
    def __init__(self, master_node_ID: int, col: SteelIShape, beam: SteelIShape, geo_transf_ID: int, mat_ID, rigid = RIGID):
        self.col = col
        self.beam = beam
        super().__init__(master_node_ID, col.d/2.0, beam.d/2.0, col.E, max(col.A, beam.A)*rigid, max(col.Iy, beam.Iy)*rigid, geo_transf_ID, mat_ID)

        self.col_section_name_tag = col.name_tag
        self.beam_section_name_tag = beam.name_tag
        self.UpdateStoredData()

class PanelZoneSteelIShapeGupta1999(PanelZoneSteelIShape):
    def __init__(self, master_node_ID: int, col: SteelIShape, beam: SteelIShape, geo_transf_ID: int, t_dp = 0, rigid=RIGID):
        self.col = col
        self.beam = beam
        mat_ID = master_node_ID
        pz_spring = Gupta1999SteelIShape(mat_ID, col, beam, t_dp)
        pz_spring.Hysteretic()

        super().__init__(master_node_ID, col, beam, geo_transf_ID, mat_ID, rigid)

class PanelZoneSteelIShapeSkiadopoulos2021(PanelZoneSteelIShape):
    def __init__(self, master_node_ID: int, col: SteelIShape, beam: SteelIShape, geo_transf_ID: int, t_dp = 0, rigid=RIGID):
        self.col = col
        self.beam = beam
        mat_ID = master_node_ID
        pz_spring = Skiadopoulos2021SteelIShape(mat_ID, col, beam, t_dp)
        pz_spring.Hysteretic()

        super().__init__(master_node_ID, col, beam, geo_transf_ID, mat_ID, rigid)





def DefinePanelZoneNodes(MasterNodeID, MidPanelZoneWidth, MidPanelZoneHeight):
	######################################################################################
	## DefinePanelZoneNodes
	######################################################################################
	## Define the remaining 10 nodes of a panel zone given the dimensions and the master node (top center one).
	## Note that the top right node is defined differently because is where the spring is
	##			Carmine Schipani, 2021
	##
	## 	MasterNodeID : 			The top center node. The conventional denomination is xy, x = Pier axis, y = Floor (ground = 1)
	##	MidPanelZoneWidth :		Half the panel zone width (that should be equal to half the column depth)
	##	MidPanelZoneHeight :	Half the panel zone height (that should be equal to half the beam depth)
	##	AxisCL :				The x coordinate of the centerline of the column
	##	FloorCL :				The y coordinate of the centerline of the beam
	##	
	#######################################################################################
    # 		PZNodeID:		12 nodes: top right xy (master), xy1 top right,						xy09, xy10 	   xy 		xy1, xy01					
    # 						clockwise 10 nodes xy01-xy10 (with double node at corners)				o-----------o-----------o			
    #						Spring at node xy1														|						|
    #		PZElemeneID:	8 elements: starting at node xy, clockwise								|						|
    #						(see function DefinePanelZoneElements for more info)					|						|
    #																								|						|
    #																						   xy08	o						o xy02
    #																								|						|
    #																								|						|
    #																								|						|
    #																								|						|
    #																								o-----------o-----------o
    #																							xy06, xy07 	   xy05  	xy03, xy04

    # Get node coord and define useful variables
    m_node = np.array(nodeCoord(MasterNodeID))
    AxisCL = m_node[0]
    FloorCL = m_node[1] - MidPanelZoneHeight

	# Convention: Node of the spring (top right) is xy1
    node(IDConvention(MasterNodeID, 1), AxisCL+MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)
	# Convention: Two notes in the corners (already defined one, xy1) clockwise from xy01 to xy10
    node(IDConvention(MasterNodeID, 1, 1), AxisCL+MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)
    node(IDConvention(MasterNodeID, 2, 1), AxisCL+MidPanelZoneWidth, FloorCL)
    node(IDConvention(MasterNodeID, 3, 1), AxisCL+MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(IDConvention(MasterNodeID, 4, 1), AxisCL+MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(IDConvention(MasterNodeID, 5, 1), AxisCL, FloorCL-MidPanelZoneHeight)
    node(IDConvention(MasterNodeID, 6, 1), AxisCL-MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(IDConvention(MasterNodeID, 7, 1), AxisCL-MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(IDConvention(MasterNodeID, 8, 1), AxisCL-MidPanelZoneWidth, FloorCL)
    node(IDConvention(MasterNodeID, 9, 1), AxisCL-MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)
    node(IDConvention(MasterNodeID, 10), AxisCL-MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)


def DefinePanelZoneElements(MasterNodeID, E, RigidA, RigidI, TransfID):
    ######################################################################################
    ## DefinePanelZoneElements
    ######################################################################################
    ## Defines the 8 panel zone elements. Warning: the algorithm for the elements ID work only for MasterNodeID of 2 digit.
    ##			Carmine Schipani, 2021
    ##
    ## 	MasterNodeID : 	The top center node. The conventional denomination is xy, x = Pier axis, y = Floor (ground = 1)
    ##	E :				Young's modulus
    ##	RigidA :		Area two ordrs bigger that the members
    ##	RigidI:			Moment of inertia two ordrs bigger that the members
    ##	TransfID :		Geometric transformation ID
    # 		PZNodeID:		12 nodes: top right xy (master), xy1 top right,						xy09, xy10 	   xy 		xy1, xy01					
    # 						clockwise 10 nodes xy01-xy10 (with double node at corners)				o-----------o-----------o			
    #						Spring at node xy1														|						|
    #		PZElemeneID:	8 elements: starting at node xy, clockwise								|						|
    #						(see function DefinePanelZoneElements for more info)					|						|
    #																								|						|
    #																						   xy08	o						o xy02
    #																								|						|
    #																								|						|
    #																								|						|
    #																								|						|
    #																								o-----------o-----------o
    #																							xy06, xy07 	   xy05  	xy03, xy04


    # Compute the ID of the nodes obeying to the convention used
    xy = MasterNodeID
    xy1 = IDConvention(xy, 1)
    xy01 = IDConvention(xy, 1, 1)
    xy02 = IDConvention(xy, 2, 1)
    xy03 = IDConvention(xy, 3, 1)
    xy04 = IDConvention(xy, 4, 1)
    xy05 = IDConvention(xy, 5, 1)
    xy06 = IDConvention(xy, 6, 1)
    xy07 = IDConvention(xy, 7, 1)
    xy08 = IDConvention(xy, 8, 1)
    xy09 = IDConvention(xy, 9, 1)
    xy10 = IDConvention(xy, 10)

    # Create element IDs using the convention: xy(a)xy(a)	with xy(a) = NodeID i and j
    #	Starting at MasterNodeID, clockwise
    if MasterNodeID > 99:
        print("Warning, convention: MasterNodeID's digits should be 2")

    ele1 = IDConvention(xy, xy1)
    ele2 = IDConvention(xy01, xy02)
    ele3 = IDConvention(xy02, xy03)
    ele4 = IDConvention(xy04, xy05)
    ele5 = IDConvention(xy05, xy06)
    ele6 = IDConvention(xy07, xy08)
    ele7 = IDConvention(xy08, xy09)
    ele8 = IDConvention(xy10, xy)

    # Create panel zone elements
    #                             ID   ndI   ndJ      A    E     I     Transf
    element("elasticBeamColumn", ele1, xy, 	 xy1,  RigidA, E, RigidI, TransfID)
    element("elasticBeamColumn", ele2, xy01, xy02, RigidA, E, RigidI, TransfID)
    element("elasticBeamColumn", ele3, xy02, xy03, RigidA, E, RigidI, TransfID)
    element("elasticBeamColumn", ele4, xy04, xy05, RigidA, E, RigidI, TransfID)
    element("elasticBeamColumn", ele5, xy05, xy06, RigidA, E, RigidI, TransfID)
    element("elasticBeamColumn", ele6, xy07, xy08, RigidA, E, RigidI, TransfID)
    element("elasticBeamColumn", ele7, xy08, xy09, RigidA, E, RigidI, TransfID)
    element("elasticBeamColumn", ele8, xy10, xy, RigidA, E, RigidI, TransfID)

    # Create element array for forther manipulations
    element_array = [[ele1, xy, xy1],
        [ele2, xy01, xy02],
        [ele3, xy02, xy03],
        [ele4, xy04, xy05],
        [ele5, xy05, xy06],
        [ele6, xy07, xy08],
        [ele7, xy08, xy09],
        [ele8, xy10, xy]]
    
    return element_array


class ElasticElement(MemberModel):
    # Class that stores funcions and material properties of a elastic element.
    # Warning: the units should be m and N
    # Convention:																							
    # 		NodeID: 		xy 			with x = pier, y = floor 											o xy7	|	
    #		PlHingeID:		xya			with x = pier, y = floor, a:	 --o xy  xy2 o-----o xy3  xy o--	|		o xy
    #																										|	
    #																										|		o xy
    #																										o xy6	|
    #		ElementID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		TrussID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		PDeltaColID:	xy(a)xy(a)	with xy(a) = NodeID i and j
    #		Spring:			xy(a)xy(a)	with xy(a) = NodeID i and j
    
    def __init__(self, iNode_ID: int, jNode_ID: int, A, E, Iy, geo_transf_ID: int):

        # Check
        if iNode_ID < 1: raise NegativeValue()
        if jNode_ID < 1: raise NegativeValue()
        if A < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if Iy < 0: raise NegativeValue()
        if geo_transf_ID < 1: raise NegativeValue()

        # Arguments
        self.iNode_ID = iNode_ID
        self.jNode_ID = jNode_ID
        self.A = A
        self.E = E
        self.Iy = Iy
        self.geo_transf_ID = geo_transf_ID

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.Initialized = False
        self.ReInit()

    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # element ID
        self.element_ID = IDConvention(self.iNode_ID, self.jNode_ID)

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "ElasticElement"], # Tag for differentiating different data
            ["element_ID", self.element_ID],
            ["section_name_tag", self.section_name_tag],
            ["A", self.A],
            ["E", self.E],
            ["Iy", self.Iy],
            ["iNode_ID", self.iNode_ID],
            ["jNode_ID", self.jNode_ID],
            ["tranf_ID", self.geo_transf_ID],
            ["Initialized", self.Initialized]]

    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the member model (optional).
        """
        print("")
        print("Requested info for ElasticElement member model, ID = {}".format(self.element_ID))
        print("Section associated {} ".format(self.section_name_tag))
        print("Area A = {} mm2".format(self.A/mm2_unit))
        print("Young modulus E = {} GPa".format(self.E/GPa_unit))
        print("Moment of inertia Iy = {} mm4".format(self.Iy/mm4_unit))
        print("Geometric transformation = {}".format(self.geo_transf_ID))
        print("")

        if plot:
            if self.Initialized:
                plot_member(np.array(self.element_array))
                if block:
                    plt.show()
            else:
                print("The ElasticElement is not initialized (node and elements not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        self.element_array = [[self.element_ID, self.iNode_ID, self.jNode_ID]]
        
        # Define element
        element("elasticBeamColumn", self.element_ID, self.iNode_ID, self.jNode_ID, self.A, self.E, self.Iy, self.geo_transf_ID)

        # Update class
        self.Initialized = True
        self.UpdateStoredData()


class ElasticElementSteelIShape(ElasticElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, section: SteelIShape, geo_transf_ID: int):
        super().__init__(iNode_ID, jNode_ID, section.A, section.E, section.Iy, geo_transf_ID)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()


class SpringBasedElement(MemberModel):
    # Class that stores funcions and material properties of a spring-based element. If vertical, iNode is bottom, if horizontal, iNode is left
    # Warning: the units should be m and N
    # Convention:																							
    # 		NodeID: 		xy 			with x = pier, y = floor 											o xy7	|	
    #		PlHingeID:		xya			with x = pier, y = floor, a:	 --o xy  xy2 o-----o xy3  xy o--	|		o xy
    #																										|	
    #																										|		o xy
    #																										o xy6	|
    #		ElementID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		TrussID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		PDeltaColID:	xy(a)xy(a)	with xy(a) = NodeID i and j
    #		Spring:			xy(a)xy(a)	with xy(a) = NodeID i and j
    
    def __init__(self, iNode_ID: int, jNode_ID: int, A, E, Iy_mod, geo_transf_ID: int, mat_ID_i = -1, mat_ID_j = -1):

        # Check
        if iNode_ID < 1: raise NegativeValue()
        if jNode_ID < 1: raise NegativeValue()
        if A < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if Iy_mod < 0: raise NegativeValue()
        if geo_transf_ID < 1: raise NegativeValue()
        if mat_ID_i != -1 and mat_ID_i < 0: raise NegativeValue()
        if mat_ID_j != -1 and mat_ID_j < 0: raise NegativeValue()
        if mat_ID_i == -1 and mat_ID_j == -1: raise NameError("No springs defined for element ID = {}".format(IDConvention(iNode_ID, jNode_ID)))

        # Arguments
        self.iNode_ID = iNode_ID
        self.jNode_ID = jNode_ID
        self.A = A
        self.E = E
        self.Iy_mod = Iy_mod
        self.geo_transf_ID = geo_transf_ID
        self.mat_ID_i = mat_ID_i 
        self.mat_ID_j = mat_ID_j 

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.Initialized = False
        self.ReInit()

    def ReInit(self):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Arguments

        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        # orientation:
        self.ele_orientation = NodesOrientation(self.iNode_ID, self.jNode_ID)
        if self.ele_orientation == "zero_length": raise ZeroLength(IDConvention(self.iNode_ID, self.jNode_ID))
        
        if self.ele_orientation == "vertical":
            if self.mat_ID_i != -1:
                self.iNode_ID_spring = IDConvention(self.iNode_ID, 6)
            else:
                self.iNode_ID_spring = self.iNode_ID
            
            if self.mat_ID_j != -1:
                self.jNode_ID_spring = IDConvention(self.jNode_ID, 7)
            else:
                self.jNode_ID_spring = self.jNode_ID
        else:
            if self.mat_ID_i != -1:
                self.iNode_ID_spring = IDConvention(self.iNode_ID, 2)
            else:
                self.iNode_ID_spring = self.iNode_ID
            
            if self.mat_ID_j != -1:
                self.jNode_ID_spring = IDConvention(self.jNode_ID, 3)
            else:
                self.jNode_ID_spring = self.jNode_ID
        # element ID
        self.element_ID = IDConvention(self.iNode_ID_spring, self.jNode_ID_spring)

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "SpringBasedElement"], # Tag for differentiating different data
            ["element_ID", self.element_ID],
            ["section_name_tag", self.section_name_tag],
            ["A", self.A],
            ["E", self.E],
            ["Iy_mod", self.Iy_mod],
            ["iNode_ID", self.iNode_ID],
            ["iNode_ID_spring", self.iNode_ID_spring],
            ["mat_ID_i", self.mat_ID_i],
            ["jNode_ID", self.jNode_ID],
            ["jNode_ID_spring", self.jNode_ID_spring],
            ["mat_ID_j", self.mat_ID_j],
            ["ele_orientation", self.ele_orientation],
            ["tranf_ID", self.geo_transf_ID],
            ["Initialized", self.Initialized]]

    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the member model (optional).
        """
        print("")
        print("Requested info for SpringBasedElement member model, ID = {}".format(self.element_ID))
        print("Section associated {} ".format(self.section_name_tag))
        print("Material model of the spring i, ID = {}".format(self.mat_ID_i))
        print("Material model of the spring j, ID = {}".format(self.mat_ID_j))
        print("Area A = {} mm2".format(self.A/mm2_unit))
        print("Young modulus E = {} GPa".format(self.E/GPa_unit))
        print("n modified moment of inertia Iy_mod = {} mm4".format(self.Iy_mod/mm4_unit))
        print("Geometric transformation = {}".format(self.geo_transf_ID))
        print("")

        if plot:
            if self.Initialized:
                plot_member(np.array(self.element_array))
                if block:
                    plt.show()
            else:
                print("The SpringBasedElement is not initialized (node and elements not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        self.element_array = [[self.element_ID, self.iNode_ID_spring, self.jNode_ID_spring]]
        
        if self.mat_ID_i != -1:
            # Define zero length element i
            node(self.iNode_ID_spring, *nodeCoord(self.iNode_ID))
            iSpring_ID = IDConvention(self.iNode_ID, self.iNode_ID_spring)
            RotationalSpring(iSpring_ID, self.iNode_ID, self.iNode_ID_spring, self.mat_ID_i)
            self.element_array.append([iSpring_ID, self.iNode_ID, self.iNode_ID_spring])
        if self.mat_ID_j != -1:
            # Define zero length element j
            node(self.jNode_ID_spring, *nodeCoord(self.jNode_ID))
            jSpring_ID = IDConvention(self.jNode_ID, self.jNode_ID_spring)
            RotationalSpring(jSpring_ID, self.jNode_ID, self.jNode_ID_spring, self.mat_ID_j)
            self.element_array.append([jSpring_ID, self.jNode_ID, self.jNode_ID_spring])
        
        # Define element
        element("elasticBeamColumn", self.element_ID, self.iNode_ID_spring, self.jNode_ID_spring, self.A, self.E, self.Iy_mod, self.geo_transf_ID)

        # Update class
        self.Initialized = True
        self.UpdateStoredData()


class SpringBasedElementSteelIShape(SpringBasedElement):
    # L_b = assumed the same for top and bottom springs
    def __init__(self, iNode_ID: int, jNode_ID: int, section: SteelIShape, geo_transf_ID: int, mat_ID_i=-1, mat_ID_j=-1):
        self.section = section
        if mat_ID_i != -1 and mat_ID_i < 0: raise NegativeValue()
        if mat_ID_j != -1 and mat_ID_j < 0: raise NegativeValue()
        if mat_ID_i == -1 and mat_ID_j == -1: raise NameError("No springs defined for element ID = {}".format(IDConvention(iNode_ID, jNode_ID)))

        if mat_ID_i != -1 and mat_ID_j != -1:
            L_0 = section.L/2
        else:
            L_0 = section.L

        super().__init__(iNode_ID, jNode_ID, section.A, section.E, section.Iy_mod, geo_transf_ID, mat_ID_i=mat_ID_i, mat_ID_j=mat_ID_j)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()

class SpringBasedElementModifiedIMKSteelIShape(SpringBasedElement):
    # L_b = assumed the same for top and bottom springs
    def __init__(self, iNode_ID: int, jNode_ID: int, section: SteelIShape, geo_transf_ID: int, new_mat_ID_i=-1, new_mat_ID_j=-1, N_G = 0, L_b = -1):
        self.section = section
        if new_mat_ID_i != -1 and new_mat_ID_i < 0: raise NegativeValue()
        if new_mat_ID_j != -1 and new_mat_ID_j < 0: raise NegativeValue()
        if new_mat_ID_i == -1 and new_mat_ID_j == -1: raise NameError("No springs defined for element ID = {}".format(IDConvention(iNode_ID, jNode_ID)))

        if new_mat_ID_i != -1 and new_mat_ID_j != -1:
            L_0 = section.L/2
        else:
            L_0 = section.L

        if new_mat_ID_i != -1:
            # Create mat i
            iSpring = ModifiedIMKSteelIShape(new_mat_ID_i, section, N_G, L_0 = L_0, L_b = L_b)
            iSpring.Bilin()

        if new_mat_ID_j != -1:
            # Create mat j
            jSpring = ModifiedIMKSteelIShape(new_mat_ID_j, section, N_G, L_0 = L_0, L_b = L_b)
            jSpring.Bilin()

        super().__init__(iNode_ID, jNode_ID, section.A, section.E, section.Iy_mod, geo_transf_ID, mat_ID_i=new_mat_ID_i, mat_ID_j=new_mat_ID_j)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()


class ForceBasedElement(MemberModel):
    # Class that stores funcions and material properties of a force-based beam-column element.
    # Warning: the units should be m and N
    # Convention:																							
    # 		NodeID: 		xy 			with x = pier, y = floor 											o xy7	|	
    #		PlHingeID:		xya			with x = pier, y = floor, a:	 --o xy  xy2 o-----o xy3  xy o--	|		o xy
    #																										|	
    #																										|		o xy
    #																										o xy6	|
    #		ElementID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		TrussID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		PDeltaColID:	xy(a)xy(a)	with xy(a) = NodeID i and j
    #		Spring:			xy(a)xy(a)	with xy(a) = NodeID i and j
    
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, geo_transf_ID: int,
        new_integration_ID = -1, Ip = 5, integration_type = "Lobatto", max_iter = MAX_ITER_INTEGRATION, tol = TOL_INTEGRATION):
    
        # Check
        if iNode_ID < 1: raise NegativeValue()
        if jNode_ID < 1: raise NegativeValue()
        if fiber_ID < 1: raise NegativeValue()
        if geo_transf_ID < 1: raise NegativeValue()
        if new_integration_ID != -1 and new_integration_ID < 1: raise NegativeValue()
        if Ip != -1 and Ip < 3: raise NegativeValue()
        if max_iter < 0: raise NegativeValue()
        if tol < 0: raise NegativeValue()

        # Arguments
        self.iNode_ID = iNode_ID
        self.jNode_ID = jNode_ID
        self.fiber_ID = fiber_ID 
        self.geo_transf_ID = geo_transf_ID
        self.Ip = Ip
        self.integration_type = integration_type
        self.max_iter = max_iter
        self.tol = tol

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.Initialized = False
        self.ReInit(new_integration_ID)


    def ReInit(self, new_integration_ID):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Precompute some members
        self.element_ID = IDConvention(self.iNode_ID, self.jNode_ID)

        # Arguments
        self.new_integration_ID = self.element_ID if new_integration_ID == -1 else new_integration_ID

        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "ForceBasedElement"], # Tag for differentiating different data
            ["element_ID", self.element_ID],
            ["section_name_tag", self.section_name_tag],
            ["Ip", self.Ip],
            ["iNode_ID", self.iNode_ID],
            ["jNode_ID", self.jNode_ID],
            ["fiber_ID", self.fiber_ID],
            ["new_integration_ID", self.new_integration_ID],
            ["integration_type", self.integration_type],
            ["tol", self.tol],
            ["max_iter", self.max_iter],
            ["tranf_ID", self.geo_transf_ID],
            ["Initialized", self.Initialized]]


    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the member model (optional).
        """
        print("")
        print("Requested info for GIFBElement member model, ID = {}".format(self.element_ID))
        print("Fiber associated, ID = {} ".format(self.fiber_ID))
        print("Integration type '{}', ID = {}".format(self.integration_type, self.new_integration_ID))
        print("Section associated {} ".format(self.section_name_tag))
        print("Number of integration points along the element Ip = {}, max iter = {}, tol = {}".format(self.Ip, self.max_iter, self.tol))
        print("Geometric transformation = {}".format(self.geo_transf_ID))
        print("")

        if plot:
            if self.Initialized:
                plot_member(np.array(self.element_array))
                if block:
                    plt.show()
            else:
                print("The ForceBasedElement is not initialized (element not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        self.element_array = [[self.element_ID, self.iNode_ID, self.jNode_ID]]

        # Define integration type
        beamIntegration(self.integration_type, self.new_integration_ID, self.fiber_ID, self.Ip)
        
        # Define element
        element('forceBeamColumn', self.element_ID, self.iNode_ID, self.jNode_ID, self.geo_transf_ID, self.new_integration_ID, '-iter', self.max_iter, self.tol)
        
        # Update class
        self.Initialized = True
        self.UpdateStoredData()


class ForceBasedElementFibersRectRCRectShape(ForceBasedElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber: FibersRectRCRectShape, geo_transf_ID: int,
        new_integration_ID=-1, Ip=5, integration_type="Lobatto", max_iter=MAX_ITER_INTEGRATION, tol=TOL_INTEGRATION):
        self.section = fiber.section
        super().__init__(iNode_ID, jNode_ID, fiber.ID, geo_transf_ID,
            new_integration_ID=new_integration_ID, Ip=Ip, integration_type=integration_type, max_iter=max_iter, tol=tol)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()


class ForceBasedElementFibersCircRCCircShape(ForceBasedElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber: FibersCircRCCircShape, geo_transf_ID: int,
        new_integration_ID=-1, Ip=5, integration_type="Lobatto", max_iter=MAX_ITER_INTEGRATION, tol=TOL_INTEGRATION):
        self.section = fiber.section
        super().__init__(iNode_ID, jNode_ID, fiber.ID, geo_transf_ID,
            new_integration_ID=new_integration_ID, Ip=Ip, integration_type=integration_type, max_iter=max_iter, tol=tol)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()


class ForceBasedElementFibersIShapeSteelIShape(ForceBasedElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber: FibersIShapeSteelIShape, geo_transf_ID: int,
        new_integration_ID=-1, Ip=5, integration_type="Lobatto", max_iter=MAX_ITER_INTEGRATION, tol=TOL_INTEGRATION):
        self.section = fiber.section
        super().__init__(iNode_ID, jNode_ID, fiber.ID, geo_transf_ID,
            new_integration_ID=new_integration_ID, Ip=Ip, integration_type=integration_type, max_iter=max_iter, tol=tol)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()    


class GIFBElement(MemberModel):
    # Class that stores funcions and material properties of a Gradient-Inelastic Flexibility-based element.
    # Warning: the units should be m and N
    # Convention:																							
    # 		NodeID: 		xy 			with x = pier, y = floor 											o xy7	|	
    #		PlHingeID:		xya			with x = pier, y = floor, a:	 --o xy  xy2 o-----o xy3  xy o--	|		o xy
    #																										|	
    #																										|		o xy
    #																										o xy6	|
    #		ElementID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		TrussID:		xy(a)xy(a)	with xy(a) = NodeID i and j
    #		PDeltaColID:	xy(a)xy(a)	with xy(a) = NodeID i and j
    #		Spring:			xy(a)xy(a)	with xy(a) = NodeID i and j
    
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, D_bars, fy, geo_transf_ID: int, 
        lambda_i = -1, lambda_j = -1, Lp = -1, Ip = -1, new_integration_ID = -1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION):
    
        # Check
        if iNode_ID < 1: raise NegativeValue()
        if jNode_ID < 1: raise NegativeValue()
        if fiber_ID < 1: raise NegativeValue()
        if D_bars < 0: raise NegativeValue()
        if fy < 0: raise NegativeValue()
        if geo_transf_ID < 1: raise NegativeValue()
        if lambda_i != -1 and lambda_i < 0: raise NegativeValue()
        if lambda_j != -1 and lambda_j < 0: raise NegativeValue()
        if lambda_i == 0 and lambda_j == 0: raise print("!!!!!!! WARNING !!!!!!! No plastic length defined for element ID = {}".format(IDConvention(iNode_ID, jNode_ID)))
        if Lp != -1 and Lp < 0: raise NegativeValue()
        if Ip != -1 and Ip < 3: raise NegativeValue()
        if new_integration_ID != -1 and new_integration_ID < 1: raise NegativeValue()
        if min_tol < 0: raise NegativeValue()
        if max_tol < 0: raise NegativeValue()
        if max_iter < 0: raise NegativeValue()

        # Arguments
        self.iNode_ID = iNode_ID
        self.jNode_ID = jNode_ID
        self.D_bars = D_bars
        self.fy = fy
        self.geo_transf_ID = geo_transf_ID
        self.fiber_ID = fiber_ID 
        self.min_tol = min_tol
        self.max_tol = max_tol
        self.max_iter = max_iter

        # Initialized the parameters that are dependent from others
        self.section_name_tag = "None"
        self.Initialized = False
        self.ReInit(lambda_i, lambda_j, Lp, Ip, new_integration_ID)

    def ReInit(self, lambda_i = -1, lambda_j = -1, Lp = -1, Ip = -1, new_integration_ID = -1):
        """Function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "deepcopy()" from the module "copy".
        """
        # Precompute some members
        iNode = np.array(nodeCoord(self.iNode_ID))
        jNode = np.array(nodeCoord(self.jNode_ID))
        self.L = np.linalg.norm(iNode-jNode)
        self.element_ID = IDConvention(self.iNode_ID, self.jNode_ID)

        # Arguments
        self.Lp = self.ComputeLp() if Lp == -1 else Lp
        self.Ip = self.ComputeIp() if Ip == -1 else Ip
        self.lambda_i = self.L/self.Lp if lambda_i == -1 else lambda_i
        self.lambda_j = self.L/self.Lp if lambda_j == -1 else lambda_j
        self.new_integration_ID = self.element_ID if new_integration_ID == -1 else new_integration_ID

        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        self.data = [["INFO_TYPE", "GIFBElement"], # Tag for differentiating different data
            ["element_ID", self.element_ID],
            ["section_name_tag", self.section_name_tag],
            ["L", self.L],
            ["D_bars", self.D_bars],
            ["fy", self.fy],
            ["Lp", self.Lp],
            ["Ip", self.Ip],
            ["iNode_ID", self.iNode_ID],
            ["lambda_i", self.lambda_i],
            ["jNode_ID", self.jNode_ID],
            ["lambda_j", self.lambda_j],
            ["fiber_ID", self.fiber_ID],
            ["new_integration_ID", self.new_integration_ID],
            ["min_tol", self.min_tol],
            ["max_tol", self.max_tol],
            ["max_iter", self.max_iter],
            ["tranf_ID", self.geo_transf_ID],
            ["Initialized", self.Initialized]]


    def ShowInfo(self, plot = False, block = False):
        """Function that show the data stored in the class in the command window and plots the member model (optional).
        """
        print("")
        print("Requested info for GIFBElement member model, ID = {}".format(self.element_ID))
        print("Fiber associated, ID = {} ".format(self.fiber_ID))
        print("Integration type 'Simpson', ID = {}".format(self.new_integration_ID))
        print("Section associated {} ".format(self.section_name_tag))
        print("Length L = {} m".format(self.L/m_unit))
        print("Diameter of the reinforcing bars D_bars = {} mm2".format(self.D_bars/mm2_unit))
        print("Reinforcing bar steel strength fy = {} MPa".format(self.fy/MPa_unit))
        print("Plastic length Lp = {} mm".format(self.Lp/mm_unit))
        print("Number of integration points along the element Ip = {}, max iter = {}, (min, max tol) = ({},{})".format(self.Ip, self.max_iter, self.min_tol, self.max_tol))
        print("Lambda_i = {} and lambda_j = {}".format(self.lambda_i, self.lambda_j))
        print("Geometric transformation = {}".format(self.geo_transf_ID))
        print("")

        if plot:
            if self.Initialized:
                plot_member(np.array(self.element_array))
                if block:
                    plt.show()
            else:
                print("The GIFBElement is not initialized (element not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        self.element_array = [[self.element_ID, self.iNode_ID, self.jNode_ID]]

        # Define integration type
        beamIntegration('Simpson', self.new_integration_ID, self.fiber_ID, self.Ip)
        
        # Define element
        element('gradientInelasticBeamColumn', self.element_ID, self.iNode_ID, self.jNode_ID, self.geo_transf_ID,
            self.new_integration_ID, self.lambda_i, self.lambda_j, self.Lp, '-iter', self.max_iter, self.min_tol, self.max_tol)
        
        # Update class
        self.Initialized = True
        self.UpdateStoredData()

    def ComputeLp(self):
        """Compute the plastic length using Paulay 1992.

        Returns:
            double: Plastic length
        """
        return (0.08*self.L/m_unit + 0.022*self.D_bars/m_unit*self.fy/MPa_unit)*m_unit


    def ComputeIp(self):
        """Compute the number of integration points with equal distance along the element. For more information, see Salehi and Sideris 2020.

        Returns:
            int: Number of integration points
        """
        tmp = math.ceil(1.5*self.L/self.Lp + 1)
        if (tmp % 2) == 0:
            return tmp + 1
        else:
            return tmp


class GIFBElementRCRectShape(GIFBElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, section: RCRectShape, geo_transf_ID: int,
        lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION):
        self.section = section
        super().__init__(iNode_ID, jNode_ID, fiber_ID, section.D_bars, section.fy, geo_transf_ID,
            lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
            min_tol=min_tol, max_tol=max_tol, max_iter=max_iter)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()


class GIFBElementFibersRectRCRectShape(GIFBElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, fib: FibersRectRCRectShape, geo_transf_ID: int,
        lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION):
        self.section = fib.section
        super().__init__(iNode_ID, jNode_ID, fib.ID, self.section.D_bars, self.section.fy, geo_transf_ID,
        lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
        min_tol=min_tol, max_tol=max_tol, max_iter=max_iter)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()


class GIFBElementRCCircShape(GIFBElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, section: RCCircShape, geo_transf_ID: int,
        lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION):
        self.section = section
        super().__init__(iNode_ID, jNode_ID, fiber_ID, section.D_bars, section.fy, geo_transf_ID,
        lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
        min_tol=min_tol, max_tol=max_tol, max_iter=max_iter)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()


class GIFBElementFibersCircRCCircShape(GIFBElement):
    def __init__(self, iNode_ID: int, jNode_ID: int, fib: FibersCircRCCircShape, geo_transf_ID: int,
    lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
    min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION):
        self.section = fib.section
        super().__init__(iNode_ID, jNode_ID, fib.ID, self.section.D_bars, self.section.fy, geo_transf_ID,
        lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
        min_tol=min_tol, max_tol=max_tol, max_iter=max_iter)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()
