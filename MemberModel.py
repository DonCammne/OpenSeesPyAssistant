"""
Module for the member model.
Carmine Schipani, 2021
"""

from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from abc import abstractmethod
from copy import copy, deepcopy
from OpenSeesPyAssistant.Section import *
from OpenSeesPyAssistant.DataManagement import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *
from OpenSeesPyAssistant.Constants import *
from OpenSeesPyAssistant.Fibers import *
from OpenSeesPyAssistant.Connections import *
from OpenSeesPyAssistant.FunctionalFeatures import *

class MemberModel(DataManagement):
    """
    Parent abstract class for the storage and manipulation of a member's information (mechanical and geometrical parameters, etc) and the initialisation in the model.

    @param DataManagement: Parent abstract class.
    """
    @abstractmethod
    def Record(self, ele_ID, name_txt: str, data_dir: str, force_rec = True, def_rec = True, time_rec = True):
        """
        Abstract method that records the forces, deformation and time of the member associated with the class.

        @param ele_ID (int): The ID of the element that will be recorded.
        @param name_txt (str): Name of the recorded data (no .txt).
        @param data_dir (str): Directory for the storage of data.
        @param force_rec (bool, optional): Option to record the forces (Fx, Fy, Mz). Defaults to True.
        @param def_rec (bool, optional): Option to record the deformation (theta) for ZeroLength element. Defaults to True.
        @param time_rec (bool, optional): Option to record time. Defaults to True.
        """
        if self.Initialized:
            if not os.path.exists(data_dir):
                print("Folder {} not found in this directory; creating one".format(data_dir))
                os.makedirs(data_dir)
            
            if time_rec:
                if force_rec:
                    recorder("Element", "-file", '{}/{}.txt'.format(data_dir, name_txt), "-time", "-ele", ele_ID, "force")
                if def_rec:
                    recorder("Element", "-file", '{}/{}.txt'.format(data_dir, name_txt), "-time", "-ele", ele_ID, "deformation")
            else:
                if force_rec:
                    recorder("Element", "-file", '{}/{}.txt'.format(data_dir, name_txt), "-ele", ele_ID, "force")
                if def_rec:
                    recorder("Element", "-file", '{}/{}.txt'.format(data_dir, name_txt), "-ele", ele_ID, "deformation")
        else:
                print("The element is not initialized (node and/or elements not created), ID = {}".format(ele_ID))
    
    @abstractmethod
    def RecordNodeDef(self, iNode_ID: int, jNode_ID: int, name_txt: str, data_dir: str, time_rec = True):
        """
        Abstract method that records the deformation and time of the member's nodes associated with the class.

        @param iNode_ID (int): ID of the node i.
        @param jNode_ID (int): ID of the node j.
        @param name_txt (str): Name of the recorded data (no .txt).
        @param data_dir (str): Directory for the storage of data.
        @param time_rec (bool, optional): Option to record time. Defaults to True.
        """
        if self.Initialized:
            if not os.path.exists(data_dir):
                print("Folder {} not found in this directory; creating one".format(data_dir))
                os.makedirs(data_dir)
            
            if time_rec:
                recorder("Node", "-file", '{}/{}.txt'.format(data_dir, name_txt), "-time", "-node", iNode_ID, jNode_ID, "-dof", 1, 2, 3, "disp")
            else:
                recorder("Node", "-file", '{}/{}.txt'.format(data_dir, name_txt), "-node", iNode_ID, jNode_ID, "-dof", 1, 2, 3, "disp")
        else:
                print("The element is not initialized (node and/or elements not created), iNode ID = {}, jNode ID = {}".format(iNode_ID, jNode_ID))
    

    def _CheckL(self):
        """
        Private abstract method to check if the length of the line member is the same (with 1 cm of tolerance) with the length defined in the section used.
        """
        iNode = np.array(nodeCoord(self.iNode_ID))
        jNode = np.array(nodeCoord(self.jNode_ID))
        L = np.linalg.norm(iNode-jNode)
        if abs(L-self.section.L) > 1*cm_unit:
            print("!!!!!!! WARNING !!!!!!! The length declared in the section name '{}' (L={} m) is different from thelength of the element associated (ID={}, L ={}m)".format(
                    self.section.name_tag, L/m_unit, self.element_ID, self.section.L/m_unit))


class PanelZone(MemberModel):
    """
	Class that handles the storage and manipulation of a panel zone's information (mechanical and geometrical parameters, etc) and the initialisation in the model.

    @param MemberModel: Parent abstract class.
    """
    def __init__(self, master_node_ID: int, mid_panel_zone_width, mid_panel_zone_height, E, A_rigid, I_rigid, geo_transf_ID: int, mat_ID: int, pin_corners = True):
        """
        Constructor of the class.

        @param master_node_ID (int): ID of the master node (central top node that should be a grid node).
        @param mid_panel_zone_width (float): Mid panel zone width.
        @param mid_panel_zone_height (float): Mid panel zone height.
        @param E (float): Young modulus.
        @param A_rigid (float): A very rigid area.
        @param I_rigid (float): A very rigid moment of inertia.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param mat_ID (int): ID of the material model for the panel zone spring.
        @param pin_corners (bool, optional): Option to pin the corners (xy03/xy04, xy06/xy07, xy09/xy10) or not. Used for RCS models. Defaults to True.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: mid_panel_zone_width needs to be positive.
        @exception NegativeValue: mid_panel_zone_height needs to be positive.
        @exception NegativeValue: E needs to be positive.
        @exception NegativeValue: A_rigid needs to be positive.
        @exception NegativeValue: I_rigid needs to be positive.
        @exception NegativeValue: geo_tranf_ID needs to be a positive integer.
        @exception NegativeValue: mat_ID needs to be a positive integer.
        """
        # Check
        if master_node_ID < 1: raise NegativeValue()
        # if master_node_ID > 99: raise WrongNodeIDConvention(master_node_ID)
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
        self.pin_corners = pin_corners

        # Initialized the parameters that are dependent from others
        self.col_section_name_tag = "None"
        self.beam_section_name_tag = "None"
        self.Initialized = False
        self.ReInit()


    def ReInit(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        # Arguments
        self.spring_ID = -1

        # Members
        if self.col_section_name_tag != "None": self.col_section_name_tag = self.col_section_name_tag + " (modified)"
        if self.beam_section_name_tag != "None": self.beam_section_name_tag = self.beam_section_name_tag + " (modified)"

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
        self.data = [["INFO_TYPE", "PanelZone"], # Tag for differentiating different data
            ["master_node_ID", self.master_node_ID],
            ["col_section_name_tag", self.col_section_name_tag],
            ["beam_section_name_tag", self.beam_section_name_tag],
            ["mat_ID", self.mat_ID],
            ["spring_ID", self.spring_ID],
            ["mid_panel_zone_width", self.mid_panel_zone_width],
            ["mid_panel_zone_height", self.mid_panel_zone_height],
            ["E", self.E],
            ["A_rigid", self.A_rigid],
            ["I_rigid", self.I_rigid],
            ["tranf_ID", self.geo_transf_ID],
            ["Initialized", self.Initialized]]


    def ShowInfo(self, plot = False, block = False):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the material model. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
        """
        print("")
        print("Requested info for Panel Zone member model, master node ID = {}".format(self.master_node_ID))
        print("Section associated, column: {} ".format(self.col_section_name_tag))
        print("Section associated, beam: {} ".format(self.beam_section_name_tag))
        print("Material model of the panel zone ID = {}".format(self.mat_ID))
        print("Spring ID = {} (if -1, not defined yet)".format(self.spring_ID))
        print("Mid panel zone width = {} mm".format(self.mid_panel_zone_width/mm_unit))
        print("Mid panel zone height = {} mm".format(self.mid_panel_zone_height/mm_unit))
        print("Young modulus E = {} GPa".format(self.E/GPa_unit))
        print("Area of the elements (rigid) = {} mm2".format(self.A_rigid/mm2_unit))
        print("Moment of inetia of the elements (strong axis, rigid) = {} mm4".format(self.I_rigid/mm4_unit))
        print("Geometric transformation = {}".format(self.geo_transf_ID))
        print("")

        if plot:
            if self.Initialized:
                plot_member(self.element_array, "Panel zone, ID = {}".format(self.master_node_ID))
                if block:
                    plt.show()
            else:
                print("The panel zone is not initialized (node and elements not created) for master node ID = {}".format(self.master_node_ID))


    def CreateMember(self):
        """
        Method that initialises the member by calling the OpenSeesPy commands through various functions.
        """
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
        self.spring_ID = IDConvention(xy1, xy01)
        RotationalSpring(self.spring_ID, xy1, xy01, self.mat_ID)
        self.element_array.append([self.spring_ID, xy1, xy01])
        self.iNode_ID = xy1
        self.jNode_ID = xy01

        # Pin connections
        if self.pin_corners:
            Pin(xy03, xy04)
            Pin(xy06, xy07)
            Pin(xy09, xy10)

        # Update class
        self.Initialized = True
        self.UpdateStoredData()


    def Record(self, name_txt: str, data_dir: str, force_rec=True, def_rec=True, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().Record(self.spring_ID, name_txt, data_dir, force_rec=force_rec, def_rec=def_rec, time_rec=time_rec)


    def RecordNodeDef(self, name_txt: str, data_dir: str, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().RecordNodeDef(self.iNode_ID, self.jNode_ID, name_txt, data_dir, time_rec=time_rec)


    def _CheckL(self):
        """
        (placeholder). No applicable for the panel zone. 
        """
        print("No length check for panel zone")


class PanelZoneSteelIShape(PanelZone):
    """
    Class that is the children of PanelZone and combine the class SteelIShape (section) to retrieve the information needed.  

    @param PanelZone: Parent class.
    """
    def __init__(self, master_node_ID: int, col: SteelIShape, beam: SteelIShape, geo_transf_ID: int, mat_ID: int, rigid = RIGID):
        """
        Constructor of the class.

        @param master_node_ID (int): ID of the master node (central top node that should be a grid node).
        @param col (SteelIShape): SteelIShape column section object.
        @param beam (SteelIShape): SteelIShape beam section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param mat_ID (int): ID of the material model for the panel zone spring.
        @param rigid (float, optional): Parameter with a value enough big to assure rigidity of one element
            but enough small to avoid convergence problem. Defaults to RIGID.
        """
        self.col = deepcopy(col)
        self.beam = deepcopy(beam)
        super().__init__(master_node_ID, col.d/2.0, beam.d/2.0, col.E, max(col.A, beam.A)*rigid, max(col.Iy, beam.Iy)*rigid, geo_transf_ID, mat_ID)

        self.col_section_name_tag = col.name_tag
        self.beam_section_name_tag = beam.name_tag
        self.UpdateStoredData()


class PanelZoneRCS(PanelZone):
    """
    WIP: Class that is the children of PanelZone and it's used for the panel zone in a RCS (RC column continous, Steel beam).
    Note that the corners are not pinned (do it manually).

    @param PanelZone: Parent class.
    """
    def __init__(self, master_node_ID: int, col: RCRectShape, beam: SteelIShape, geo_transf_ID: int, mat_ID: int, rigid = RIGID):
        """
        Constructor of the class.

        @param master_node_ID (int): ID of the master node (central top node that should be a grid node).
        @param col (RCRectShape): RCRectShape column section object.
        @param beam (SteelIShape): SteelIShape beam section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param mat_ID (int): ID of the material model for the panel zone spring.
        @param rigid (float, optional): Parameter with a value enough big to assure rigidity of one element
            but enough small to avoid convergence problem. Defaults to RIGID.
        """
        self.col = deepcopy(col)
        self.beam = deepcopy(beam)
        super().__init__(master_node_ID, col.d/2.0, beam.d/2.0, beam.E, max(col.A, beam.A)*rigid, max(col.Iy, beam.Iy)*rigid, geo_transf_ID, mat_ID, False)

        self.col_section_name_tag = col.name_tag
        self.beam_section_name_tag = beam.name_tag
        self.UpdateStoredData()


class PanelZoneSteelIShapeGupta1999(PanelZoneSteelIShape):
    """
    Class that is the children of PanelZoneSteelIShape and automatically create the spring material model Gupta 1999 (ID = master_node_ID).

    @param PanelZoneSteelIShape: Parent class.
    """
    def __init__(self, master_node_ID: int, col: SteelIShape, beam: SteelIShape, geo_transf_ID: int, t_dp = 0, rigid=RIGID):
        """
        Constructor of the class.

        @param master_node_ID (int): ID of the master node (central top node that should be a grid node).
        @param col (SteelIShape): SteelIShape column section object.
        @param beam (SteelIShape): SteelIShape beam section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param t_dp (float, optional): Doubler plate thickness. Defaults to 0.
        @param rigid (float, optional): Parameter with a value enough big to assure rigidity of one element
            but enough small to avoid convergence problem. Defaults to RIGID.
        """
        self.col = deepcopy(col)
        self.beam = deepcopy(beam)
        mat_ID = master_node_ID
        pz_spring = Gupta1999SteelIShape(mat_ID, col, beam, t_dp)
        pz_spring.Hysteretic()

        super().__init__(master_node_ID, col, beam, geo_transf_ID, mat_ID, rigid)


class PanelZoneSteelIShapeSkiadopoulos2021(PanelZoneSteelIShape):
    """
    Class that is the children of PanelZoneSteelIShape and automatically create the spring material model Skiadopoulos 2021 (ID = master_node_ID).

    @param PanelZoneSteelIShape: Parent class.
    """
    def __init__(self, master_node_ID: int, col: SteelIShape, beam: SteelIShape, geo_transf_ID: int, t_dp = 0, rigid=RIGID):
        """
        Constructor of the class.

        @param master_node_ID (int): ID of the master node (central top node that should be a grid node).
        @param col (SteelIShape): SteelIShape column section object.
        @param beam (SteelIShape): SteelIShape beam section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param t_dp (float, optional): Doubler plate thickness. Defaults to 0.
        @param rigid (float, optional): Parameter with a value enough big to assure rigidity of one element
            but enough small to avoid convergence problem. Defaults to RIGID.
        """
        self.col = deepcopy(col)
        self.beam = deepcopy(beam)
        mat_ID = master_node_ID
        pz_spring = Skiadopoulos2021SteelIShape(mat_ID, col, beam, t_dp)
        pz_spring.Hysteretic()

        super().__init__(master_node_ID, col, beam, geo_transf_ID, mat_ID, rigid)


def DefinePanelZoneNodes(MasterNodeID: int, MidPanelZoneWidth, MidPanelZoneHeight):
    """
    Function that defines the remaining 10 nodes of a panel zone given the dimensions and the master node (top center one).
    ID convention for the panel zone: \n
    		PZNodeID:		12 nodes: top right 1xy (master), 1xy1 top right,					1xy09,1xy10     1xy      1xy1,1xy01 \n					
    						clockwise 10 nodes xy01-xy10 (with double node at corners)				o-----------o-----------o		\n
    						Spring at node 1xy1														|						|       \n
    		PZElemeneID:	8 elements: starting at node 1xy, clockwise								|						|       \n
    						(see function DefinePanelZoneElements for more info)					|						|       \n
    																								|						|       \n
    																						  1xy08 o						o 1xy02 \n
    																								|						|       \n
    																								|						|       \n
    																								|						|       \n
    																								|						|       \n
    																								o-----------o-----------o       \n
    																							1xy06,1xy07	   1xy05  	1xy03,1xy04 \n
	    Note that the top right node is defined differently because is where the spring is.

    @param MasterNodeID (int): ID of the master node (central top node that should be a grid node).
    @param MidPanelZoneWidth (float): Mid panel zone width.
    @param MidPanelZoneHeight (float): Mid panel zone height.

    """
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
    """
    Function that defines the 8 panel zone elements. For the ID convention, see DefinePanelZoneNodes.

    @param MasterNodeID (int): ID of the master node (central top node that should be a grid node).
    @param E (float): Young modulus.
    @param RigidA (float): A very rigid area.
    @param RigidI (float): A very rigid moment of inertia.
    @param TransfID (int): The geometric transformation (for more information, see OpenSeesPy documentation).

    @returns list: List of lists, wth each list containing the ID of the element, of node i and node j.
    """
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
    # if MasterNodeID > 99:
    #     print("Warning, convention: MasterNodeID's digits should be 2")

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
    """
	Class that handles the storage and manipulation of a elastic element's information (mechanical and geometrical parameters, etc) and the initialisation in the model.

    @param MemberModel: Parent abstract class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, A, E, Iy, geo_transf_ID: int, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param A (float): Area of the member.
        @param E (float): Young modulus.
        @param Iy (float): Second moment of inertia (strong axis).
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: A needs to be positive.
        @exception NegativeValue: E needs to be positive.
        @exception NegativeValue: Iy needs to be positive.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        """
        # Check
        if iNode_ID < 1: raise NegativeValue()
        if jNode_ID < 1: raise NegativeValue()
        if A < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if Iy < 0: raise NegativeValue()
        if geo_transf_ID < 1: raise NegativeValue()
        if ele_ID != -1 and ele_ID < 1: raise NegativeValue()

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
        self.ReInit(ele_ID = -1)

    def ReInit(self, ele_ID = -1):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.

        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # element ID
        self.element_ID = IDConvention(self.iNode_ID, self.jNode_ID) if ele_ID == -1 else ele_ID

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
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
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the material model. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
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
                plot_member(self.element_array, "Elastic Element, ID = {}".format(self.element_ID))
                if block:
                    plt.show()
            else:
                print("The ElasticElement is not initialized (node and elements not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        """
        Method that initialises the member by calling the OpenSeesPy commands through various functions.
        """
        self.element_array = [[self.element_ID, self.iNode_ID, self.jNode_ID]]
        
        # Define element
        element("elasticBeamColumn", self.element_ID, self.iNode_ID, self.jNode_ID, self.A, self.E, self.Iy, self.geo_transf_ID)

        # Update class
        self.Initialized = True
        self.UpdateStoredData()


    def Record(self, name_txt: str, data_dir: str, force_rec=True, def_rec=True, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().Record(self.element_ID, name_txt, data_dir, force_rec=force_rec, def_rec=def_rec, time_rec=time_rec)


    def RecordNodeDef(self, name_txt: str, data_dir: str, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().RecordNodeDef(self.iNode_ID, self.jNode_ID, name_txt, data_dir, time_rec=time_rec)


class ElasticElementSteelIShape(ElasticElement):
    """
    Class that is the children of ElasticElement and combine the class SteelIShape (section) to retrieve the information needed.  

    @param ElasticElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, section: SteelIShape, geo_transf_ID: int, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param section (SteelIShape): SteelIShape section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(section)
        super().__init__(iNode_ID, jNode_ID, section.A, section.E, section.Iy, geo_transf_ID, ele_ID)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class SpringBasedElement(MemberModel):
    """
	Class that handles the storage and manipulation of a spring-based element's information (mechanical and geometrical parameters, etc) and the initialisation in the model.

    @param MemberModel: Parent abstract class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, A, E, Iy_mod, geo_transf_ID: int, mat_ID_i = -1, mat_ID_j = -1, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param A (float): Area of the member.
        @param E (float): Young modulus.
        @param Iy_mod (float): Second moment of inertia (strong axis).
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param mat_ID_i (int, optional): ID of the material model for the spring in the node i (if present). Defaults to -1.
        @param mat_ID_j (int, optional): ID of the material model for the spring in the node j (if present). Defaults to -1.
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: A needs to be positive.
        @exception NegativeValue: E needs to be positive.
        @exception NegativeValue: Iy_mod needs to be positive.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer, if different from -1.
        @exception NegativeValue: ID needs to be a positive integer, if different from -1.
        @exception NameError: at least one spring needs to be defined.
        @exception NegativeValue: ID needs to be a positive integer, if different from -1.
        """
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
        if ele_ID != -1 and ele_ID < 0: raise NegativeValue()

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
        self.ReInit(ele_ID)


    def ReInit(self, ele_ID = -1):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        # orientation:
        self.ele_orientation = NodesOrientation(self.iNode_ID, self.jNode_ID)
        if self.ele_orientation == "zero_length": raise ZeroLength(IDConvention(self.iNode_ID, self.jNode_ID))
        
        if self.mat_ID_i != -1:
            self.iNode_ID_spring = OffsetNodeIDConvention(self.iNode_ID, self.ele_orientation, "i")
        else:
            self.iNode_ID_spring = self.iNode_ID
        
        if self.mat_ID_j != -1:
            self.jNode_ID_spring = OffsetNodeIDConvention(self.jNode_ID, self.ele_orientation, "j")
        else:
            self.jNode_ID_spring = self.jNode_ID
                
        # element ID
        self.element_ID = IDConvention(self.iNode_ID_spring, self.jNode_ID_spring) if ele_ID == -1 else ele_ID

        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
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
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the material model. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
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
                plot_member(self.element_array, "SpringBased Element, ID = {}".format(self.element_ID))
                if block:
                    plt.show()
            else:
                print("The SpringBasedElement is not initialized (node and elements not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        """
        Method that initialises the member by calling the OpenSeesPy commands through various functions.
        """
        self.element_array = [[self.element_ID, self.iNode_ID_spring, self.jNode_ID_spring]]
        if self.mat_ID_i != -1:
            # Define zero length element i
            node(self.iNode_ID_spring, *nodeCoord(self.iNode_ID))
            self.iSpring_ID = IDConvention(self.iNode_ID, self.iNode_ID_spring)
            RotationalSpring(self.iSpring_ID, self.iNode_ID, self.iNode_ID_spring, self.mat_ID_i)
            self.element_array.append([self.iSpring_ID, self.iNode_ID, self.iNode_ID_spring])
        if self.mat_ID_j != -1:
            # Define zero length element j
            node(self.jNode_ID_spring, *nodeCoord(self.jNode_ID))
            self.jSpring_ID = IDConvention(self.jNode_ID, self.jNode_ID_spring)
            RotationalSpring(self.jSpring_ID, self.jNode_ID, self.jNode_ID_spring, self.mat_ID_j)
            self.element_array.append([self.jSpring_ID, self.jNode_ID, self.jNode_ID_spring])
        
        # Define element
        element("elasticBeamColumn", self.element_ID, self.iNode_ID_spring, self.jNode_ID_spring, self.A, self.E, self.Iy_mod, self.geo_transf_ID)

        # Update class
        self.Initialized = True
        self.UpdateStoredData()


    def Record(self, spring_or_element: str, name_txt: str, data_dir: str, force_rec=True, def_rec=True, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        if spring_or_element == "element":
            super().Record(self.element_ID, name_txt, data_dir, force_rec=force_rec, def_rec=def_rec, time_rec=time_rec)
        elif spring_or_element == "spring_i":
            if self.mat_ID_i == -1:
                print("Spring i recorded not present in element ID = {}".format(self.element_ID))
            else:
                super().Record(self.iSpring_ID, name_txt, data_dir, force_rec=force_rec, def_rec=def_rec, time_rec=time_rec)
        elif spring_or_element == "spring_j":
            if self.mat_ID_j == -1:
                print("Spring j recorded not present in element ID = {}".format(self.element_ID))
            else:
                super().Record(self.jSpring_ID, name_txt, data_dir, force_rec=force_rec, def_rec=def_rec, time_rec=time_rec)
        else:
            print("No recording option with: '{}' with element ID: {}".format(spring_or_element, self.element_ID))
    

    def RecordNodeDef(self, name_txt: str, data_dir: str, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().RecordNodeDef(self.iNode_ID, self.jNode_ID, name_txt, data_dir, time_rec=time_rec)


class SpringBasedElementSteelIShape(SpringBasedElement):
    """
    Class that is the children of SpringBasedElement and combine the class SteelIShape (section) to retrieve the information needed.  
    L_b is assumed the same for top and bottom springs.

    @param SpringBasedElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, section: SteelIShape, geo_transf_ID: int, mat_ID_i=-1, mat_ID_j=-1, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param section (SteelIShape): SteelIShape section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param mat_ID_i (int, optional): ID of the material model for the spring in the node i (if present). Defaults to -1.
        @param mat_ID_j (int, optional): ID of the material model for the spring in the node j (if present). Defaults to -1.
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NameError: at least one spring needs to be defined.
        @exception NegativeValue: ID needs to be a positive integer.
        """
        self.section = deepcopy(section)
        if mat_ID_i != -1 and mat_ID_i < 0: raise NegativeValue()
        if mat_ID_j != -1 and mat_ID_j < 0: raise NegativeValue()
        if mat_ID_i == -1 and mat_ID_j == -1: raise NameError("No springs defined for element ID = {}".format(IDConvention(iNode_ID, jNode_ID)))
        if ele_ID != -1 and ele_ID < 0: raise NegativeValue()

        super().__init__(iNode_ID, jNode_ID, section.A, section.E, section.Iy_mod, geo_transf_ID, mat_ID_i=mat_ID_i, mat_ID_j=mat_ID_j, ele_ID=ele_ID)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class SpringBasedElementModifiedIMKSteelIShape(SpringBasedElement):
    """
    Class that is the children of SpringBasedElement and combine the class SteelIShape (section) to retrieve the information needed.  
    If there are two springs and the inflection point not in the middle, use two spring elements, connect them rigidly in the inflection point with one spring each in the extremes.
    L_b is assumed the same for top and bottom springs.

    @param SpringBasedElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, section: SteelIShape, geo_transf_ID: int, new_mat_ID_i=-1, new_mat_ID_j=-1,
        N_G = 0, L_0 = -1, L_b = -1, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param section (SteelIShape): SteelIShape section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param new_mat_ID_i (int, optional): New ID for the definition of the material model for the spring in the node i.
            If -1 is passed, the class generate no material model and no spring. If 0 is passed, no i spring. Defaults to -1.
        @param new_mat_ID_j (int, optional): New ID for the definition of the material model for the spring in the node j.
            If -1 is passed, the class generate no material model and no spring. If 0 is passed, no j spring. Defaults to -1.
        @param N_G (float, optional): Axial load. Defaults to 0.
        @param L_0 (float, optional): Distance from the maximal moment to zero. Defaults to -1, e.g. computed in __init__().
        @param L_b (float, optional): Maximal unbraced lateral buckling length. Defaults to -1, e.g. computed in __init__().
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NameError: at least one spring needs to be defined.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception ZeroLength: The two nodes are superimposed.
        """
        self.section = deepcopy(section)
        if new_mat_ID_i != -1 and new_mat_ID_i < 0: raise NegativeValue()
        if new_mat_ID_j != -1 and new_mat_ID_j < 0: raise NegativeValue()
        if new_mat_ID_i == 0 and new_mat_ID_j == 0: raise NameError("No springs imposed for element ID = {}. Use ElasticElement instead".format(IDConvention(iNode_ID, jNode_ID)))
        if ele_ID != -1 and ele_ID < 0: raise NegativeValue()

        if L_0 == -1:
            if new_mat_ID_i != 0 and new_mat_ID_j != 0:
                L_0 = section.L/2
            else:
                L_0 = section.L
        L_b = L_0 if L_b == -1 else L_b

        # auto assign ID for material of springs
        ele_orientation = NodesOrientation(iNode_ID, jNode_ID)
        if ele_orientation == "zero_length": raise ZeroLength(IDConvention(iNode_ID, jNode_ID))
        if new_mat_ID_i != 0 and new_mat_ID_i == -1:
            new_mat_ID_i = OffsetNodeIDConvention(iNode_ID, ele_orientation, "i")
        if new_mat_ID_j != 0 and new_mat_ID_j == -1:
            new_mat_ID_j = OffsetNodeIDConvention(jNode_ID, ele_orientation, "j")

        if new_mat_ID_i != 0:
            # Create mat i
            iSpring = ModifiedIMKSteelIShape(new_mat_ID_i, section, N_G, L_0 = L_0, L_b = L_b)
            iSpring.Bilin()

        if new_mat_ID_j != 0:
            # Create mat j
            jSpring = ModifiedIMKSteelIShape(new_mat_ID_j, section, N_G, L_0 = L_0, L_b = L_b)
            jSpring.Bilin()

        new_mat_ID_i = -1 if new_mat_ID_i == 0 else new_mat_ID_i
        new_mat_ID_j = -1 if new_mat_ID_j == 0 else new_mat_ID_j

        super().__init__(iNode_ID, jNode_ID, section.A, section.E, section.Iy_mod, geo_transf_ID, mat_ID_i=new_mat_ID_i, mat_ID_j=new_mat_ID_j, ele_ID=ele_ID)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class ForceBasedElement(MemberModel):
    """
	Class that handles the storage and manipulation of a force-based element's information (mechanical and geometrical parameters, etc) and the initialisation in the model.

    @param MemberModel: Parent abstract class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, geo_transf_ID: int,
        new_integration_ID = -1, Ip = 5, integration_type = "Lobatto", max_iter = MAX_ITER_INTEGRATION, tol = TOL_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fiber_ID (int): ID of the fiber section.
        @param geo_transf_ID (int): The geometric transformation (for more information, see OpenSeesPy documentation).
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param integration_type (str, optional): Integration type. FOr more information, see OpenSeesPy documentation.
            Defaults to "Lobatto".
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param tol (float, optional): Tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer, if different from -1.
        @exception NegativeValue: Ip needs to be a positive integer bigger than 3, if different from -1.
        @exception NegativeValue: max_iter needs to be a positive integer.
        @exception NegativeValue: tol needs to be positive.
        @exception NegativeValue: ID needs to be a positive integer, if different from -1.
        """
        # Check
        if iNode_ID < 1: raise NegativeValue()
        if jNode_ID < 1: raise NegativeValue()
        if fiber_ID < 1: raise NegativeValue()
        if geo_transf_ID < 1: raise NegativeValue()
        if new_integration_ID != -1 and new_integration_ID < 1: raise NegativeValue()
        if Ip != -1 and Ip < 3: raise NegativeValue()
        if max_iter < 0: raise NegativeValue()
        if tol < 0: raise NegativeValue()
        if ele_ID != -1 and ele_ID < 1: raise NegativeValue()

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
        self.ReInit(new_integration_ID, ele_ID)


    def ReInit(self, new_integration_ID, ele_ID = -1):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.

        @param new_integration_ID (int): ID of the integration technique.
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        # Precompute some members
        self.element_ID = IDConvention(self.iNode_ID, self.jNode_ID) if ele_ID == -1 else ele_ID

        # Arguments
        self.new_integration_ID = self.element_ID if new_integration_ID == -1 else new_integration_ID

        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
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
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the material model. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
        """
        print("")
        print("Requested info for ForceBasedElement member model, ID = {}".format(self.element_ID))
        print("Fiber associated, ID = {} ".format(self.fiber_ID))
        print("Integration type '{}', ID = {}".format(self.integration_type, self.new_integration_ID))
        print("Section associated {} ".format(self.section_name_tag))
        print("Number of integration points along the element Ip = {}, max iter = {}, tol = {}".format(self.Ip, self.max_iter, self.tol))
        print("Geometric transformation = {}".format(self.geo_transf_ID))
        print("")

        if plot:
            if self.Initialized:
                plot_member(self.element_array, "ForceBased Element, ID = {}".format(self.element_ID))
                if block:
                    plt.show()
            else:
                print("The ForceBasedElement is not initialized (element not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        """
        Method that initialises the member by calling the OpenSeesPy commands through various functions.
        """
        self.element_array = [[self.element_ID, self.iNode_ID, self.jNode_ID]]

        # Define integration type
        beamIntegration(self.integration_type, self.new_integration_ID, self.fiber_ID, self.Ip)
        
        # Define element
        element('forceBeamColumn', self.element_ID, self.iNode_ID, self.jNode_ID, self.geo_transf_ID, self.new_integration_ID, '-iter', self.max_iter, self.tol)
        
        # Update class
        self.Initialized = True
        self.UpdateStoredData()


    def Record(self, name_txt: str, data_dir: str, force_rec=True, def_rec=True, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().Record(self.element_ID, name_txt, data_dir, force_rec=force_rec, def_rec=def_rec, time_rec=time_rec)
    

    def RecordNodeDef(self, name_txt: str, data_dir: str, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().RecordNodeDef(self.iNode_ID, self.jNode_ID, name_txt, data_dir, time_rec=time_rec)


class ForceBasedElementFibersRectRCRectShape(ForceBasedElement):
    """
    Class that is the children of ForceBasedElement and combine the class FibersRectRCRectShape (fiber section) to retrieve the information needed.  

    @param ForceBasedElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber: FibersRectRCRectShape, geo_transf_ID: int,
        new_integration_ID=-1, Ip=5, integration_type="Lobatto", max_iter=MAX_ITER_INTEGRATION, tol=TOL_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fiber (FibersRectRCRectShape): FibersRectRCRectShape fiber section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param integration_type (str, optional): Integration type. FOr more information, see OpenSeesPy documentation.
            Defaults to "Lobatto".
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param tol (float, optional): Tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(fiber.section)
        super().__init__(iNode_ID, jNode_ID, fiber.ID, geo_transf_ID,
            new_integration_ID=new_integration_ID, Ip=Ip, integration_type=integration_type, max_iter=max_iter, tol=tol, ele_ID= ele_ID)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class ForceBasedElementFibersCircRCCircShape(ForceBasedElement):
    """
    Class that is the children of ForceBasedElement and combine the class FibersCircRCCircShape (fiber section) to retrieve the information needed.  

    @param ForceBasedElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber: FibersCircRCCircShape, geo_transf_ID: int,
        new_integration_ID=-1, Ip=5, integration_type="Lobatto", max_iter=MAX_ITER_INTEGRATION, tol=TOL_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fiber (FibersCircRCCircShape): FibersCircRCCircShape fiber section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param integration_type (str, optional): Integration type. FOr more information, see OpenSeesPy documentation.
            Defaults to "Lobatto".
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param tol (float, optional): Tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(fiber.section)
        super().__init__(iNode_ID, jNode_ID, fiber.ID, geo_transf_ID,
            new_integration_ID=new_integration_ID, Ip=Ip, integration_type=integration_type, max_iter=max_iter, tol=tol, ele_ID=ele_ID)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class ForceBasedElementFibersIShapeSteelIShape(ForceBasedElement):
    """
    Class that is the children of ForceBasedElement and combine the class FibersIShapeSteelIShape (fiber section) to retrieve the information needed.  

    @param ForceBasedElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber: FibersIShapeSteelIShape, geo_transf_ID: int,
        new_integration_ID=-1, Ip=5, integration_type="Lobatto", max_iter=MAX_ITER_INTEGRATION, tol=TOL_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fiber (FibersIShapeSteelIShape): FibersIShapeSteelIShape fiber section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param integration_type (str, optional): Integration type. FOr more information, see OpenSeesPy documentation.
            Defaults to "Lobatto".
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param tol (float, optional): Tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(fiber.section)
        super().__init__(iNode_ID, jNode_ID, fiber.ID, geo_transf_ID,
            new_integration_ID=new_integration_ID, Ip=Ip, integration_type=integration_type, max_iter=max_iter, tol=tol, ele_ID=ele_ID)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()   
        # Check length
        self._CheckL() 


class GIFBElement(MemberModel):
    """
	Class that handles the storage and manipulation of a Gradient-Inelastic Flexibility-based element's information
        (mechanical and geometrical parameters, etc) and the initialisation in the model.
    The integration technique is Simpson. For more information, see Sideris and Salehi 2016, 2017 and 2020.

    @param MemberModel: Parent abstract class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, D_bars, fy, geo_transf_ID: int, 
        lambda_i = -1, lambda_j = -1, Lp = -1, Ip = -1, new_integration_ID = -1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fiber_ID (int): ID of the fiber section.
        @param D_bars (float): Diameter of the vertical reinforcing bars.
        @param fy (float): Yield stress of the reinforcing bars.
        @param geo_transf_ID (int): The geometric transformation (for more information, see OpenSeesPy documentation).
        @param lambda_i (float, optional): Fraction of beam length over the plastic hinge length at end i (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end i.
        @param lambda_j (float, optional): Fraction of beam length over the plastic hinge length at end j (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end j.
        @param Lp (float, optional): Plastic hinge length. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param min_tol (float, optional): Minimal tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param max_tol (float, optional): Maximal tolerance for the integration convergence. Defaults to TOL_INTEGRATION*1e4.
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.

        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: D_bars needs to be positive.
        @exception NegativeValue: fy needs to be positive.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: lambda_i needs to be positive.
        @exception NegativeValue: lambda_j needs to be positive.
        @exception NegativeValue: No plastic length defined.
        @exception NegativeValue: Lp needs to be positive, if different from -1.
        @exception NegativeValue: Ip needs to be a positive integer bigger than 3, if different from -1.
        @exception NegativeValue: ID needs to be a positive integer.
        @exception NegativeValue: min_tol needs to be positive.
        @exception NegativeValue: max_tol needs to be positive.
        @exception NegativeValue: max_iter needs to be a positive integer.
        @exception NegativeValue: ID needs to be a positive integer, if different from -1.
        """
        # Check
        if iNode_ID < 1: raise NegativeValue()
        if jNode_ID < 1: raise NegativeValue()
        if fiber_ID < 1: raise NegativeValue()
        if D_bars < 0: raise NegativeValue()
        if fy < 0: raise NegativeValue()
        if geo_transf_ID < 1: raise NegativeValue()
        if lambda_i != -1 and lambda_i < 0: raise NegativeValue()
        if lambda_j != -1 and lambda_j < 0: raise NegativeValue()
        if lambda_i == 0 and lambda_j == 0: print("!!!!!!! WARNING !!!!!!! No plastic length defined for element ID = {}".format(IDConvention(iNode_ID, jNode_ID)))
        if Lp != -1 and Lp < 0: raise NegativeValue()
        if Ip != -1 and Ip < 3: raise NegativeValue()
        if new_integration_ID != -1 and new_integration_ID < 1: raise NegativeValue()
        if min_tol < 0: raise NegativeValue()
        if max_tol < 0: raise NegativeValue()
        if max_iter < 0: raise NegativeValue()
        if ele_ID != -1 and ele_ID < 0: raise NegativeValue()

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
        self.ReInit(lambda_i, lambda_j, Lp, Ip, new_integration_ID, ele_ID)

    def ReInit(self, lambda_i = -1, lambda_j = -1, Lp = -1, Ip = 5, new_integration_ID = -1, ele_ID = -1):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.

        @param lambda_i (float, optional): Fraction of beam length over the plastic hinge length at end i (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end i.
        @param lambda_j (float, optional): Fraction of beam length over the plastic hinge length at end j (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end j.
        @param Lp (float, optional): Plastic hinge length. Defaults to -1, e.g. computed here.
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        # Precompute some members
        iNode = np.array(nodeCoord(self.iNode_ID))
        jNode = np.array(nodeCoord(self.jNode_ID))
        self.L = np.linalg.norm(iNode-jNode)
        self.element_ID = IDConvention(self.iNode_ID, self.jNode_ID) if ele_ID == -1 else ele_ID

        # Arguments
        self.Lp = self.ComputeLp() if Lp == -1 else Lp
        self.Ip = self.ComputeIp() if Ip == -1 else Ip
        self.lambda_i = self.Lp/self.L if lambda_i == -1 else lambda_i
        self.lambda_j = self.Lp/self.L if lambda_j == -1 else lambda_j
        self.new_integration_ID = self.element_ID if new_integration_ID == -1 else new_integration_ID

        # Members
        if self.section_name_tag != "None": self.section_name_tag = self.section_name_tag + " (modified)"
        
        # Data storage for loading/saving
        self.UpdateStoredData()


    # Methods
    def UpdateStoredData(self):
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        """
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
        """
        Implementation of the homonym abstract method.
        See parent class DataManagement for detailed information.
        
        @param plot (bool, optional): Option to show the plot of the material model. Defaults to False.
	    @param block (bool, optional): Option to wait the user command 'plt.show()' (avoiding the stop of the program everytime that a plot should pop up). Defaults to False.
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
                plot_member(self.element_array, "GIFB Element, ID = {}".format(self.element_ID))
                if block:
                    plt.show()
            else:
                print("The GIFBElement is not initialized (element not created), ID = {}".format(self.element_ID))


    def CreateMember(self):
        """
        Method that initialises the member by calling the OpenSeesPy commands through various functions.
        """
        self.element_array = [[self.element_ID, self.iNode_ID, self.jNode_ID]]

        # Define integration type
        beamIntegration('Simpson', self.new_integration_ID, self.fiber_ID, self.Ip)
        
        # Define element TODO: Dr. Salehi: lambda useless
        element('gradientInelasticBeamColumn', self.element_ID, self.iNode_ID, self.jNode_ID, self.geo_transf_ID,
            self.new_integration_ID, self.lambda_i, self.lambda_j, self.Lp, '-iter', self.max_iter, self.min_tol, self.max_tol)
        
        # Update class
        self.Initialized = True
        self.UpdateStoredData()


    def Record(self, name_txt: str, data_dir: str, force_rec=True, def_rec=True, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().Record(self.element_ID, name_txt, data_dir, force_rec=force_rec, def_rec=def_rec, time_rec=time_rec)


    def RecordNodeDef(self, name_txt: str, data_dir: str, time_rec=True):
        """
        Implementation of the homonym abstract method.
        See parent class MemberModel for detailed information.
        """
        super().RecordNodeDef(self.iNode_ID, self.jNode_ID, name_txt, data_dir, time_rec=time_rec)


    def ComputeLp(self):
        """
        Method that computes the plastic length using Paulay 1992.

        @returns double: Plastic length
        """
        return (0.08*self.L/m_unit + 0.022*self.D_bars/m_unit*self.fy/MPa_unit)*m_unit


    def ComputeIp(self):
        """
        Compute the number of integration points with equal distance along the element. For more information, see Salehi and Sideris 2020.

        @returns int: Number of integration points
        """
        tmp = math.ceil(1.5*self.L/self.Lp + 1)
        if (tmp % 2) == 0:
            return tmp + 1
        else:
            return tmp


class GIFBElementRCRectShape(GIFBElement):
    """
    Class that is the children of GIFBElement and combine the class RCRectShape (section) to retrieve the information needed.  

    @param GIFBElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, section: RCRectShape, geo_transf_ID: int,
        lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fiber_ID (int): ID of the fiber section.
        @param section (RCRectShape): RCRectShape section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param lambda_i (float, optional): Fraction of beam length over the plastic hinge length at end i (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end i.
        @param lambda_j (float, optional): Fraction of beam length over the plastic hinge length at end j (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end j.
        @param Lp (float, optional): Plastic hinge length. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param min_tol (float, optional): Minimal tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param max_tol (float, optional): Maximal tolerance for the integration convergence. Defaults to TOL_INTEGRATION*1e4.
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(section)
        super().__init__(iNode_ID, jNode_ID, fiber_ID, section.D_bars, section.fy, geo_transf_ID,
            lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
            min_tol=min_tol, max_tol=max_tol, max_iter=max_iter, ele_ID = ele_ID)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class GIFBElementFibersRectRCRectShape(GIFBElement):
    """
    Class that is the children of GIFBElement and combine the class FibersRectRCRectShape (fiber section) to retrieve the information needed.  

    @param GIFBElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fib: FibersRectRCRectShape, geo_transf_ID: int,
        lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fib (FibersRectRCRectShape): FibersRectRCRectShape fiber section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param lambda_i (float, optional): Fraction of beam length over the plastic hinge length at end i (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end i.
        @param lambda_j (float, optional): Fraction of beam length over the plastic hinge length at end j (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end j.
        @param Lp (float, optional): Plastic hinge length. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param min_tol (float, optional): Minimal tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param max_tol (float, optional): Maximal tolerance for the integration convergence. Defaults to TOL_INTEGRATION*1e4.
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(fib.section)
        super().__init__(iNode_ID, jNode_ID, fib.ID, self.section.D_bars, self.section.fy, geo_transf_ID,
            lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
            min_tol=min_tol, max_tol=max_tol, max_iter=max_iter, ele_ID = ele_ID)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class GIFBElementRCCircShape(GIFBElement):
    """
    Class that is the children of GIFBElement and combine the class RCCircShape (section) to retrieve the information needed.  

    @param GIFBElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fiber_ID: int, section: RCCircShape, geo_transf_ID: int,
        lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fiber_ID (int): ID of the fiber section.
        @param section (RCCircShape): RCCircShape section object.
        @param geo_transf_ID (int): The geometric transformation (for more information, see OpenSeesPy documentation).
        @param lambda_i (float, optional): Fraction of beam length over the plastic hinge length at end i (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end i.
        @param lambda_j (float, optional): Fraction of beam length over the plastic hinge length at end j (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end j.
        @param Lp (float, optional): Plastic hinge length. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param min_tol (float, optional): Minimal tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param max_tol (float, optional): Maximal tolerance for the integration convergence. Defaults to TOL_INTEGRATION*1e4.
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(section)
        super().__init__(iNode_ID, jNode_ID, fiber_ID, section.D_bars, section.fy, geo_transf_ID,
            lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
            min_tol=min_tol, max_tol=max_tol, max_iter=max_iter, ele_ID = ele_ID)
        self.section_name_tag = section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()


class GIFBElementFibersCircRCCircShape(GIFBElement):
    """
    Class that is the children of GIFBElement and combine the class FibersCircRCCircShape (fiber section) to retrieve the information needed.  

    @param GIFBElement: Parent class.
    """
    def __init__(self, iNode_ID: int, jNode_ID: int, fib: FibersCircRCCircShape, geo_transf_ID: int,
        lambda_i=-1, lambda_j=-1, Lp=-1, Ip=-1, new_integration_ID=-1,
        min_tol = TOL_INTEGRATION, max_tol = TOL_INTEGRATION*1e4, max_iter = MAX_ITER_INTEGRATION, ele_ID = -1):
        """
        Constructor of the class.

        @param iNode_ID (int): ID of the first end node.
        @param jNode_ID (int): ID of the second end node.
        @param fib (FibersCircRCCircShape): FibersCircRCCircShape fiber section object.
        @param geo_transf_ID (int): A geometric transformation (for more information, see OpenSeesPy documentation).
        @param lambda_i (float, optional): Fraction of beam length over the plastic hinge length at end i (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end i.
        @param lambda_j (float, optional): Fraction of beam length over the plastic hinge length at end j (0 = no plastic hinge).
            Defaults to -1, e.g. plastic hinge in the end j.
        @param Lp (float, optional): Plastic hinge length. Defaults to -1, e.g. computed in ReInit().
        @param Ip (int, optional): Number of integration points (min. 3). Defaults to 5.
        @param new_integration_ID (int, optional): ID of the integration technique. Defaults to -1, e.g. computed in ReInit().
        @param min_tol (float, optional): Minimal tolerance for the integration convergence. Defaults to TOL_INTEGRATION (Units).
        @param max_tol (float, optional): Maximal tolerance for the integration convergence. Defaults to TOL_INTEGRATION*1e4.
        @param max_iter (int, optional): Maximal number of iteration to reach the integretion convergence. Defaults to MAX_ITER_INTEGRATION (Units).
        @param ele_ID (int, optional): Optional ID of the element. Defaults to -1, e.g. use IDConvention to define it.
        """
        self.section = deepcopy(fib.section)
        super().__init__(iNode_ID, jNode_ID, fib.ID, self.section.D_bars, self.section.fy, geo_transf_ID,
            lambda_i=lambda_i, lambda_j=lambda_j, Lp=Lp, Ip=Ip, new_integration_ID=new_integration_ID,
            min_tol=min_tol, max_tol=max_tol, max_iter=max_iter, ele_ID = ele_ID)
        self.section_name_tag = self.section.name_tag
        self.UpdateStoredData()
        # Check length
        self._CheckL()




