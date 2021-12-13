# Module with the member model
#   Carmine Schipani, 2021

from numpy.lib.function_base import append
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from abc import abstractmethod
from copy import copy, deepcopy
import openseespy.postprocessing.internal_plotting_functions as ipltf
import openseespy.postprocessing.Get_Rendering as opsplt
from OpenSeesPyHelper.Section import *
from OpenSeesPyHelper.DataManagement import *
from OpenSeesPyHelper.ErrorHandling import *
from OpenSeesPyHelper.Units import *
from OpenSeesPyHelper.Fibers import *
from OpenSeesPyHelper.Connections import *

# Member model
class MemberModel(DataManagement):
    pass

class PanelZone(MemberModel):
    # Class that stores funcions and material properties of a panel zone.
    # Warning: the units should be m and N
    # Warning: you need to define the spring element (node xy1, xy01)
    
    def __init__(self, master_node_ID: int, mid_panel_zone_width, mid_panel_zone_height, E, A_rigid, I_rigid, geo_transf_ID: int, mat_ID = -1):

        # Check
        if master_node_ID < 1: raise NegativeValue()
        if master_node_ID > 99: raise WrongNodeIDConvention(master_node_ID)
        if mid_panel_zone_width < 0: raise NegativeValue()
        if mid_panel_zone_height < 0: raise NegativeValue()
        if E < 0: raise NegativeValue()
        if A_rigid < 0: raise NegativeValue()
        if I_rigid < 0: raise NegativeValue()
        if geo_transf_ID > 1: raise NegativeValue()
        if mat_ID != -1 and mat_ID < 0: raise NegativeValue()

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
                plotMember(np.array(self.element_array))
                if block:
                    plt.show()
            else:
                print("The panel zone is not initialized (node and elements not created) for master node ID = {}".format(self.master_node_ID))


    def CreateMember(self):
        # Define nodes
        DefinePanelZoneNodes(self.master_node_ID, self.mid_panel_zone_width, self.mid_panel_zone_height)
        xy1 = self.master_node_ID*10+1
        xy01 = self.master_node_ID*100+1 
        xy03 = self.master_node_ID*100+3 
        xy04 = self.master_node_ID*100+4 
        xy06 = self.master_node_ID*100+6 
        xy07 = self.master_node_ID*100+7 
        xy09 = self.master_node_ID*100+9 
        xy10 = self.master_node_ID*100+10

        # Define rigid elements
        self.element_array = DefinePanelZoneElements(self.master_node_ID, self.E, self.A_rigid, self.I_rigid, self.geo_transf_ID)
        
        # Define zero length element
        if self.mat_ID == -1:
            print("Warning: the definition of the zero length element of the panel zone using these two nodes ({}, {}) has to be done by the user".format(xy1, xy01))
        else:
            spring_ID = xy1*10000+xy01
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
    def __init__(self, master_node_ID: int, col: SteelIShape, beam: SteelIShape, geo_transf_ID: int, mat_ID = -1, RIGID = 100.0):
        super().__init__(master_node_ID, col.d/2.0, beam.d/2.0, col.E, max(col.A, beam.A)*RIGID, max(col.Iy, beam.Iy)*RIGID, geo_transf_ID, mat_ID)

        self.col_section_name_tag = col.name_tag
        self.beam_section_name_tag = beam.name_tag
        self.UpdateStoredData()






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
    node(MasterNodeID*10+1, AxisCL+MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)
	# Convention: Two notes in the corners (already defined one, xy1) clockwise from xy01 to xy10
    node(MasterNodeID*100+1, AxisCL+MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)
    node(MasterNodeID*100+2, AxisCL+MidPanelZoneWidth, FloorCL)
    node(MasterNodeID*100+3, AxisCL+MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(MasterNodeID*100+4, AxisCL+MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(MasterNodeID*100+5, AxisCL, FloorCL-MidPanelZoneHeight)
    node(MasterNodeID*100+6, AxisCL-MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(MasterNodeID*100+7, AxisCL-MidPanelZoneWidth, FloorCL-MidPanelZoneHeight)
    node(MasterNodeID*100+8, AxisCL-MidPanelZoneWidth, FloorCL)
    node(MasterNodeID*100+9, AxisCL-MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)
    node(MasterNodeID*100+10, AxisCL-MidPanelZoneWidth, FloorCL+MidPanelZoneHeight)


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
    xy1 = MasterNodeID*10+1
    xy01 = MasterNodeID*100+1 
    xy02 = MasterNodeID*100+2 
    xy03 = MasterNodeID*100+3 
    xy04 = MasterNodeID*100+4 
    xy05 = MasterNodeID*100+5 
    xy06 = MasterNodeID*100+6 
    xy07 = MasterNodeID*100+7 
    xy08 = MasterNodeID*100+8 
    xy09 = MasterNodeID*100+9 
    xy10 = MasterNodeID*100+10 

    # Create element IDs using the convention: xy(a)xy(a)	with xy(a) = NodeID i and j
    #	Starting at MasterNodeID, clockwise
    if MasterNodeID > 99:
        raise NameError("MasterNodeID's digits should be 2")

    ele1 = xy*1000+xy1
    ele2 = xy01*10000+xy02
    ele3 = xy02*10000+xy03
    ele4 = xy04*10000+xy05
    ele5 = xy05*10000+xy06
    ele6 = xy07*10000+xy08
    ele7 = xy08*10000+xy09
    ele8 = xy10*100+xy


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


# Public functions
def plotMember(element_array: np.ndarray, show_element_ID = True, show_node_ID = True):
    ele_style = {'color':'black', 'linewidth':1, 'linestyle':'-'}
    node_style = {'color':'black', 'marker':'o', 'facecolor':'black','linewidth':0.}
    node_style_animation = {'color':'black', 'marker':'o','markersize':2., 'linewidth':0.} 

    node_text_style = {'fontsize':8, 'fontweight':'regular', 'color':'green'} 
    ele_text_style = {'fontsize':8, 'fontweight':'bold', 'color':'darkred'} 
    track_node = {}

    if show_element_ID:
        show_e_ID = 'yes'
    else:
        show_e_ID = 'no'

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    for ele in element_array:
        eleTag = int(ele[0])
        Nodes =ele[1:]
        
        if len(Nodes) == 2:
            # 2D element
            iNode = np.array(nodeCoord(Nodes[0].item()))
            jNode = np.array(nodeCoord(Nodes[1].item()))
            ipltf._plotBeam2D(iNode, jNode, ax, show_e_ID, eleTag, "solid")
            ax.scatter(*iNode, **node_style)
            ax.scatter(*jNode, **node_style)
            if show_node_ID:
                if abs(sum(iNode - jNode)) > 1e-6:
                    # beam-col
                    __plt_node(Nodes[0], track_node, iNode, ax, node_text_style)
                    __plt_node(Nodes[1], track_node, jNode, ax, node_text_style, h_align='right', v_align='bottom')
                else:
                    # zerolength
                    __plt_node(Nodes[0], track_node, iNode, ax, node_text_style)
                    __plt_node(Nodes[1], track_node, jNode, ax, node_text_style, h_align='right', v_align='bottom')
        else:
            print("Too many nodes in this elemnet (see shell elements)")
        
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    plt.axis('equal')

def __plt_node(nodeID: int, track_node: dict, NodeXY, ax, node_text_style, x_off = 0, y_off = 0, h_align = 'left', v_align='top'):
    if not nodeID in track_node:
        track_node[nodeID] = True
        ax.text(NodeXY[0]+x_off, NodeXY[1]+y_off, nodeID,**node_text_style, horizontalalignment=h_align, verticalalignment=v_align)
