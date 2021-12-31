# Module with the different types of connections
#   Carmine Schipani, 2021

from openseespy.opensees import *
import numpy as np
import matplotlib.pyplot as plt
from OpenSeesPyAssistant.ErrorHandling import *

def RigidSupport(NodeID: int):
    """Function that fixes the x, y movements and the rotation of one node.

    Args:
        NodeID (int): ID of the node to be fixed
    """
    if NodeID < 1: raise NegativeValue()

    fix(NodeID, 1, 1, 1)


def Pin(NodeRID: int, NodeCID: int):
    """Function that constrains the translational DOF with a multi-point constraint.

    Args:
        NodeRID (int): Node ID which will be retained by the multi-point constraint
        NodeCID (int): Node ID which will be constrained by the multi-point constraint
    """
    if NodeCID == NodeRID: raise WrongArgument()
    if NodeRID < 1: raise NegativeValue()
    if NodeCID < 1: raise NegativeValue()

	# Constrain the translational DOF with a multi-point constraint
	#   		retained constrained DOF_1 DOF_2
    equalDOF(NodeRID, NodeCID, 1, 2)


def RotationalSpring(ElementID: int, NodeRID: int, NodeCID: int, MatID: int, Rigid = False):
    """Function that defines a uniaxial material spring and constrains the translations DOFs of the spring.

    Args:
        ElementID (int): ID of the zerolength element that models the spring
        NodeRID (int): Node ID which will be retained by the multi-point constraint
        NodeCID (int): Node ID which will be constrained by the multi-point constraint
        MatID (int): ID of the material model chosen
        Rigid (bool, optional): Optional argument that transforms the joint in a completely rigid connection. Defaults to False.
    """
    if ElementID < 1: raise NegativeValue()
    if NodeCID < 1: raise NegativeValue()
    if NodeRID < 1: raise NegativeValue()
    if NodeCID == NodeRID: raise WrongArgument()
    if MatID < 1: raise NegativeValue()

    if not Rigid:
    	# Zero length element (spring)
        element("zeroLength", ElementID, NodeRID, NodeCID, "-mat", MatID, "-dir", 6)
	    
        # Constrain the translational DOF with a multi-point constraint	
        Pin(NodeRID, NodeCID)
    else:
        equalDOF(NodeRID, NodeCID, 1, 2, 3)



	
