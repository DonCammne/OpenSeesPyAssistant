"""Module with different functions useful when defining boundary conditions (fix) or connections (pin, rigid or springs). \n
Carmine Schipani, 2021
"""

from openseespy.opensees import *
from OpenSeesPyAssistant.ErrorHandling import *


def RigidSupport(NodeID: int):
    """
    Function that fixes the x, y movements and the rotation of one node.

    @param NodeID (int): ID of the node to be fixed

    @exception NegativeValue: The ID of NodeID needs to be a positive integer.
    """
    if NodeID < 1: raise NegativeValue()

    fix(NodeID, 1, 1, 1)


def Pin(NodeRID: int, NodeCID: int):
    """
    Function that constrains the translational DOF with a multi-point constraint.

    @param NodeRID (int): Node ID which will be retained by the multi-point constraint
    @param NodeCID (int): Node ID which will be constrained by the multi-point constraint

    @exception WrongArgument: The IDs passed needs to be different.
    @exception NegativeValue: The ID of NodeRID needs to be a positive integer.
    @exception NegativeValue: The ID of NodeCID needs to be a positive integer.
    """
    if NodeCID == NodeRID: raise WrongArgument()
    if NodeRID < 1: raise NegativeValue()
    if NodeCID < 1: raise NegativeValue()

	# Constrain the translational DOF with a multi-point constraint
	#   		retained constrained DOF_1 DOF_2
    equalDOF(NodeRID, NodeCID, 1, 2)


def RotationalSpring(ElementID: int, NodeRID: int, NodeCID: int, MatID: int, Rigid = False):
    """
    Function that defines a zero-length spring and constrains the translations DOFs of the spring. Can be used also to create rigid connections.

    @param ElementID (int): ID of the zerolength element that models the spring
    @param NodeRID (int): Node ID which will be retained by the multi-point constraint
    @param NodeCID (int): Node ID which will be constrained by the multi-point constraint
    @param MatID (int): ID of the material model chosen
    @param Rigid (bool, optional): Optional argument that transforms the joint in a completely rigid connection. Defaults to False.

    @exception NegativeValue: The ID of ElementID needs to be a positive integer.
    @exception NegativeValue: The ID of NodeCID needs to be a positive integer.
    @exception NegativeValue: The ID of NodeRID needs to be a positive integer.
    @exception WrongArgument: The IDs of the nodes passed needs to be different.
    @exception NegativeValue: The ID of MatID needs to be a positive integer.
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



	
