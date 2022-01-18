"""Module dedicated to the error handling. \n
Carmine Schipani, 2021
"""


class ZeroDivision(Exception):
    """Exception class for the "zero division" error.
    """
    pass

class WrongArgument(Exception):
    """Exception class for the "input of a wrong argument" error.
    """
    pass

class NegativeValue(Exception):
    """Exception class for the "negative value (argument or result)" error.
    """
    pass

class PositiveValue(Exception):
    """Exception class for the "positive value (argument or result)" error.
    """
    pass

class WrongDimension(Exception):
    """Exception class for the "wrong array dimensions" error.
    """
    pass

class InconsistentGeometry(Exception):
    """Exception class for the "inconsistent geometry" error.
    """
    pass

class MemberFailure(Exception):
    """Exception class for the "member failure" error.
    """
    pass

class WrongNodeIDConvention(Exception):
    """Exception class for the "wrong node ID convention definition" error.
    """
    def __init__(self, node):
        self.node = node

class NoApplicability(Exception):
    """Exception class for the "no applicability of formula of theory" error.
    """
    pass

class ZeroLength(Exception):
    """Exception class for the "zero length element (non intentional)" error.
    """
    def __init__(self, element):
        self.element = element