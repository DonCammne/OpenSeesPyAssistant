# Module that handles the exceptions and errors of the library
#   Carmine Schipani, 2021

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

