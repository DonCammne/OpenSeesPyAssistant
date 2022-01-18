"""
Module with the parent abstract class DataManagement. \n
Carmine Schipani, 2021
"""

from abc import ABC, abstractmethod
from OpenSeesPyAssistant.ErrorHandling import *
import numpy as np


class DataManagement(ABC):
    """
    Abstract parent class for data management.
    Using the associated MATLAB class \n
    LOAD_CLASS.m \n
    for the postprocessing in MATLAB, allowing for simpler and more reliable data management because the parameters
    from the OpenSeesPy analysis are imported automatically. 
    """

    def SaveData(self, f):
        """
        Function that lists in the command window and saves in a opened file text "f" the data from the "self" class that calls it. 
        Example: call this function after this line: \n 
        with open(FileName, 'w') as f:

        @param f (io.TextIOWrapper): Opened file to write into

        @exception WrongDimension: The number of lists in the list self.data needs to be 2
        """
        if len(self.data[0]) != 2: raise WrongDimension() 
        
        delimiter = "##############################" # 30 times #
        col_delimiter = "\t" # tab
        for data_line in self.data:
            f.write('\n')
            for col in data_line:
                if type(col) == np.ndarray:
                    tmp_str = np.array_str(col, max_line_width = np.inf)
                else:
                    tmp_str = str(col)
                f.write(tmp_str)
                f.write(col_delimiter)
        f.write('\n')
        f.write('NEW INFO SECTION DELIMITER \t')
        f.write(delimiter)

    @abstractmethod
    def ShowInfo(self):
        """
        Abstract method that shows the data stored in the class in the command window.
        In some cases, it's possible to plot some information (for example the curve of the material model).
        """
        pass

    @abstractmethod
    def ReInit(self):
        """
        Abstract method that computes the value of the parameters with respect of the arguments. \n
        Use after changing the value of argument inside the class (to update the values accordingly). \n
        This function can be very useful in combination with the function "deepcopy()" from the module "copy". \n
        Be careful that the parameter self.Initialized is also copied, thus it is safer to copy the class before the method that calls the actual OpenSees commands (and initialise the object).
        """
        pass

    @abstractmethod
    def UpdateStoredData(self):
        """
        Abstract method used to define and update the self.data member variable. \n
        This member variable (self.data) is a list of lists with 2 entries (info_name and info_value)
            and for each list is stored a different member variable of the class. \n
        Useful to debug the model, export data, copy object.
        """
        pass