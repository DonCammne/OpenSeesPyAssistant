# Module for the parent class that manage data
#         Carmine Schipani, 2021

from abc import ABC, abstractmethod
from OpenSeesPyHelper.ErrorHandling import WrongDimension
import numpy as np

class DataManagement(ABC):
    """Abstract parent class for data management. Using the associated MATLAB 
    LOAD_CLASS.m
    class, the postprocessing in MATLAB is simplified and reliable because the parameters from the OpenSeesPy analysis are imported automatically. 
    """

    def SaveData(self, f):
        """Function that list in the command window and saves in a opened file text "f" the data from the "self" class that calls it. 
        Normally, call this function after this line: \n 
        with open(FileName, 'w') as f:

        Args:
            f (io.TextIOWrapper): Opened file to write into
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
        """Abstract function that show the data stored in the class in the command window.
        """
        pass

    @abstractmethod
    def ReInit(self):
        """Abstract function that computes the value of the parameters that are computed with respect of the arguments.
        Use after changing the value of argument inside the class (to update the values accordingly). 
        This function can be very useful in combination with the function "copy()" from the module "copy".
        """
        pass

    @abstractmethod
    def UpdateStoredData(self):
        """Abstract function used to define a self.data member.
        This member is a list of lists with 2 columns (info_name and info_value) and at each row is stored a different member of the class.
        """
        pass