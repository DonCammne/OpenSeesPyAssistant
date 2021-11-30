from abc import ABC, abstractmethod
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

        delimiter = "##############################" # 30 times #
        for data_line in self.data:
            f.write('\n')
            if type(data_line) == np.ndarray:
                tmp_str = np.array_str(data_line, max_line_width = np.inf)
            else:
                tmp_str = str(data_line)
            f.write(tmp_str)
        f.write('\n')
        f.write(delimiter)

    @abstractmethod
    def ShowInfo(self):
        """Abstract function that show the data stored in the class in the command window.
        """
        pass
    