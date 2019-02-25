"""
.. module:: Matrix
   :synopsis: This module implements the matrix class.
"""

# Third-party modules
import pandas as pd
import numpy as np


class Hic:
    """
    .. class:: Hic

        This class groups informations about a Hi-C matrix.

    Attributes:
        sparse (Pandas DataFrame): A dataframe containing the Hi-C sparse matrix.
        resolution (int): The resolution of the matrix.
        start (int):
        end (int):
        size (int):
        matrix (numpy array):
    """

    def __init__(self, filename, resolution):
        self.sparse = pd.read_csv(filename, sep='\t', header=None)
        self.sparse.columns = ['chr', 'base_1', 'base_2', 'value']
        self.resolution = resolution
        self.start = int(self.sparse['base_1'].iloc[0])
        self.end = int(self.sparse['base_2'].iloc[-1])
        self.size = int((self.end - self.start) / self.resolution)
        self.matrix = np.empty((self.size, self.size), dtype=object)

    def get_value(self, chr_num, base_1, base_2):
        """
            Search in a Pandas Dataframe (the Hi-C/NGS matrix) the unique row decribed
            by two given nucleotid bases and return the value corresponding.

            Args:
                base_1 (int): Index of the first nucleotid base.
                base_2 (int): Index of the second nucleotid base.
            Returns:
                float: The corresponding value if the row exists or 0 if not.
        """
        # The row could be find by two combinaisons (base_1 base_2 or base_2 base_1)
        if not isinstance(base_1, int) or not isinstance(base_2, int):
            raise TypeError("Args must be integer.")

        filtered_matrix_1 = self.sparse[(self.sparse['chr'] == chr_num) &
                                        (self.sparse['base_1'] == base_1) &
                                        (self.sparse['base_2'] == base_2)]
        filtered_matrix_2 = self.sparse[(self.sparse['chr'] == chr_num) &
                                        (self.sparse['base_1'] == base_2) &
                                        (self.sparse['base_2'] == base_1)]

        # The row do exist in one of the two filtered dataframes :
        if not filtered_matrix_1.empty:
            return float(filtered_matrix_1['value'])
        elif not filtered_matrix_2.empty:
            return float(filtered_matrix_2['value'])
        # The row does not exist :
        return 0

    def set_matrix(self, chr_num):
        """
            Function test
        """
        for index_i, base_i in enumerate(range(self.start, self.end, self.resolution)):
            for index_j, base_j in enumerate(range(index_i, self.end, self.resolution)):
                self.matrix[index_i][index_j] = self.sparse.get_value(chr_num, base_i, base_j)
