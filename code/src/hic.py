"""
.. module:: Matrix
   :synopsis: This module implements the matrix class.
"""

# Third-party modules
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix


class Hic:
    """
    .. class:: Hic

        This class groups informations about a Hi-C matrix.

    Attributes:
        sparse (Pandas DataFrame): A data frame containing the Hi-C sparse matrix.
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
        self.size = int(self.end / self.resolution)

    def set_matrix(self):
        """
            Create a sparse matrix in a very fast way using scipy.sparse module
            and convert it into a numpy array.
        """
        # All the values are stored in a numpy array
        data = np.array(self.sparse.value)
        # base_1 and base_2 columns must be converted to index by dividing by the resolution number
        # This step is necesary for the creation of the sparse matrix with scipy.sparse
        row = np.array((self.sparse.base_1 / self.resolution).astype(int))
        col = np.array((self.sparse.base_2 / self.resolution).astype(int))
        # Creation of the sparse matrix and conversion into a numpy array
        self.matrix = coo_matrix((data, (row, col)), shape=(self.size+1, self.size+1)).toarray()

    def get_value(self, chr_num, base_1, base_2):
        """
            Search in a Pandas Dataframe (the Hi-C/NGS matrix) the unique row described
            by two given nucleotide bases and return the value corresponding.

            Args:
            	chr_num (str): Name pf the chromosome.
                base_1 (int): Index of the first nucleotide base.
                base_2 (int): Index of the second nucleotide base.
            Returns:
                float: The corresponding value if the row exists or 0 if not.
        """
        # The row could be find by two combinations (base_1 base_2 or base_2 base_1)
        if not isinstance(base_1, int) or not isinstance(base_2, int):
            raise TypeError("Args must be integer.")

        filtered_matrix_1 = self.sparse[(self.sparse['chr'] == chr_num) &
                                        (self.sparse['base_1'] == base_1) &
                                        (self.sparse['base_2'] == base_2)]
        filtered_matrix_2 = self.sparse[(self.sparse['chr'] == chr_num) &
                                        (self.sparse['base_1'] == base_2) &
                                        (self.sparse['base_2'] == base_1)]

        # The row do exist in one of the two filtered data frames :
        if not filtered_matrix_1.empty:
            return float(filtered_matrix_1['value'])
        elif not filtered_matrix_2.empty:
            return float(filtered_matrix_2['value'])
        # The row does not exist :
        return 0
