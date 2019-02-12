"""
.. module:: HiC
   :synopsis: This module implements the HiC class.
"""

# Third-party modules
import pandas as pd


class Hic:
    """
    .. class:: Hic

        This class groups informations about a HiC matrix.

    Attributes:
        matrix (Pandas DataFrame): HiC matrix.
    """

    def __init__(self, filename):
        self.matrix = pd.read_csv(filename, sep='\t', header=None)
        self.matrix.columns = ['chr', 'base_1', 'base_2', 'value']

    def get_value(self, base_1, base_2):
        """
            Search in a Pandas Dataframe (the HiC matrix) the unique row decribed
            by two given nucleotid bases and return the value corresponding.

            Args:
                base_1(int): Index of the first nucleotid base.
                base_2(int): Index of the second nucleotid base.
            Returns:
                value(float): The corresponding value if the row exists or 0 if not.
        """
        # The row could be find by two positions (base_1 base_2 or base_2 base_1)
        if not isinstance(base_1, int) or not isinstance(base_2, int):
            raise TypeError("Args must be integer.")

        filtered_matrix_1 = self.matrix[(self.matrix['base_1'] == base_1) &
                                        (self.matrix['base_2'] == base_2)]
        filtered_matrix_2 = self.matrix[(self.matrix['base_1'] == base_2) &
                                        (self.matrix['base_2'] == base_1)]

        # The row exists in one of the two filtered dataframes :
        if not filtered_matrix_1.empty:
            return float(filtered_matrix_1['value'])
        elif not filtered_matrix_2.empty:
            return float(filtered_matrix_2['value'])
        # The row does not exist :
        return 0
