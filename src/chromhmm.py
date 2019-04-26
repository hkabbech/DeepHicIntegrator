"""
.. module:: ChromHMM
   :synopsis: This module implements the Hic class.
"""

# Third-party modules
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix


class ChromHMM:
    """
    .. class:: Hic
        This class groups informations about a ChromHMM file.

    Attributes:
        TO DO
    """

    def __init__(self, chromhmm_filename, chrom):
        self.chromhmm_filename = chromhmm_filename
        self.chrom = chrom
        self.dataframe = None


    def set_df(self):
        """
            TO DO
            Preprocessing
        """
        chromhmm_df = pd.read_csv(self.chromhmm_filename, sep='\t',
                                  header=None, usecols=[0, 1, 2, 3])
        chromhmm_df.columns = ['chr', 'base_1', 'base_2', 'state']
        chromhmm_df = chromhmm_df.replace({"1_Active_Promoter": 0/10,
                                           "2_Weak_Promoter": 1/10,
                                           "3_Poised_Promoter": 2/10,
                                           "4_Strong_Enhancer": 3/10,
                                           "5_Strong_Enhancer": 3/10,
                                           "6_Weak_Enhancer": 4/10,
                                           "7_Weak_Enhancer": 4/10,
                                           "8_Insulator": 5/10,
                                           "9_Txn_Transition": 6/10,
                                           "10_Txn_Elongation": 6/10,
                                           "11_Weak_Txn": 7/10,
                                           "12_Repressed": 8/10,
                                           "13_Heterochrom/lo": 9/10,
                                           "14_Repetitive/CNV": 10/10,
                                           "15_Repetitive/CNV": 10/10
                                          })
        chromhmm_df = chromhmm_df[chromhmm_df.chr == "chr"+str(self.chrom)]
        self.dataframe = chromhmm_df


    def set_matrix(self):
        """
            TO DO
        """
        data = np.array(self.dataframe.state)
        row = np.array((self.dataframe.base_1 / 5000).astype(int))
        col = np.array((self.dataframe.base_2 / 5000).astype(int))

        coo_matrix((data, (row, col)))
