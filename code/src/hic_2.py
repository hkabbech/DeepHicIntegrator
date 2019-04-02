"""
.. module:: Hic
   :synopsis: This module implements the matrix class.
"""

# Third-party modules
import csv
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math as m


REDS = cm.get_cmap('Reds', 300)
WHITE = np.array([1, 1, 1, 1])
NEW_CMP = ListedColormap(np.vstack((WHITE, REDS(np.linspace(0, 1, 300)))))

class Hic:
    """
    .. class:: Hic
        This class groups informations about a Hi-C matrix.

    Attributes:
        sparse (Pandas DataFrame): A data frame containing the Hi-C sparse matrix.
        matrix (numpy array):
    """

    def __init__(self, cooler, chrom):
        self.cooler = cooler
        self.resolution = self.cooler.info['bin-size']
        self.chrom = chrom
        self.matrix = None
        self.sub_matrices = None
        self.predicted_sub_matrices = None
        self.reconstructed_matrix = None

    def calculate_cum_length(self):
        chroms = self.cooler.chroms()[:]
        cum_length = []
        for i, size in enumerate(list(chroms.length)):
            if i == 0:
                cum_length.append(size)
            else:
                cum_length.append(size + cum_length[i-1])
        chroms['cum_length'] = cum_length
        return(chroms)

    def set_matrix(self):
        """
            Create a sparse matrix in a very fast way using scipy.sparse module
            and convert it into a numpy array.
        """
        chroms = self.calculate_cum_length()
        if self.chrom == 1:
            bin_1 = 0
        else:
            bin_1 = m.ceil(chroms[self.chrom-2:self.chrom-1]['cum_length']/self.resolution)
        bin_2 = m.ceil(chroms[self.chrom-1:self.chrom]['cum_length']/self.resolution)
        # Creation of the sparse matrix and conversion into a numpy array
        matrix = self.cooler.matrix(balance=False, sparse=True)[bin_1:bin_2, bin_1:bin_2].toarray()
        ## 0s are interpreted as missing values by Keras
        # matrix = matrix + 1
        # The matrix is symetric then we keep only the upper triangular matrix
        matrix = np.triu(matrix)
        # Conversion into float32
        matrix = np.float32(matrix)
        # Log scale to visualize better the matrix in plot
        matrix = np.log2(matrix+1)
        # Rescaling of the values in range 0-1
        # (min-max scaling method)
        matrix = (matrix - matrix.min()) / (matrix.max() - matrix.min())
        # # Lines are removing or added (0s)
        # if lines < 0:
        #     matrix = matrix[abs(lines):, abs(lines):]
        # elif lines > 0
        #     matrix = np.insert(matrix, 0, 0, axis = 0)
        #     matrix = np.insert(matrix, 0, 0, axis = 1)
        self.matrix = matrix

    def set_sub_matrices(self, n, channel=1):
        sub_matrices_list = []
        for i in range(0, self.matrix.shape[1], n):
            for j in range(0, self.matrix.shape[1], n):
                sub_matrix = self.matrix[i:i+n, j:j+n]
                if sub_matrix.shape != (n, n):
                    break
                if sub_matrix.sum() != 0:
                    sub_matrices_list.append(sub_matrix)
        sub_matrices = np.array(sub_matrices_list)
        sub_matrices = sub_matrices.reshape(-1, n, n, channel)


    def plot_matrix(self, chr_type, output_path):
        """
        TO DO doctrings
        """
        if chr_type == "predicted":
            matrix = self.reconstructed_matrix
        else:
            matrix = self.matrix
        fig = plt.figure(figsize=(12, 12))
        axes = plt.subplot(111, aspect='equal')
        img = axes.matshow(matrix, cmap=NEW_CMP)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="2%", pad=0.15)
        plt.colorbar(img, cax=cax)
        plt.subplots_adjust(left=0.07, bottom=0, right=0.95, top=0.91, wspace=0, hspace=0)
        axes.set_title('{} {} Hi-C matrix'.format(chr_type, self.chr_name), fontsize=25)
        fig.savefig('{}/{}_{}.png'.format(output_path, self.chr_name, chr_type))


