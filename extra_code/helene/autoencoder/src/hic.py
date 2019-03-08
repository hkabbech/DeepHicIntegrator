"""
.. module:: Hic
   :synopsis: This module implements the matrix class.
"""

# Third-party modules
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable


REDS = cm.get_cmap('Reds', 300)
WHITE = np.array([1, 1, 1, 1])
NEW_CMP = ListedColormap(np.vstack((WHITE, REDS(np.linspace(0, 1, 300)))))

class Hic:
    """
    .. class:: Hic
        This class groups informations about a Hi-C matrix.

    Attributes:
        sparse (Pandas DataFrame): A data frame containing the Hi-C sparse matrix.
        resolution (int): The resolution of the matrix.
        size (int):
        matrix (numpy array):
    """

    def __init__(self, filename, resolution):
        self.df = pd.read_csv(filename, sep='\t', header=None)
        self.df.columns = ['chr', 'base_1', 'base_2', 'value']
        self.chr_name = ", ".join(self.df.chr.unique())
        self.resolution = resolution
        self.size = int(int(max(self. df['base_2'])) / self.resolution)
        self.matrix = None
        self.sub_matrices = None
        self.predicted_sub_matrices = None
        self.reconstructed_matrix = None

    def set_matrix(self, added_lines, deleted_lines):
        """
            Create a sparse matrix in a very fast way using scipy.sparse module
            and convert it into a numpy array.
        """
        # All the values are stored in a numpy array
        # The matrix will be symetric, so the values must be added twice
        values = self.df.value
        data = np.array(values.append(values))
        # base_1 and base_2 columns must be converted to index by dividing by the resolution number
        # This step is necesary for the creation of the sparse matrix with scipy.sparse
        base_1_index = ((self.df.base_1 / self.resolution)+added_lines).astype(int)
        base_2_index = ((self.df.base_2 / self.resolution)+added_lines).astype(int)
        row = np.array(base_1_index.append(base_2_index))
        col = np.array(base_2_index.append(base_1_index))
        # Creation of the sparse matrix and conversion into a numpy array
        matrix = coo_matrix((data, (row, col)),
                            shape=(self.size+added_lines+1, self.size+added_lines+1)
                           ).toarray()
        matrix = matrix[deleted_lines:, deleted_lines:]
        matrix = np.float32(matrix)
        matrix = np.log(matrix+1)
        # Rescaling of the values in range 0-1
        # (min-max scaling method)
        self.matrix = (matrix - matrix.min()) / (matrix.max() - matrix.min())

    def set_sub_matrices(self, nrow, ncol, nele=1):
        """
            Reshape of the sparse matrix in order to obtain N matrices of size nrow*ncol.

            Args:
                nrow (int): Number of rows for a sub-matrix
                ncol (int): Number of columns a sub-matrix
                nele (int): Number of element in a list (Default=1)
        """
        #self.sub_matrices = self.matrix.reshape(-1, nrow, ncol, nele)
        # Another less fancy method
        sub_matrices_list = []
        for i in range(0, self.matrix.shape[1], nrow):
            for j in range(0, self.matrix.shape[1], ncol):
                sub_matrix = self.matrix[i:i+nrow, j:j+ncol]
                sub_matrices_list.append(sub_matrix)
        sub_matrices = np.array(sub_matrices_list)
        self.sub_matrices = sub_matrices.reshape(-1, nrow, ncol, nele)

    def set_predicted_sub_matrices(self, predicted_sub_matrices):
        """
        TO DO doctrings
        """
        self.predicted_sub_matrices = predicted_sub_matrices

    def set_reconstructed_matrix(self, N):
        """
        TO DO doctrings
        """
        nb_sub_matrices = int(self.matrix.shape[0]/N)
        j = 1
        for _, sub_matrix in enumerate(self.predicted_sub_matrices):
            if j == nb_sub_matrices:
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                    reconstructed_matrix = np.concatenate((reconstructed_matrix, line), axis=0)
                except:
                    reconstructed_matrix = line
                j = 1
                del line
            else:
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                except:
                    line = sub_matrix
                j += 1
        self.reconstructed_matrix = reconstructed_matrix.reshape(reconstructed_matrix.shape[0],
                                                                 reconstructed_matrix.shape[1])

    def plot_matrix(self, chr_type="true"):
        """
        TO DO doctrings
        """
        if chr_type == "predicted":
            matrix = self.reconstructed_matrix
        else:
            matrix = self.matrix
        fig = plt.figure(figsize=(12, 12))
        ax = plt.subplot(111, aspect='equal')
        im = ax.matshow(matrix, cmap=NEW_CMP)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.15)
        plt.colorbar(im, cax=cax)
        plt.subplots_adjust(left=0.07, bottom=0, right=0.95, top=0.91, wspace=0, hspace=0)
        ax.set_title('{} {} Hi-C matrix'.format(chr_type, self.chr_name), fontsize=25)
        fig.savefig('{}_{}.png'.format(chr_type, self.chr_name))

    def plot_ten_sub_matrices(self, indices_list, chr_type="true"):
        """
        TO DO doctrings
        """
        if chr_type == "predicted":
            sub_matrices = self.predicted_sub_matrices
        else:
            sub_matrices = self.sub_matrices
        plt.figure(figsize=(22, 4))
        for i, sub_matrix_index in enumerate(indices_list):
            plt.subplot(2, 10, i+1)
            plt.imshow(sub_matrices[sub_matrix_index, ..., 0], cmap=NEW_CMP)
            plt.title("submatrix n°{}".format(sub_matrix_index))
            plt.savefig('{}_{}_10_submatrices.png'.format(chr_type, self.chr_name))
