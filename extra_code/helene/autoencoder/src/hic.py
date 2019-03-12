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

    def __init__(self, filename):
        self.df = pd.read_csv(filename, sep='\t', header=None)
        self.df.columns = ['chr', 'base_1', 'base_2', 'value']
        self.chr_name = ", ".join(self.df.chr.unique())
        self.matrix = None
        self.sub_matrices = None
        self.predicted_sub_matrices = None
        self.reconstructed_matrix = None

    def set_matrix(self, resolution, added_lines, deleted_lines):
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
        base_1_index = ((self.df.base_1 / resolution)+added_lines).astype(int)
        base_2_index = ((self.df.base_2 / resolution)+added_lines).astype(int)
        row = np.array(base_1_index.append(base_2_index))
        col = np.array(base_2_index.append(base_1_index))
        # Creation of the sparse matrix and conversion into a numpy array
        size = int(max(max(self.df['base_2']), max(self.df['base_1'])) / resolution)
        matrix = coo_matrix((data, (row, col)),
                            shape=(size+added_lines+1, size+added_lines+1)
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

    def set_reconstructed_matrix(self, side):
        """
        TO DO doctrings
        """
        nb_sub_matrices = int(self.matrix.shape[0] / side)
        j = 1
        for _, sub_matrix in enumerate(self.predicted_sub_matrices):
            if j == nb_sub_matrices:
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                    reconstructed_matrix = np.concatenate((reconstructed_matrix, line), axis=0)
                except NameError:
                    reconstructed_matrix = line
                j = 1
                del line
            else:
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                except NameError:
                    line = sub_matrix
                j += 1
        self.reconstructed_matrix = reconstructed_matrix.reshape(reconstructed_matrix.shape[0],
                                                                 reconstructed_matrix.shape[1])
    def save_reconstructed_matrix(self, output_path, resolution):
        """
        TO DO doctrings
        """
        sparse = coo_matrix(np.triu(self.reconstructed_matrix))
        with open(output_path+'/reconstructed_'+self.chr_name+'.txt', 'w') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerows(zip([self.chr_name]*len(sparse.row), sparse.row*resolution,
                                 sparse.col*resolution, sparse.data))

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

    def plot_ten_sub_matrices(self, chr_type, indices_list, output_path):
        """
        TO DO doctrings
        """
        if chr_type == "predicted":
            sub_matrices = self.predicted_sub_matrices
        else:
            sub_matrices = self.sub_matrices
        plt.figure(figsize=(18, 2))
        for i, sub_matrix_index in enumerate(indices_list):
            plt.subplot(1, 10, i+1)
            plt.imshow(sub_matrices[sub_matrix_index, ..., 0], cmap=NEW_CMP)
            plt.title("submatrix n°{}".format(sub_matrix_index))
        plt.subplots_adjust(left=0.03, right=0.98, wspace=0.3)
        plt.savefig('{}/10submatrices_{}_{}_.png'.format(output_path, self.chr_name, chr_type))


# N = 30
# chr20 = 180*180, chr10 = 84*84
#
# N = 60 
# chr20 =  90*90 , chr10 = 42*42
#
# N = 120
# chr20 =  45*45 , chr10 = 21*21
#
# N = 180
# chr20 =  30*30 , chr10 = 14*14
#
# N = 360
# chr20 =  15*15 , chr10 =  7*7 
#