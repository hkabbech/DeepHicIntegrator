"""
.. module:: PredictHic
   :synopsis: This module implements the PredictHic class.
"""

# Third-party modules
import csv
import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Parent Hic class
from src.hic import Hic


class PredictHic(Hic):
    """
    .. class:: PredictHic
        This class herits attributes and methods from Hic class and
        groups informations about a predicted Hi-C matrix.

    Attributes:
        predicted_sub_matrices (numpy array): Sub-matrices predicted by a deep-learning model
        predicted_matrix (numpy array): The reconstructed matrix from the predicted sub-matrices
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.predicted_sub_matrices = None
        self.predicted_matrix = None

    def set_predicted_sub_matrices(self, predicted_sub_matrices):
        """
            Set the sub-matrices predicted by a deep-learning model.

            Args:
                predicted_sub_matrices(Numpy array): The predicted sub-matrices
        """
        self.predicted_sub_matrices = predicted_sub_matrices

    def construct_predicted_matrix(self):
        """
            Construction and set of the predicted Hi-C matrix from the predicted sub-matrices.
        """
        white_sub_matrix = np.zeros(shape=(self.side, self.side, 1))
        nb_sub_matrices = int(self.matrix.shape[0] / self.side)
        nb_white = 0
        j = 1
        for _, sub_matrix in enumerate(self.predicted_sub_matrices):
            if j == nb_sub_matrices - 1:
                # The current line is concatenated with the previous lines
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                    predicted_matrix = np.concatenate((predicted_matrix, line), axis=0)
                except NameError:
                    predicted_matrix = line
                nb_white += 1
                j = 1
                del line
            else:
                # A new sub-matrix is concatenated with the current line
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                except NameError:
                    if nb_white == 0:
                        line = sub_matrix
                    else:
                        # The matrix is upper triangular, we have to add white sub-matrices
                        line = white_sub_matrix
                        for _ in range(1, nb_white):
                            line = np.concatenate((line, white_sub_matrix), axis=1)
                            j += 1
                        line = np.concatenate((line, sub_matrix), axis=1)
                        j += 1
                j += 1
        self.predicted_matrix = predicted_matrix.reshape(predicted_matrix.shape[0],
                                                         predicted_matrix.shape[1])

    def plot_predicted_sub_matrices(self, color_map, output_path, index_list):
        """
            40 predicted random sub-matrices are plotted in a file.

            Args:
                color_map(matplotlib.colors.ListedColormap): Color map for the plot
                output_path(str): Path of the output plot
                index_list(list): List of the 40 sub-matrix indexes to plot
        """
        fig, axes = plt.subplots(4, 10, figsize=(24, 11))
        fig.suptitle('Predicted chr{} Hi-C sub-matrices'.format(self.chrom), fontsize=20)
        fig.subplots_adjust(left=0.03, right=0.98, wspace=0.3, hspace = 0.4)
        i = 0
        for ax, index in zip(axes.flat, self.predicted_sub_matrices[index_list, ..., 0]):
            ax.imshow(index, cmap=color_map)
            ax.set_title("submatrix n°{}".format(index_list[i]))
            i += 1
        plt.savefig('{}/submatrices_chr{}_predicted.png'.format(output_path, self.chrom))

    def plot_predicted_matrix(self, color_map, output_path):
        """
            The reconstructed and predicted Hi-C matrix is plotted in a file.

            Args:
                color_map(matplotlib.colors.ListedColormap): Color map for the plot
                output_path(str): Path of the output plot
        """
        fig = plt.figure(figsize=(12, 12))
        axes = plt.subplot(111, aspect='equal')
        img = axes.matshow(self.predicted_matrix, cmap=color_map)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="2%", pad=0.15)
        plt.colorbar(img, cax=cax)
        plt.subplots_adjust(left=0.07, bottom=0, right=0.95, top=0.91, wspace=0, hspace=0)
        axes.set_title('Predicted chr{} Hi-C matrix'.format(self.chrom), fontsize=25)
        fig.savefig('{}/chr{}_predicted.png'.format(output_path, self.chrom))

    def write_predicted_matrix(self, threshold, output_path):
        """
            The reconstructed and predicted Hi-C matrix is saved in a sparse matrix file.

            Args:
                threshold(float): The predicted values under the threshold will be set to 0
                output_path(str): Path of the output plot
        """
        # Prediction under threshold value are set to 0
        self.predicted_matrix[self.predicted_matrix < threshold] = 0
        # Creation of the sparse matrix
        sparse = coo_matrix(self.predicted_matrix)
        with open(output_path+'/predicted_chr'+str(self.chrom)+'.txt', 'w') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerows(zip(['chr'+str(self.chrom)]*len(sparse.row),
                                 sparse.row*self.resolution,
                                 sparse.col*self.resolution, sparse.data))
