"""
.. module:: Autoencoder
   :synopsis: This module implements the Autoencoder class.
"""

import random as rd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Interpolation:
    """
    .. class:: Interpolation
        This class groups attributes and functions about the Autoencoder.
    """

    def __init__(self, alpha):
        self.alpha = alpha
        self.interpolated_ls = None
        self.decoded_interpolated_ls = None
        self.interpolated_predicted_img = None
        self.integrated_matrix = None

    def interpolate_predicted_img(self, predicted_hic, predicted_ngs):
        """
            Image space interpolation
        """
        self.interpolated_predicted_img = predicted_hic*(1-self.alpha) + predicted_ngs*self.alpha

    def interpolate_latent_spaces(self, hic_latent_spaces, ngs_latent_spaces):
        """
            Latent space interpolation
        """
        self.interpolated_ls = hic_latent_spaces*(1-self.alpha) + ngs_latent_spaces*self.alpha

    def set_decoded_latent_spaces(self, decoder):
        """
            Latent space interpolation
        """
        self.decoded_interpolated_ls = decoder.predict(self.interpolated_ls)

    def plot_random_interpolation(self, hic, predicted_hic, predicted_ngs,
                                  color_map, size_img, path):
        """
            TO DO
        """
        submatrix = rd.randint(0, hic.matrix.shape[0])

        fig, axes = plt.subplots(2, 3, figsize=(24, 11))
        fig.suptitle('Interpolation of the sub-matrix n°{}'.format(submatrix), fontsize=20)
        fig.subplots_adjust(left=0.03, right=0.98, wspace=0.3, hspace=0.4)
        for i in range(2):
            axes[i, 0].imshow(predicted_hic[submatrix].reshape(size_img, size_img),
                              cmap=color_map)
            axes[i, 0].set_title("Hi-C")
            axes[i, 0].axis('off')

            axes[i, 2].imshow(predicted_ngs[submatrix].reshape(size_img, size_img),
                              cmap=color_map)
            axes[i, 2].set_title("NGS")
            axes[i, 2].axis('off')

        axes[0, 1].imshow(self.interpolated_predicted_img[submatrix].reshape(size_img, size_img),
                          cmap=color_map)
        axes[0, 1].set_title("Interpolated predicted img")
        axes[0, 1].axis('off')

        axes[1, 1].imshow(self.decoded_interpolated_ls[submatrix].reshape(size_img, size_img),
                          cmap=color_map)
        axes[1, 1].set_title("Interpolated latent-spaces")
        axes[1, 1].axis('off')

        plt.savefig('{}/submatrix_{}.png'.format(path, submatrix))
        plt.close()

        plt.figure(figsize=(20, 8))
        plt.imshow(self.interpolated_ls[submatrix].reshape((hic.side//4) * 128, hic.side//4).T,
                   cmap=color_map)
        plt.title("Latent space")
        plt.axis('off')
        plt.savefig('{}/latentSpace_submatrix_{}.png'.format(path, submatrix))
        plt.close()

    def construct_integrated_matrix(self, hic):
        """
            Construction and set of the predicted Hi-C matrix from the predicted sub-matrices.
        """
        white = np.zeros(shape=(hic.side, hic.side, 1))
        line_limit = int(hic.matrix.shape[0] / hic.side)
        nb_sub_matrices = 1
        pred_sub_matrix_ind = 0

        for ind in range(hic.total_sub_matrices):
            if ind in hic.white_sub_matrices_ind:
                # The current sub-matrix is white
                sub_matrix = white
            else:
                sub_matrix = self.decoded_interpolated_ls[pred_sub_matrix_ind]
                pred_sub_matrix_ind += 1

            if nb_sub_matrices == line_limit:
                # The current line is concatenated with the previous lines
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                    integrated_matrix = np.concatenate((integrated_matrix, line), axis=0)
                except NameError:
                    integrated_matrix = line
                nb_sub_matrices = 1
                del line
            else:
                # A new sub-matrix is concatenated with the current line
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                except NameError:
                    line = sub_matrix
                nb_sub_matrices += 1

        self.integrated_matrix = integrated_matrix.reshape(integrated_matrix.shape[0],
                                                           integrated_matrix.shape[1])

    def plot_integrated_matrix(self, color_map, output_path):
        """
            The reconstructed and predicted Hi-C matrix is plotted in a file.

            Args:
                color_map(matplotlib.colors.ListedColormap): Color map for the plot
                output_path(str): Path of the output plot
        """
        fig = plt.figure(figsize=(12, 12))
        axes = plt.subplot(111, aspect='equal')
        img = axes.matshow(self.integrated_matrix, cmap=color_map)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="2%", pad=0.15)
        plt.colorbar(img, cax=cax)
        plt.subplots_adjust(left=0.07, bottom=0, right=0.95, top=0.91, wspace=0, hspace=0)
        axes.set_title('Predicted chr{} Hi-C matrix'.format('x'), fontsize=25)
        fig.savefig('{}/chr{}_integrated.png'.format(output_path, 'x'))
        plt.close()
