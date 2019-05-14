"""
.. module:: Interpolation
   :synopsis: This module implements the Interpolation class.
"""

# Third-party modules
import csv
import random as rd
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Interpolation:
    """
    .. class:: Interpolation
        This class groups attributes and functions about the Interpolation of matrices

    Attributes:
        alphas (list): ...
        TO DO
    """

    def __init__(self, alphas):
        self.alphas = alphas
        self.interpolated_submatrices = []
        self.integrated_matrix = []

    def construct_integrated_matrix(self, hic):
        """
            Construction and set of the predicted Hi-C matrix from the predicted sub-matrices.
        """
        white = np.zeros(shape=(hic.side, hic.side, 1))
        line_limit = int(hic.matrix.shape[0] / hic.side)

        integrated_matrix_list = []
        for _, matrix in enumerate(self.interpolated_submatrices):
            nb_sub_matrices = 1
            pred_sub_matrix_ind = 0
            for ind in range(hic.total_sub_matrices):
                # if ind in hic.white_sub_matrices_ind:
                    # The current sub-matrix is white
                    # sub_matrix = white
                # else:
                sub_matrix = matrix[pred_sub_matrix_ind]
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

            integrated_matrix_list.append(integrated_matrix.reshape(integrated_matrix.shape[0],
                                                                    integrated_matrix.shape[1]))
            del integrated_matrix

        self.integrated_matrix = integrated_matrix_list

    def write_predicted_sparse_matrix(self, hic, path, threshold=0.0001):
        """
            The reconstructed and predicted Hi-C matrix is saved in a sparse matrix file.

            Args:
                threshold(float): The predicted values under the threshold will be set to 0
                output_path(str): Path of the output plot
        """
        for i, integrated_matrix in enumerate(self.integrated_matrix):
            # Prediction under threshold value are set to 0
            # matrix[matrix < threshold] = 0
            # Creation of the sparse matrix
            sparse = coo_matrix(integrated_matrix)
            with open('{}/integrated_{}.bed'.format(path, i), 'w') as file:
                writer = csv.writer(file, delimiter='\t')
                writer.writerows(zip(['chr'+str(hic.chrom_num)]*len(sparse.row),
                                     sparse.row*hic.resolution,
                                     sparse.col*hic.resolution, sparse.data))

    def plot_integrated_matrix(self, hic, color_map, path):
        """
            The reconstructed and predicted Hi-C matrix is plotted in a file.

            Args:
                color_map(matplotlib.colors.ListedColormap): Color map for the plot
                output_path(str): Path of the output plot
        """
        for i, integrated_matrix in enumerate(self.integrated_matrix):
            fig = plt.figure(figsize=(12, 12))
            axes = plt.subplot(111, aspect='equal')
            img = axes.matshow(integrated_matrix, cmap=color_map)
            divider = make_axes_locatable(axes)
            cax = divider.append_axes("right", size="2%", pad=0.15)
            plt.colorbar(img, cax=cax)
            plt.subplots_adjust(left=0.07, bottom=0, right=0.95, top=0.91, wspace=0, hspace=0)
            axes.set_title('Integrated matrix (chr{}, resolution={}, alpha={})'\
                           .format(hic.chrom_num, hic.resolution, self.alphas[i]), fontsize=25)
            fig.savefig('{}/integrated_{}.png'.format(path, i))
            plt.close()

    def plot_interpolated_submatrices(self, hic, index_list, color_map, path):
        """
            The reconstructed and predicted Hi-C matrix is plotted in a file.

            Args:
                color_map(matplotlib.colors.ListedColormap): Color map for the plot
                output_path(str): Path of the output plot
        """
        for k, interpolated_submatrix in enumerate(self.interpolated_submatrices):
            fig, axes = plt.subplots(4, 10, figsize=(24, 11))
            fig.suptitle('Integrated sub-matrices (chr{}, resolution={}, alpha={})'\
                           .format(hic.chrom_num, hic.resolution, self.alphas[k]), fontsize=20)
            fig.subplots_adjust(left=0.03, right=0.98, wspace=0.3, hspace=0.4)
            i = 0
            for axe, index in zip(axes.flat, interpolated_submatrix[index_list, ..., 0]):
                axe.imshow(index, cmap=color_map)
                axe.set_title("submatrix n°{}".format(index_list[i]))
                i += 1
            plt.savefig('{}/interpolated_submatrices_{}.png'.format(path, k))
            plt.close()

    # def plot_random_interpolation(self, hic, color_map, path):
    #     """
    #         TO DO
    #     """
    #     submatrix = rd.randint(0, self.interpolated_pred_submatrices[0].shape[0])

    #     fig, axes = plt.subplots(2, 5, figsize=(24, 11))
    #     fig.suptitle('Interpolation of the sub-matrix n°{} (chr{}, resolution={})'
    #                  .format(submatrix, hic.chrom_num, hic.resolution), fontsize=20)
    #     fig.subplots_adjust(left=0.03, right=0.98, wspace=0.3, hspace=0.4)

    #     i = 0
    #     for axe, integrated_submatrices in zip(axes.flat, self.interpolated_pred_submatrices+\
    #                                            self.decoded_interpolated_ls):
    #         axe.imshow(integrated_submatrices[submatrix].reshape(hic.side, hic.side),
    #                    cmap=color_map)
    #         axe.axis('off')
    #         # axe.set_title('alpha={}'.format(self.alphas[i]))
    #         i += 1

    #     plt.savefig('{}/intergrated_submatrix_{}.png'.format(path, submatrix))
    #     plt.close()


class NormalInterpolation(Interpolation):
    """
    .. class:: InterpolationInLatentSpace
        This class groups attributes and functions about the Autoencoder.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def interpolate_predicted_img(self, hist_marks, predicted_hic):
        """
            Image space interpolation
        """
        for hist_mark in hist_marks.values():
            try:
                predicted_hm += (hist_mark.predicted_sub_matrices * (1/len(hist_marks)))
            except NameError:
                predicted_hm = (hist_mark.predicted_sub_matrices * (1/len(hist_marks)))
        for alpha in self.alphas:
            self.interpolated_submatrices.append(predicted_hic*(1-alpha) + predicted_hm*alpha)


class InterpolationInLatentSpace(Interpolation):
    """
    .. class:: InterpolationInLatentSpace
        This class groups attributes and functions about the Autoencoder.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.interpolated_ls = []

    def interpolate_latent_spaces(self, hist_marks, hic_latent_spaces):
        """
            Latent space interpolation
        """
        for hist_mark in hist_marks.values():
            try:
                hm_latent_spaces += (hist_mark.latent_spaces * (1/len(hist_marks)))
            except NameError:
                hm_latent_spaces = (hist_mark.latent_spaces * (1/len(hist_marks)))
        for alpha in self.alphas:
            self.interpolated_ls.append(hic_latent_spaces*(1-alpha) + hm_latent_spaces*alpha)

    def set_decoded_latent_spaces(self, decoder):
        """
            Latent space interpolation
        """
        for _, latent_spaces in enumerate(self.interpolated_ls):
            self.interpolated_submatrices.append(decoder.predict(latent_spaces))

    def plot_interpolate_ls(self, submatrix, hic, color_map, path):
        """
            TO DO
        """
        num_alpha = 2 # alpha=0.5
        submatrix = rd.randint(0, self.interpolated_ls[num_alpha].shape[0])
        plt.figure(figsize=(20, 8))
        plt.imshow(self.interpolated_ls[num_alpha][submatrix].reshape((hic.side//4) * 32,
                                                                       hic.side//4).T,
                   cmap=color_map)
        plt.title("Latent space (alpha={})".format(self.alphas[num_alpha]))
        plt.axis('off')
        plt.savefig('{}/latent_space_submatrix_{}.png'.format(path, submatrix))
        plt.close()
