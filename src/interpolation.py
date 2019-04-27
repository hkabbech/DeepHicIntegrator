"""
.. module:: Autoencoder
   :synopsis: This module implements the Autoencoder class.
"""

import random as rd
import matplotlib.pyplot as plt

class Interpolation:
    """
    .. class:: Interpolation
        This class groups attributes and functions about the Autoencoder.
    """

    def __init__(self, alpha, hic_latent_spaces, predicted_hic,
                 chromhmm_latent_spaces, predicted_chromhmm, decoder):
        self.alpha = alpha
        self.hic_ls = hic_latent_spaces
        self.predicted_hic = predicted_hic
        self.chromhmm_ls = chromhmm_latent_spaces
        self.predicted_chromhmm = predicted_chromhmm
        self.decoder = decoder
        self.interpolated_ls = None
        self.decoded_interpolated_ls = None
        self.interpolated_predicted_img = None
        self.reconstructed_matrix = None

    def interpolate_predicted_img(self):
        """
            Image space interpolation
        """
        self.interpolated_predicted_img = self.predicted_hic*(1-self.alpha) +\
                                          self.predicted_chromhmm*self.alpha

    def interpolate_latent_spaces(self):
        """
            Latent space interpolation
        """
        self.interpolated_ls = self.hic_ls*(1-self.alpha) + self.chromhmm_ls*self.alpha

    def set_decoded_latent_spaces(self):
        """
            Latent space interpolation
        """
        self.decoded_interpolated_ls = self.decoder.predict(self.interpolated_ls)

    def plot_random_interpolation(self, color_map, size_img, path):
        """
            TO DO
        """
        submatrix = rd.randint(0, self.predicted_hic.shape[0])

        fig, axes = plt.subplots(2, 3, figsize=(24, 11))
        fig.suptitle('Interpolation of the sub-matrix n°{}'.format(submatrix), fontsize=20)
        # fig.subplots_adjust(left=0.03, right=0.98, wspace=0.3, hspace=0.4)
        for i in range(2):
            axes[i, 0].imshow(self.predicted_hic[submatrix].reshape(size_img, size_img),
                              cmap=color_map)
            axes[i, 0].set_title("Hi-C")
            axes[i, 0].axis('off')

            axes[i, 2].imshow(self.predicted_chromhmm[submatrix].reshape(size_img, size_img),
                              cmap=color_map)
            axes[i, 2].set_title("ChromHMM")
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
        plt.imshow(self.interpolated_ls[i].reshape(15 * 128, 15).T, cmap=color_map)
        plt.title("Latent space")
        plt.axis('off')
        plt.savefig('{}/latentSpace_submatrix_{}.png'.format(path, submatrix))
        plt.close()

    def construct_predicted_matrix(self, predicted_hic):
        """
            Construction and set of the predicted Hi-C matrix from the predicted sub-matrices.
        """
        white = np.zeros(shape=(predicted_hic.side, predicted_hic.side, 1))
        line_limit = int(predicted_hic.matrix.shape[0] / predicted_hic.side)
        nb_sub_matrices = 1
        pred_sub_matrix_ind = 0

        for ind in range(predicted_hic.total_sub_matrices):
            if ind in predicted_hic.white_sub_matrices_ind:
                # The current sub-matrix is white
                sub_matrix = white
            else:
                sub_matrix = self.decoded_interpolated_ls[pred_sub_matrix_ind]
                pred_sub_matrix_ind += 1

            if nb_sub_matrices == line_limit:
                # The current line is concatenated with the previous lines
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                    predicted_matrix = np.concatenate((predicted_matrix, line), axis=0)
                except NameError:
                    predicted_matrix = line
                nb_sub_matrices = 1
                del line
            else:
                # A new sub-matrix is concatenated with the current line
                try:
                    line = np.concatenate((line, sub_matrix), axis=1)
                except NameError:
                    line = sub_matrix
                nb_sub_matrices += 1

        self.reconstructed_matrix = predicted_matrix.reshape(predicted_matrix.shape[0],
                                                             predicted_matrix.shape[1])

    # def plot_sub_matrices(self, color_map, output_path, index_list):
    #     """
    #         40 random sub-matrices are plotted in a file.

    #         Args:
    #             color_map(matplotlib.colors.ListedColormap): Color map for the plot
    #             output_path(str): Path of the output plot
    #             index_list(list): List of the 40 sub-matrix indexes to plot
    #     """
    #     fig, axes = plt.subplots(4, 10, figsize=(24, 11))
    #     fig.suptitle('True chr{} Hi-C sub-matrices'.format(self.chrom), fontsize=20)
    #     fig.subplots_adjust(left=0.03, right=0.98, wspace=0.3, hspace=0.4)
    #     i = 0
    #     for axe, index in zip(axes.flat, self.sub_matrices[index_list, ..., 0]):
    #         axe.imshow(index, cmap=color_map)
    #         axe.set_title("submatrix n°{}".format(index_list[i]))
    #         i += 1
    #     plt.savefig('{}/interpolated_pred_submatrices_chr{}_true.png'.format(output_path, self.chrom))
    #     plt.close()

    # def plot_matrix(self, color_map, output_path):
    #     """
    #         The Hi-C matrix is plotted in a file.

    #         Args:
    #             color_map(matplotlib.colors.ListedColormap): Color map for the plot
    #             output_path(str): Path of the output plot
    #     """
    #     fig = plt.figure(figsize=(12, 12))
    #     axes = plt.subplot(111, aspect='equal')
    #     img = axes.matshow(self.matrix, cmap=color_map)
    #     divider = make_axes_locatable(axes)
    #     cax = divider.append_axes("right", size="2%", pad=0.15)
    #     plt.colorbar(img, cax=cax)
    #     plt.subplots_adjust(left=0.07, bottom=0, right=0.95, top=0.91, wspace=0, hspace=0)
    #     axes.set_title('True chr{} Hi-C matrix'.format(self.chrom), fontsize=25)
    #     fig.savefig('{}/chr{}_true.png'.format(output_path, self.chrom))
    #     plt.close()
