"""
.. module:: Interpolation
   :synopsis: This module implements the Interpolation class.
"""

# Third-party modules
import csv
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Interpolation:
    """
    .. class:: Interpolation
        This class groups attributes and functions which aim to construct, write in a sparse matrix
        and plot two or several interpolated matrices.

    Attributes:
        alphas (list): List of float values to use for the interpolation (alpha parameter)
        interpolated_submatrices (list): List of all the interpolated sub-matrices. Each item in
                                         the list contains an interpolation with a different alpha.
        integrated_matrix (list): List of all the integrated (interpolated) reconstructed matrices.
                                  Each item in the list contains an interpolation with a different
                                  alpha.
    """

    def __init__(self, alphas):
        self.alphas = alphas
        self.interpolated_submatrices = []
        self.integrated_matrix = []

    def construct_integrated_matrix(self, hic):
        """
            Construction of the whole integrated matrices from the interpolated sub-matrices.

            Args:
                hic(Hic(Matrix) object): Hi-C matrix
        """
        line_limit = int(hic.matrix.shape[0] / hic.side)

        integrated_matrix_list = []
        for _, matrix in enumerate(self.interpolated_submatrices):
            nb_sub_matrices = 1
            for ind in range(hic.total_sub_matrices):
                sub_matrix = matrix[ind]
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
            The integrated matrices are saved in sparse matrix files for each alpha value.

            Args:
                hic(Hic(Matrix) object): Hi-C matrix
                path(str): Path of the output
                threshold(float): The values under the threshold will be set to 0
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
            The integrated matrices are plotted for each alpha value.

            Args:
                hic(Hic(Matrix) object): Hi-C matrix
                color_map(matplotlib.colors.ListedColormap): Color map
                path(str): Path of the output plot
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
            axes.axis('off')
            fig.savefig('{}/integrated_{}.pdf'.format(path, i))
            plt.close()

    def plot_interpolated_submatrices(self, hic, index_list, color_map, path):
        """
            40 random integrated sub-matrices are plotted for each alpha value.

            Args:
                hic(Hic(Matrix) object): Hi-C matrix
                index_list(list): List of the 40 sub-matrix indexes to plot
                color_map(matplotlib.colors.ListedColormap): Color map
                path(str): Path of the output plot
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
                axe.axis('off')
                i += 1
            plt.savefig('{}/interpolated_submatrices_{}.pdf'.format(path, k))
            plt.close()


class NormalInterpolation(Interpolation):
    """
    .. class:: InterpolationInLatentSpace
        This class inherits the Interpolation class and interpolate sub-matrices in the pixel space
        (= without the use of encoder and decoder).


    Attributes:
        alphas (list): List of float values to use for the interpolation (alpha parameter)
        interpolated_submatrices (list): List of all the interpolated sub-matrices. Each item in
                                         the list contains an interpolation with a different alpha.
        integrated_matrix (list): List of all the integrated (interpolated) reconstructed matrices.
                                  Each item in the list contains an interpolation with a different
                                  alpha.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def interpolate_predicted_img(self, hist_marks, predicted_hic):
        """
            Double linear interpolation of the predicted sub-matrices of the Hi-C and histone marks.

            Args:
                hist_marks(dict): Dictionary containing all histone mark HistoneMark objects.
                predicted_hic(numpy array): Predicted sub-matrices of the Hi-C
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
        This class inherits the Interpolation class and interpolate sub-matrices in the latent space
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.interpolated_ls = []

    def interpolate_latent_spaces(self, hist_marks, hic_latent_spaces):
        """
            Double linear interpolation of the latent spaces of the Hi-C and histone marks.

            Args:
                hist_marks(dict): Dictionary containing all histone mark HistoneMark objects.
                predicted_hic(numpy array): Predicted sub-matrices of the Hi-C
        """
        for hist_mark in hist_marks.values():
            try:
                hm_latent_spaces += (hist_mark.latent_spaces * (1/len(hist_marks)))
            except NameError:
                hm_latent_spaces = (hist_mark.latent_spaces * (1/len(hist_marks)))
        for alpha in self.alphas:
            self.interpolated_ls.append(hic_latent_spaces*(1-alpha) + hm_latent_spaces*alpha)

    def set_decoded_latent_spaces(self, decoder, side):
        """
            The interpolated latent spaces are decoded.

            Args:
                decoder(keras model object): Hi-C matrix
                side(int): Square side
        """
        for _, latent_spaces in enumerate(self.interpolated_ls):
            decoded_ls = np.array(np.zeros(shape=(latent_spaces.shape[0], side, side, 1)))
            for i, ls in enumerate(latent_spaces):
                if ls.sum() == 0:
                    continue
                decoded_ls[i] = decoder.predict(ls.reshape(1, ls.shape[0], ls.shape[1],
                                                           ls.shape[2]))
            self.interpolated_submatrices.append(decoded_ls)
