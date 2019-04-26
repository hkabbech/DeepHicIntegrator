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
        plt.imshow(self.interpolated_ls[i].reshape(8 * 128, 8).T, cmap=color_map)
        plt.title("Latent space")
        plt.axis('off')
        plt.savefig('{}/latentSpace_submatrix_{}.png'.format(path, submatrix))
        plt.close()
