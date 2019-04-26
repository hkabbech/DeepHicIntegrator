"""
.. module:: Autoencoder
   :synopsis: This module implements the Autoencoder class.
"""

# Third-party modules
from contextlib import redirect_stdout
from sklearn.model_selection import train_test_split
from keras.layers import Input, Conv2D, MaxPooling2D, UpSampling2D
from keras.models import Model
from keras.optimizers import RMSprop
import matplotlib.pyplot as plt


class Autoencoder:
    """
    .. class:: Autoencoder
        This class groups attributes and functions about the Autoencoder.

    Attributes:
        chr_train(Hic object): Instance of the HicPredict class - Chromosome for training the model
        model(Model object): Instance of the Model keras class - Model of the Autoencoder
        trained_model(Model object): Instance of the Model keras class - The trained Model

    """

    def __init__(self, chr_train, img_size, epochs, batch_size):
        self.chr_train = chr_train
        self.img_size = img_size
        self.epochs = epochs
        self.batch_size = batch_size
        self.encoder = None
        self.decoder = None
        self.cae = None
        self.trained_cae = None

    def set_models(self):
        """
            Layers of the Autoencoder network.

            Args:
                input_img(Input class): Input of the network

            Returns:
                decoded(Conv2D keras layer): The last decoded output layer of the Autoencoder
        """

        input_img = Input(shape=(self.img_size, self.img_size, 1))

        tmp = Conv2D(32, (3, 3), activation='relu', padding='same')(input_img)
        tmp = MaxPooling2D((2, 2), padding='same')(tmp)
        tmp = Conv2D(64, (3, 3), activation='relu', padding='same')(tmp)
        tmp = MaxPooling2D((2, 2), padding='same')(tmp)
        tmp = Conv2D(128, (3, 3), activation='relu', padding='same')(tmp)
        encoded = MaxPooling2D((2, 2), padding='same')(tmp)

        conv_decoded_1 = Conv2D(128, (3, 3), activation='relu', padding='same')
        up_decoded_1 = UpSampling2D((2, 2))
        conv_decoded_2 = Conv2D(64, (3, 3), activation='relu', padding='same')
        up_decoded_2 = UpSampling2D((2, 2))
        conv_decoded_3 = Conv2D(32, (3, 3), activation='relu')
        up_decoded_3 = UpSampling2D((2, 2))
        conv_decoded_4 = Conv2D(1, (3, 3), activation='sigmoid', padding='same')

        self.encoder = Model(input_img, encoded)

        tmp = conv_decoded_1(encoded)
        tmp = up_decoded_1(tmp)
        tmp = conv_decoded_2(tmp)
        tmp = up_decoded_2(tmp)
        tmp = conv_decoded_3(tmp)
        tmp = up_decoded_3(tmp)
        decoded_autoencoder = conv_decoded_4(tmp)

        self.cae = Model(input_img, decoded_autoencoder)

        latent_space_input = Input(shape=(8, 8, 128))
        tmp = conv_decoded_1(latent_space_input)
        tmp = up_decoded_1(tmp)
        tmp = conv_decoded_2(tmp)
        tmp = up_decoded_2(tmp)
        tmp = conv_decoded_3(tmp)
        tmp = up_decoded_3(tmp)
        decoded = conv_decoded_4(tmp)

        self.decoder = Model(latent_space_input, decoded)


    def compile(self, loss_function='mean_squared_error', metrics_list=['accuracy']):
        """
            Compilation of the Autoencoder model.

            Args:
                loss_function(str): The loss function to use
                matrics_list(list): The metrics to generate during compilation
        """
        self.cae.compile(loss=loss_function, optimizer=RMSprop(), metrics=metrics_list)
        self.cae.summary()

    def train(self):
        """
            Training of the Autoencoder model and set of the trained_model attribute.
        """
        train_x, valid_x, train_ground, valid_ground = train_test_split(self.chr_train.sub_matrices,
                                                                        self.chr_train.sub_matrices,
                                                                        test_size=0.2,
                                                                        random_state=13)

        self.trained_cae = self.cae.fit(train_x, train_ground, verbose=1,
                                        batch_size=self.batch_size, epochs=self.epochs,
                                        validation_data=(valid_x, valid_ground))

    def plot_loss_curve(self, path):
        """
            Plot of the loss curve in a file.

            Args:
                epochs(int): Number of epochs
                path(str): Path of the output plot
        """
        plt.plot(range(self.epochs), self.trained_cae.history['loss'], 'bo', label='Training')
        plt.plot(range(self.epochs), self.trained_cae.history['val_loss'], 'b', label='Validation')
        plt.title('Training and validation loss')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.savefig(path+'training_validation_loss.png')
        plt.close()

    def plot_accuracy_curve(self, path):
        """
            Plot of the accuracy curve in a file.

            Args:
                epochs(int): Number of epochs
                path(str): Path of the output plot
        """
        plt.plot(range(self.epochs), self.trained_cae.history['acc'], 'bo', label='Training')
        plt.plot(range(self.epochs), self.trained_cae.history['val_acc'], 'b', label='Validation')
        plt.title('Training and validation accuracy')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.savefig(path+'training_validation_acc.png')
        plt.close()

    def save_model(self, path):
        """
            The model is saved in h5 and JSON formats, the weights are saved in h5 format and
            the summary of the model (the network layers) is saved in a txt file.

            Args:
                path(str): Path to the different files
        """
        # Model in HDF5 format
        self.cae.save(path+'model.h5')

        # Serialize weights to HDF5
        self.cae.save_weights(path+'model_weights.h5')

        # Model in JSON format
        with open(path+'model.json', 'w') as json_file:
            json_file.write(self.cae.to_json())

        # Model Summary
        with open(path+'model_summary.txt', 'w') as file:
            with redirect_stdout(file):
                self.cae.summary()
