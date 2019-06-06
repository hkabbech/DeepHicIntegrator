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
        This class groups attributes and functions for building compiling, training and saving the
        parameters of an Autoencoder.

    Attributes:
        chr_train (Hic(Matrix) object): Chromosome for training the model
        img_size (int): Square side of a numpy array sub-matrix
        epochs (int): Number of epochs for training the cae
        batch_size (int): Batch size
        cae (keras model object): Model of the whole convolutional autoencoder (cae)
        trained_cae (keras model object): The trained cae
        encoder (keras model object): Encoder part of the cae
        decoder (keras model object): Decoder part of the cae
    """

    def __init__(self, chr_train, img_size, epochs, batch_size):
        self.chr_train = chr_train
        self.img_size = img_size
        self.epochs = epochs
        self.batch_size = batch_size
        self.cae = None
        self.trained_cae = None
        self.encoder = None
        self.decoder = None

    def set_models(self):
        """
            Building of the convolutional autoencoder and setting of the encoder, decoder and cae
            attributes.

            Args:
                input_img(Input class): Input image for the encoder
        """
        # Encoder part
        input_img = Input(shape=(self.img_size, self.img_size, 1))
        layer = Conv2D(128, (3, 3), activation='relu', padding='same')(input_img)
        layer = MaxPooling2D((2, 2), padding='same')(layer)
        layer = Conv2D(64, (3, 3), activation='relu', padding='same')(layer)
        layer = MaxPooling2D((2, 2), padding='same')(layer)
        latent = Conv2D(32, (3, 3), activation='relu', padding='same')(layer)
        self.encoder = Model(input_img, latent, name='Encoder')
        print("\nEncoder Model")
        self.encoder.summary()
        # Decoder part
        input_ls = Input(shape=(self.img_size/4, self.img_size/4, 32))
        layer = Conv2D(32, (3, 3), activation='relu', padding='same')(input_ls)
        layer = UpSampling2D((2, 2))(layer)
        layer = Conv2D(64, (3, 3), activation='relu', padding='same')(layer)
        layer = UpSampling2D((2, 2))(layer)
        output = Conv2D(1, (3, 3), activation='sigmoid', padding='same')(layer)
        self.decoder = Model(input_ls, output, name='Decoder')
        print("\nDecoder Model")
        self.decoder.summary()

        self.cae = Model(input_img, self.decoder(self.encoder(input_img)), name='Autoencoder')
        print("\nAutoencoder Model")
        self.cae.summary()

    def compile(self, loss_function='mean_squared_error', metrics_list=['accuracy']):
        """
            Compilation of the cae model.

            Args:
                loss_function(str): The loss function
                matrics_list(list): The metrics to generate during compilation
        """
        self.cae.compile(loss=loss_function, optimizer=RMSprop(), metrics=metrics_list)

    def train(self):
        """
            Training of the cae model and setting of the trained_cae attribute.
        """
        train_x, valid_x, train_ground, valid_ground = train_test_split(self.chr_train,
                                                                        self.chr_train,
                                                                        test_size=0.2,
                                                                        random_state=13)

        self.trained_cae = self.cae.fit(train_x, train_ground, verbose=1,
                                        batch_size=self.batch_size, epochs=self.epochs,
                                        validation_data=(valid_x, valid_ground))

    def save_model(self, path):
        """
            The model (h5 and JSON formats), the weights (h5 format) and the summary of the models
            are saved in a txt file.

            Args:
                path(str): Path to the different files
        """
        # Model in HDF5 format
        self.encoder.save(path+'encoder.h5')
        self.decoder.save(path+'decoder.h5')
        self.cae.save(path+'autoencoder.h5')

        # Serialize weights to HDF5
        self.encoder.save_weights(path+'encoder_weights.h5')
        self.decoder.save_weights(path+'decoder_weights.h5')
        self.cae.save_weights(path+'autoencoder_weights.h5')

        # Model in JSON format
        with open(path+'model.json', 'w') as json_file:
            json_file.write(self.cae.to_json())

        # Model Summary
        with open(path+'model_summary.txt', 'w') as file:
            with redirect_stdout(file):
                print('Autoencoder Model')
                self.cae.summary()
                print('\n\nEncoder Model')
                self.encoder.summary()
                print('\n\nDecoder Model')
                self.decoder.summary()

    def plot_loss_curve(self, path):
        """
            Plot of the loss curve in a file.

            Args:
                path(str): Path of the output plot
        """
        plt.plot(range(self.epochs), self.trained_cae.history['loss'], 'bo', label='Training')
        plt.plot(range(self.epochs), self.trained_cae.history['val_loss'], 'b', label='Validation')
        plt.title('Training and validation loss')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.savefig(path+'/training_validation_loss.pdf')
        plt.close()
        with open(path+'/training_loss.tsv', 'w') as file:
            for i, loss in enumerate(self.trained_cae.history['loss']):
                file.write('{}\t{}\n'.format(i+1, loss))
        with open(path+'/validation_loss.tsv', 'w') as file:
            for i, loss in enumerate(self.trained_cae.history['val_loss']):
                file.write('{}\t{}\n'.format(i+1, loss))

    def plot_accuracy_curve(self, path):
        """
            Plot of the accuracy curve in a file.

            Args:
                path(str): Path of the output plot
        """
        plt.plot(range(self.epochs), self.trained_cae.history['acc'], 'bo', label='Training')
        plt.plot(range(self.epochs), self.trained_cae.history['val_acc'], 'b', label='Validation')
        plt.title('Training and validation accuracy')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.savefig(path+'training_validation_acc.pdf')
        plt.close()
