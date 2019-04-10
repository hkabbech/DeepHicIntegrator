"""
.. module:: Autoencoder
   :synopsis: This module implements the Autoencoder class.
"""

# Third-party modules
from contextlib import redirect_stdout
from sklearn.model_selection import train_test_split
from keras.layers import Conv2D, MaxPooling2D, UpSampling2D
from keras.models import Model
from keras.optimizers import RMSprop
import matplotlib.pyplot as plt


def build_network(input_img):
    """
        Layers of the Autoencoder network.

        Args:
            input_img(Input class): Input of the network

        Returns:
            decoded(Conv2D keras layer): The last decoded output layer of the Autoencoder
    """
    # Encoder
    conv1 = Conv2D(32, (3, 3), activation='relu', padding='same')(input_img)
    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
    conv2 = Conv2D(64, (3, 3), activation='relu', padding='same')(pool1)
    pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)
    conv3 = Conv2D(128, (3, 3), activation='relu', padding='same')(pool2)
    # Decoder
    conv4 = Conv2D(128, (3, 3), activation='relu', padding='same')(conv3)
    up1 = UpSampling2D((2, 2))(conv4)
    conv5 = Conv2D(64, (3, 3), activation='relu', padding='same')(up1)
    up2 = UpSampling2D((2, 2))(conv5)
    decoded = Conv2D(1, (3, 3), activation='sigmoid', padding='same')(up2)
    return decoded


class Autoencoder:
    """
    .. class:: Autoencoder
        This class groups attributes and functions about the Autoencoder.

    Attributes:
        chr_train(Hic object): Instance of the HicPredict class - Chromosome for training the model
        model(Model object): Instance of the Model keras class - Model of the Autoencoder
        trained_model(Model object): Instance of the Model keras class - The trained Model

    """

    def __init__(self, input_img, chr_train):
        self.chr_train = chr_train
        self.model = Model(input_img, build_network(input_img))
        self.trained_model = None

    def compile_model(self, loss_function='mean_squared_error', metrics_list=['accuracy']):
        """
            Compilation of the Autoencoder model.

            Args:
                loss_function(str): The loss function to use
                matrics_list(list): The metrics to generate during compilation
        """
        self.model.compile(loss=loss_function, optimizer=RMSprop(),
                           metrics=metrics_list)
        self.model.summary()

    def train_model(self, batch_size, epochs):
        """
            Training of the Autoencoder model and set of the trained_model attribute.

            Args:
                batch_size(int): Size of a batch
                epochs(int): Number of epochs
        """
        train_x, valid_x, train_ground, valid_ground = train_test_split(self.chr_train.sub_matrices,
                                                                        self.chr_train.sub_matrices,
                                                                        test_size=0.2,
                                                                        random_state=13)

        self.trained_model = self.model.fit(train_x, train_ground, verbose=1,
                                            batch_size=batch_size, epochs=epochs,
                                            validation_data=(valid_x, valid_ground))

    def plot_loss_curve(self, epochs, path):
        """
            Plot of the loss curve in a file.

            Args:
                epochs(int): Number of epochs
                path(str): Path of the output plot
        """
        plt.plot(range(epochs), self.trained_model.history['loss'], 'bo', label='Training')
        plt.plot(range(epochs), self.trained_model.history['val_loss'], 'b', label='Validation')
        plt.title('Training and validation loss')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.savefig(path+'training_validation_loss.png')
        plt.close()

    def plot_accuracy_curve(self, epochs, path):
        """
            Plot of the accuracy curve in a file.

            Args:
                epochs(int): Number of epochs
                path(str): Path of the output plot
        """
        plt.plot(range(epochs), self.trained_model.history['acc'], 'bo', label='Training')
        plt.plot(range(epochs), self.trained_model.history['val_acc'], 'b', label='Validation')
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
        self.model.save(path+'model.h5')

        # Serialize weights to HDF5
        self.model.save_weights(path+'model_weights.h5')

        # Model in JSON format
        with open(path+'model.json', 'w') as json_file:
            json_file.write(self.model.to_json())

        # Model Summary
        with open(path+'model_summary.txt', 'w') as file:
            with redirect_stdout(file):
                self.model.summary()

    def save_parameters(self, path, resolution, side, chr_test, epochs, batch_size, time):
        """
            All the parameters of the training and test of the model are saved in a log file.

            Args:
                path(str): Path of the log file
        """
        with open(path+'parameters.log', 'w') as file:
            file.write('Hi-C parameters:\n Resolution: {}\n Size sub-matrices: {}*{}\n\n'
                       .format(resolution, side, side))
            file.write('Train:\n Chromosome: {}\n Shape matrix: {}\n Shape sub_matrices: {}\n\n'
                       .format(self.chr_train.chrom, self.chr_train.matrix.shape,
                               self.chr_train.sub_matrices.shape))
            file.write('Test:\n Chromosome: {}\n Shape matrix: {}\n Shape sub_matrices: {}\n\n'
                       .format(chr_test.chrom, chr_test.matrix.shape, chr_test.sub_matrices.shape))
            file.write('Autoencoder parameters:\n Epochs: {}\n Batch size: {}\n'
                       .format(epochs, batch_size))
            file.write(' Optimizer: {}\n Loss function: {}\n Metrics: {}\n\n'
                       .format(self.model.optimizer, self.model.loss, self.model.metrics))
            file.write('Running time: {}'.format(time))
