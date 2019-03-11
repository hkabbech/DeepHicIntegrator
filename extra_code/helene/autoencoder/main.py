#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./main.py <train_hic> <test_hic> [--resolution INT] [--square_side INT] [--epochs INT]
                                         [--batch_size INT] [--inchannel INT] [--cpu INT]
                                         [--output PATH]

    Arguments:
        <train_hic>                      Path to the Hi-C matrix file.
        <test_hic>                       Path to the Hi-C matrix file.

    Options:
        -h, --help                      Show this.
        -r, INT, --resolution INT       Resolution used for the sparse matrix.
                                        [default: 25000]
        -s INT, --square_side INT       [default: 60]
        -e INT, --epochs INT            [default: 50]
        -b INT, --batch_size INT        [default: 128]
        -i INT, --inchannel INT         [default: 1]
        -c INT, --cpu INT               Number of cpus to use for parallelisation.
                                        By default using all available (0).
                                        [default: 0]
        -o PATH, --output PATH          [default: output/]
"""


# Third-party modules
# from multiprocessing import cpu_count#, Pool
from datetime import datetime
# from schema import Schema, And, Use, SchemaError
import os
from contextlib import redirect_stdout
from docopt import docopt
from sklearn.model_selection import train_test_split
from keras.layers import Input, Conv2D, MaxPooling2D, UpSampling2D
from keras.models import Model
from keras.optimizers import RMSprop
import matplotlib.pyplot as plt



# Local modules
from src.hic import Hic

def autoencoder_network(input_img):
    """
    TO DO doctrings
    """
    #encoder
    conv1 = Conv2D(32, (3, 3), activation='relu', padding='same')(input_img)
    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
    conv2 = Conv2D(64, (3, 3), activation='relu', padding='same')(pool1)
    pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)
    conv3 = Conv2D(128, (3, 3), activation='relu', padding='same')(pool2)

    #decoder
    conv4 = Conv2D(128, (3, 3), activation='relu', padding='same')(conv3)
    up1 = UpSampling2D((2, 2))(conv4)
    conv5 = Conv2D(64, (3, 3), activation='relu', padding='same')(up1)
    up2 = UpSampling2D((2, 2))(conv5)
    decoded = Conv2D(1, (3, 3), activation='sigmoid', padding='same')(up2)
    return decoded

def plot_loss(autoencoder_train):
    """
    TO DO doctrings
    """
    plt.figure()
    plt.plot(range(EPOCHS), autoencoder_train.history['loss'], 'bo', label='Training loss')
    plt.plot(range(EPOCHS), autoencoder_train.history['val_loss'], 'b', label='Validation loss')
    plt.title('Training and validation loss')
    plt.legend()
    plt.savefig(OUTPUT_PATH+'plots/training_validation_loss.png')

def save_summary(model):
    """
    TO DO doctrings
    """
    with open(OUTPUT_PATH+'files/model_summary.txt', 'w') as file:
        with redirect_stdout(file):
            model.summary()

def save_parameters():
    """
    TO DO doctrings
    """
    with open(OUTPUT_PATH+'files/parameters.log', 'w') as file:
        file.write('Hi-C parameters:\n Resolution: {}\n Size sub-matrices: {}*{}\n\n'
                   .format(RESOLUTION, N, N))
        file.write('Train:\n Filename: {}\n Added lines: {}\n Deleted lines: {}\n\n'
                   .format(TRAIN_FILENAME, ADDED_LINES_TRAIN, DELETED_LINES_TRAIN))
        file.write('Test:\n Filename: {}\n Added lines: {}\n Deleted lines: {}\n\n'
                   .format(TEST_FILENAME, ADDED_LINES_TEST, DELETED_LINES_TEST))
        file.write('Autoencoder parameters:\n Epochs: {}\n Batch size: {}\n InChannel: {}\n\n'
                   .format(EPOCHS, BATCH_SIZE, INCHANNEL, TRAIN_FILENAME, TEST_FILENAME))
        file.write('Running time: {}'.format(TIME_TOTAL))

if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    # check_args()
    TRAIN_FILENAME = ARGS['<train_hic>']
    TEST_FILENAME = ARGS['<test_hic>']
    RESOLUTION = int(ARGS['--resolution'])
    N = int(ARGS['--square_side'])
    EPOCHS = int(ARGS['--epochs'])
    BATCH_SIZE = int(ARGS['--batch_size'])
    INCHANNEL = int(ARGS['--inchannel'])
    # NB_PROC = cpu_count() if int(ARGS["--cpu"]) == 0 else int(ARGS["--cpu"])
    OUTPUT_PATH = ARGS['--output']
    os.makedirs(OUTPUT_PATH+'plots/', exist_ok=True)
    os.makedirs(OUTPUT_PATH+'files/', exist_ok=True)


    ### Autoencoder Training
    ########################

    TRAIN = Hic(TRAIN_FILENAME)
    # ADDED_LINES_TRAIN, DELETED_LINES_TRAIN = 0, 21
    ADDED_LINES_TRAIN, DELETED_LINES_TRAIN = 0, 0
    TRAIN.set_matrix(RESOLUTION, ADDED_LINES_TRAIN, DELETED_LINES_TRAIN)
    TRAIN.set_sub_matrices(N, N)

    TRAIN_X, VALID_X, TRAIN_GROUND, VALID_GROUND = train_test_split(TRAIN.sub_matrices,
                                                                    TRAIN.sub_matrices,
                                                                    test_size=0.2, random_state=13)
    INPUT_IMG = Input(shape=(N, N, INCHANNEL))
    AUTOENCODER = Model(INPUT_IMG, autoencoder_network(INPUT_IMG))
    AUTOENCODER.compile(loss='mean_squared_error', optimizer=RMSprop())
    AUTOENCODER_TRAIN = AUTOENCODER.fit(TRAIN_X, TRAIN_GROUND, batch_size=BATCH_SIZE, epochs=EPOCHS,
                                        verbose=1, validation_data=(VALID_X, VALID_GROUND))
    plot_loss(AUTOENCODER_TRAIN)
    save_summary(AUTOENCODER)


    ### Test
    ########
    TEST = Hic(TEST_FILENAME)
    # ADDED_LINES_TEST, DELETED_LINES_TEST = 1, 0
    ADDED_LINES_TEST, DELETED_LINES_TEST = 0, 0
    TEST.set_matrix(RESOLUTION, ADDED_LINES_TEST, DELETED_LINES_TEST)
    #TEST.set_matrix(RESOLUTION, added_lines=1, deleted_lines=0)
    TEST.set_matrix(RESOLUTION, added_lines=0, deleted_lines=0)
    TEST.set_sub_matrices(N, N)
    TEST.set_predicted_sub_matrices(AUTOENCODER.predict(TEST.sub_matrices))
    TEST.set_reconstructed_matrix(N)
    TEST.save_reconstructed_matrix(OUTPUT_PATH, RESOLUTION)

    # INDICES_LIST = [85, 86, 258, 311, 312, 313, 502, 624, 908, 1203]
    INDICES_LIST = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    TEST.plot_ten_sub_matrices("true", INDICES_LIST, OUTPUT_PATH)
    TEST.plot_ten_sub_matrices("predicted", INDICES_LIST, OUTPUT_PATH)
    TEST.plot_matrix("true", OUTPUT_PATH)
    TEST.plot_matrix("predicted", OUTPUT_PATH)

    TIME_TOTAL = datetime.now() - START_TIME
    print("\nTotal runtime: {} seconds".format(TIME_TOTAL))
    save_parameters()
