#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./main.py <HIC_FILE> [--resolution INT] [--square_side INT] [--epochs INT]
                             [--batch_size INT] [--output PATH] [--train INT] [--test INT]

    Arguments:
        <HIC_FILE>                      Path to the Hi-C matrix file
                                        (.hic format).

    Options:
        -h, --help                      Show this.
        -r, INT, --resolution INT       Resolution used. [default: 25000]
        -a INT, --train INT             [default: 1]
        -t INT, --test INT              [default: 20]
        -n INT, --square_side INT       [default: 60]
        -e INT, --epochs INT            [default: 50]
        -b INT, --batch_size INT        [default: 128]
        -o PATH, --output PATH          [default: output/]

"""


# Third-party modules
from datetime import datetime
# from schema import Schema, And, Use, SchemaError
import os
import random as rd
from contextlib import redirect_stdout
import numpy as np
from docopt import docopt
from tensorflow.python.client import timeline
import tensorflow as tf
from sklearn.model_selection import train_test_split
from keras.layers import Input, Conv2D, MaxPooling2D, UpSampling2D
from keras.models import Model
from keras.optimizers import RMSprop
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import cooler
from hic2cool import hic2cool_convert


# Local modules
from src.hic import Hic
from src.predict_hic import PredictHic

def autoencoder_network(input_img):
    """
    TO DO doctrings
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

if __name__ == "__main__":

    START_TIME = datetime.now()

    ### PARSE COMMAND LINE
    ######################
    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    # check_args()
    RESOLUTION = int(ARGS['--resolution'])
    HIC_FILE = ARGS['<HIC_FILE>']
    COOL_FILE = HIC_FILE.split('.hic')[0]+'.cool'
    # Conversion from .hic to .cool
    if not os.path.exists(COOL_FILE):
        hic2cool_convert(HIC_FILE, COOL_FILE, RESOLUTION)
    CHROM_TRAIN = int(ARGS['--train'])
    CHROM_TEST = int(ARGS['--test'])
    N = int(ARGS['--square_side'])
    EPOCHS = int(ARGS['--epochs'])
    BATCH_SIZE = int(ARGS['--batch_size'])
    OUTPUT_PATH = ARGS['--output']
    os.makedirs(OUTPUT_PATH+'model/', exist_ok=True)

    # HIC_FILE = 'data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic'
    # COOL_FILE = HIC_FILE.split('.hic')[0]+'.cool'
    # RESOLUTION = 25000
    # N = 60
    # EPOCHS = 50
    # BATCH_SIZE = 128
    # OUTPUT_PATH = 'results/03_04_2019_cooler/'
    # os.makedirs(OUTPUT_PATH+'model/', exist_ok=True)

    REDS = cm.get_cmap('Reds', 300)
    CMP = ListedColormap(np.vstack((np.array([1, 1, 1, 1]), REDS(np.linspace(0, 1, 300)))))


    ### TRAINING OF THE MODEL
    #########################
    TRAIN = Hic(cooler.Cooler(COOL_FILE), CHROM_TRAIN, N)
    TRAIN.set_matrix()
    TRAIN.plot_matrix(CMP, OUTPUT_PATH+'model')
    TRAIN.set_sub_matrices()

    TRAIN_X, VALID_X, TRAIN_GROUND, VALID_GROUND = train_test_split(TRAIN.sub_matrices,
                                                                    TRAIN.sub_matrices,
                                                                    test_size=0.2, random_state=13)
    INPUT_IMG = Input(shape=(N, N, 1))
    AUTOENCODER = Model(INPUT_IMG, autoencoder_network(INPUT_IMG))
    AUTOENCODER.compile(loss='mean_squared_error', optimizer=RMSprop(),
                        metrics=['accuracy'])
    AUTOENCODER.summary()
    AUTOENCODER_TRAIN = AUTOENCODER.fit(TRAIN_X, TRAIN_GROUND, batch_size=BATCH_SIZE, epochs=EPOCHS,
                                        verbose=1, validation_data=(VALID_X, VALID_GROUND))

    ### TEST OF THE MODEL
    #####################
    TEST = PredictHic(cooler.Cooler(COOL_FILE), CHROM_TEST, N)
    TEST.set_matrix()
    TEST.set_sub_matrices()
    TEST.set_predicted_sub_matrices(AUTOENCODER.predict(TEST.sub_matrices))
    TEST.construct_predicted_matrix()
    TEST.write_predicted_matrix(0.1, OUTPUT_PATH)

    RANDOM_INDEX_LIST = rd.sample(range(0, TEST.sub_matrices.shape[0]), 10)
    TEST.plot_ten_sub_matrices(CMP, OUTPUT_PATH, RANDOM_INDEX_LIST)
    TEST.plot_ten_predicted_submatrices(CMP, OUTPUT_PATH, RANDOM_INDEX_LIST)
    TEST.plot_matrix(CMP, OUTPUT_PATH)
    TEST.plot_predicted_matrix(CMP, OUTPUT_PATH)


    ### SAVE PARAMETERS
    ###################
    RUN_OPTIONS = tf.RunOptions(trace_level=tf.RunOptions.FULL_TRACE)
    RUN_METADATA = tf.RunMetadata()

    # Model Summary
    with open(OUTPUT_PATH+'model/model_summary.txt', 'w') as file:
        with redirect_stdout(file):
            AUTOENCODER.summary()

    # Timeline in JSON format
    with open(OUTPUT_PATH+'model/timeline.json', 'w') as file:
        tl = timeline.Timeline(RUN_METADATA.step_stats)
        file.write(tl.generate_chrome_trace_format())

    # Model in JSON format
    with open(OUTPUT_PATH+'model/model.json', 'w') as json_file:
        json_file.write(AUTOENCODER.to_json())

    # Model in HDF5 format
    AUTOENCODER.save(OUTPUT_PATH+'model/model.h5')

    # Serialize weights to HDF5
    AUTOENCODER.save_weights(OUTPUT_PATH+'model/model_weights.h5')

    # Loss plot
    plt.figure()
    plt.plot(range(EPOCHS), AUTOENCODER_TRAIN.history['loss'], 'b', label='Training')
    plt.plot(range(EPOCHS), AUTOENCODER_TRAIN.history['val_loss'], 'r', label='Validation')
    plt.title('Training and validation loss')
    plt.xlabel('epochs')
    plt.ylabel('loss')
    plt.legend()
    plt.savefig(OUTPUT_PATH+'model/training_validation_loss.png')
    plt.close()

    # Accuracy plot
    plt.figure()
    plt.plot(range(EPOCHS), AUTOENCODER_TRAIN.history['acc'], 'b', label='Training')
    plt.plot(range(EPOCHS), AUTOENCODER_TRAIN.history['val_acc'], 'r', label='Validation')
    plt.title('Training and validation accuracy')
    plt.xlabel('epochs')
    plt.ylabel('loss')
    plt.legend()
    plt.savefig(OUTPUT_PATH+'model/training_validation_acc.png')
    plt.close()

    # Function parameters
    TIME = datetime.now() - START_TIME
    with open(OUTPUT_PATH+'model/parameters.log', 'w') as file:
        file.write('Hi-C parameters:\n Resolution: {}\n Size sub-matrices: {}*{}\n\n'
                   .format(RESOLUTION, N, N))
        file.write('Train:\n Chromosome: {}\n Shape matrix: {}\n Shape sub_matrices: {}\n\n'
                   .format(CHROM_TRAIN, TRAIN.matrix.shape, TRAIN.sub_matrices.shape))
        file.write('Test:\n Chromosome: {}\n Shape matrix: {}\n Shape sub_matrices: {}\n\n'
                   .format(CHROM_TEST, TEST.matrix.shape, TEST.sub_matrices.shape))
        file.write('Autoencoder parameters:\n Epochs: {}\n Batch size: {}\n\n'
                   .format(EPOCHS, BATCH_SIZE))
        file.write('Running time: {}'.format(TIME))

    print("\nTotal runtime: {}.{} seconds".format(TIME.seconds, TIME.microseconds))
