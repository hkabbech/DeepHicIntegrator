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
from multiprocessing import cpu_count#, Pool
from datetime import datetime
#from schema import Schema, And, Use, SchemaError
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

def train_autoencoder():
    """
    TO DO doctrings
    """
    train_x, valid_x, train_ground, valid_ground = train_test_split(TRAIN.sub_matrices, TRAIN.sub_matrices, test_size=0.2,
                                                                    random_state=13)
    input_image = Input(shape=(N, N, INCHANNEL))

    autoencoder = Model(input_image, autoencoder_network(input_image))
    autoencoder.compile(loss='mean_squared_error', optimizer=RMSprop())
    autoencoder.summary()
    autoencoder_train = autoencoder.fit(train_x, train_ground, batch_size=BATCH_SIZE, epochs=EPOCHS,
                                        verbose=1, validation_data=(valid_x, valid_ground))
    plt.figure()
    plt.plot(range(EPOCHS), autoencoder_train.history['loss'], 'bo', label='Training loss')
    plt.plot(range(EPOCHS), autoencoder_train.history['val_loss'], 'b', label='Validation loss')
    plt.title('Training and validation loss')
    plt.legend()
    plt.savefig('training_validation_loss.png')

    TEST.set_predicted_sub_matrices(autoencoder.predict(TEST))


if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    # check_args()
    TRAIN = Hic(ARGS['<train_hic>'], int(ARGS['--resolution']))
    TEST = Hic(ARGS['<test_hic>'], int(ARGS['--resolution']))
    N = int(ARGS['--square_side'])
    EPOCHS = int(ARGS['--epochs'])
    BATCH_SIZE = int(ARGS['--batch_size'])
    INCHANNEL = int(ARGS['--inchannel'])
    NB_PROC = cpu_count() if int(ARGS["--cpu"]) == 0 else int(ARGS["--cpu"])

    TRAIN.set_matrix(added_lines=0, deleted_lines=21)
    TRAIN.set_sub_matrices(N, N)

    TEST.set_matrix(added_lines=1, deleted_lines=0)
    TEST.set_sub_matrices(N, N)

    train_autoencoder()


    ### Plots
    #########
    INDICES_LIST = [85, 86, 258, 311, 312, 313, 502, 624, 908, 1203]
    TRAIN.plot_ten_sub_matrices(INDICES_LIST)
    TRAIN.plot_matrix()
    TEST.plot_ten_sub_matrices(INDICES_LIST, chr_type="predicted")
    TEST.plot_matrix(chr_type="predicted")

    print("\nTotal runtime: {} seconds".format(datetime.now() - START_TIME))
