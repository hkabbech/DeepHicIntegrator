#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./main.py <HIC_FILE> [--resolution INT] [--square_side INT] [--epochs INT]
                             [--batch_size INT] [--output PATH]

    Arguments:
        <HIC_FILE>                       Path to the Hi-C matrix file
                                        (.hic format).

    Options:
        -h, --help                      Show this.
        -r, INT, --resolution INT       Resolution used. [default: 25000]
        -n INT, --square_side INT       [default: 60]
        -e INT, --epochs INT            [default: 50]
        -b INT, --batch_size INT        [default: 128]
        -o PATH, --output PATH          [default: output/]
"""


# Third-party modules
from datetime import datetime
# from schema import Schema, And, Use, SchemaError
import os
from docopt import docopt
from hic2cool import hic2cool_convert
import cooler


# Local modules
from src.hic import Hic

# def autoencoder_network(input_img):
#     """
#     TO DO doctrings
#     """
#     # Encoder
#     conv1 = Conv2D(32, (3, 3), activation='relu', padding='same')(input_img)
#     pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
#     conv2 = Conv2D(64, (3, 3), activation='relu', padding='same')(pool1)
#     pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)
#     conv3 = Conv2D(128, (3, 3), activation='relu', padding='same')(pool2)
#     # Decoder
#     conv4 = Conv2D(128, (3, 3), activation='relu', padding='same')(conv3)
#     up1 = UpSampling2D((2, 2))(conv4)
#     conv5 = Conv2D(64, (3, 3), activation='relu', padding='same')(up1)
#     up2 = UpSampling2D((2, 2))(conv5)
#     decoded = Conv2D(1, (3, 3), activation='sigmoid', padding='same')(up2)
#     return decoded

if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    # check_args()
    HIC_FILE = ARGS['<HIC_FILE>']
    COOL_FILE = HIC_FILE.split('.hic')[0]+'.cool'
    RESOLUTION = int(ARGS['--resolution'])
    N = int(ARGS['--square_side'])
    EPOCHS = int(ARGS['--epochs'])
    BATCH_SIZE = int(ARGS['--batch_size'])
    OUTPUT_PATH = ARGS['--output']
    os.makedirs(OUTPUT_PATH+'model/', exist_ok=True)

    if not os.path.exists(COOL_FILE):
        ## using 0 triggers a multi-res output (.mcool format)
        # hic2cool_convert(HIC_FILE, COOL_FILE, 0)
        hic2cool_convert(HIC_FILE, COOL_FILE, RESOLUTION)


    TRAIN = Hic(cooler.Cooler(COOL_FILE), chrom=10)
    TRAIN.set_matrix()
    TRAIN.set_sub_matrices(N)

    





    TIME = datetime.now() - START_TIME
    print("\nTotal runtime: {}.{} seconds".format(TIME.seconds, TIME.microseconds))
