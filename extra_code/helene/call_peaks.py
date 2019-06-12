#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    Call peaks in a matrix.
    Autor: Hélène Kabbech, Institute of Pathologie, Medical Center University of Göttingen (Germany)

    Usage:
        ./call_peaks.py <FILE> [--resolution INT] [--percentile INT] [--diagonale_cutoff INT]
                                                   [--output PATH]
                                                   [--help]

    Arguments:
        <FILE>                              Path of the Hi-C matrix file (.hic format)

    Options:
        -r INT, --resolution INT            Resolution representing the number of pair-ended reads
                                            spanning between a pair of bins. [default: 100000]
        -p INT, --percentile INT            Percentile threshold. [default: 98]
        -d INT, --diagonale_cutoff INT      Diagonale cutoff. [default: 800000]
        -o PATH, --output PATH              Output path [default: results/]
        -h, --help                          Show this

"""

# Third-party modules
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
from docopt import docopt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == "__main__":

    RED_CMP = ListedColormap(np.vstack((np.array([1, 1, 1, 1]),
                                        cm.get_cmap('Reds', 300)(np.linspace(0, 1, 300)))))

    # Arguments
    ARGS = docopt(__doc__)
    FILENAME = ARGS['<FILE>']
    RESOLUTION = int(ARGS['--resolution'])
    PERCENTILE = int(ARGS['--percentile'])
    DIAGONALE_CUTOFF = int(ARGS['--diagonale_cutoff'])


    # Integrated matrix in a dataframe
    DF = pd.read_csv(FILENAME, header=None, sep='\t')
    DF.columns = ['chr', 'base_1', 'base_2', 'value']

    # Conversion from pandas to numpy array
    DATA = DF['value']
    ROW = ((DF['base_1'] / RESOLUTION)).astype(int)
    COL = ((DF['base_2'] / RESOLUTION)).astype(int)
    SIZE = int(max(max(DF['base_2']), max(DF['base_1'])) / RESOLUTION)
    MATRIX = coo_matrix((DATA, (ROW, COL)), shape=(SIZE+1, SIZE+1)).toarray()

    # Call peaks
    THRESHOLD = np.percentile(DF['value'], PERCENTILE)
    DF_THRESHOLDED = DF[DF['value'] > THRESHOLD]
    for i in range(0, DIAGONALE_CUTOFF+RESOLUTION, RESOLUTION):
        DF_THRESHOLDED = DF_THRESHOLDED[DF_THRESHOLDED['base_2'] != DF_THRESHOLDED['base_1']+i]
    DF_THRESHOLDED.to_csv('{}_call_peaks.bed'.format(FILENAME.split('.')[0]), sep='\t', index=None)


    ROW_THRESHOLDED = ((DF_THRESHOLDED['base_1'] / RESOLUTION)).astype(int)
    COL_THRESHOLDED = ((DF_THRESHOLDED['base_2'] / RESOLUTION)).astype(int)

    # Plot
    FIG = plt.figure(figsize=(12, 12))
    AXES = plt.subplot(111, aspect='equal')
    # Plot of the integrated matrix
    IMG = AXES.matshow(MATRIX, cmap=RED_CMP)
    # Circles on the integrated matrix
    AXES.scatter(COL_THRESHOLDED, ROW_THRESHOLDED, s=20, alpha=1, facecolors='none',
                 edgecolors='black')
    DIVIDER = make_axes_locatable(AXES)
    CAX = DIVIDER.append_axes("right", size="2%", pad=0.15)
    plt.colorbar(IMG, cax=CAX)
    plt.subplots_adjust(left=0.07, bottom=0, right=0.95, top=0.91, wspace=0, hspace=0)
    AXES.axis('off')
    plt.show()
    plt.savefig('{}_call_peaks.pdf'.format(FILENAME.split('.')[0]))
    plt.close()
