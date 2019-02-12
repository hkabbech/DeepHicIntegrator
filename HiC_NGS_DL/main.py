#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./main.py <file> <path> [--cpu INT]

    Arguments:
        <file>                          Path to the Hi-C matrix file.
        <path>                          Path to the folder containing all NGS data relative
                                        to the HiC matrix.

    Options:
        -h, --help                      Show this.
        -c INT, --cpu INT               Number of cpus to use for parallelisation. By default
                                        using all available (0).
                                        [default: 0]
"""


__authors__ = "Eduardo Gade Gusmão and Hélène Kabbech"


# Third-party modules
from multiprocessing import Pool, cpu_count
from datetime import datetime
from schema import Schema, And, Use, SchemaError
from docopt import docopt
import os
import pandas as pd

# Local modules
from src.hic import Hic


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command linSe.
    """
    schema = Schema({
        '<file>': Use(open),
        '<path>': os.path.exists,
        '--cpu': And(Use(int), lambda n: 0 <= n <= cpu_count(),
                     error='--cpus should be integer 1 <= N <= ' + str(cpu_count()))
        })
    try:
        schema.validate(ARGS)
    except SchemaError as err:
      exit(err)


if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()
    NGS = ARGS['<path>']
    CPU = ARGS['--cpu']
    HIC = Hic(ARGS['<file>'])

    # Example How to get a value from 2 given bases :
    base_1, base_2 = 75000, 350000
    print(HIC.get_value(base_1, base_2))

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
