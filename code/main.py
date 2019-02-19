#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./main.py <file_1> <file_2> [--cpu INT]

    Arguments:
        <file_1>                          Path to the Hi-C matrix file.
        <file_2>                          Path to the NGS matrix file.

    Options:
        -h, --help                      Show this.
        -c INT, --cpu INT               Number of cpus to use for parallelisation. By default
                                        using all available (0).
                                        [default: 0]
"""


__authors__ = "Eduardo Gade Gusmão and Hélène Kabbech"


# Third-party modules
from multiprocessing import cpu_count#, Pool
from datetime import datetime
import os
from schema import Schema, And, Use, SchemaError
from docopt import docopt

# Local modules
from src.matrix import Matrix


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        '<file_1>': Use(open),
        '<file_2>': Use(open),
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
    HIC = Matrix(ARGS['<file_1>'])
    NGS = Matrix(ARGS['<file_2>'])
    CPU = ARGS['--cpu']
    

    # Example How to get a value from 2 given bases :
    BASE_1, BASE_2 = 75000, 350000
    print(HIC.get_value(BASE_1, BASE_2))

     # Example How to get a value from 2 given bases :
    BASE_1, BASE_2 = 23050000, 58100000
    print(NGS.get_value(BASE_1, BASE_2))

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
