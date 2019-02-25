#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./main.py <hic_file> <ngs_file> [--resolution INT] [--cpu INT]

    Arguments:
        <hic_file>                      Path to the Hi-C matrix file.
        <ngs_file>                      Path to the NGS matrix file.

    Options:
        -h, --help                      Show this.
        -r, INT, --resolution INT       Resolution used for the sparse matrix.
                                        [default: 25000]
        -c INT, --cpu INT               Number of cpus to use for parallelisation.
                                        By default using all available (0).
                                        [default: 0]
"""


__authors__ = "Eduardo Gade Gusmão and Hélène Kabbech"


# Third-party modules
from multiprocessing import cpu_count#, Pool
from datetime import datetime
from schema import Schema, And, Use, SchemaError
from docopt import docopt

# Local modules
from src.hic import Hic


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        '<hic_file>': Use(open),
        '<ngs_file>': Use(open),
        '--resolution': And(Use(int), lambda n: n >= 0,
                            error='--resolution=INT should be a positive integer'),
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
    HIC = Hic(ARGS['<hic_file>'], int(ARGS['--resolution']))
    NGS = ARGS['<ngs_file>']
    NB_PROC = cpu_count() if int(ARGS["--cpu"]) == 0 else int(ARGS["--cpu"])



    # Example How to get a value from 2 given bases :
    BASE_1, BASE_2 = 75000, 350000
    print(HIC.get_value('chr20', BASE_1, BASE_2))

    # from datetime import datetime
    # START_TIME = datetime.now()
    # h.set_matrix()
    # print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
