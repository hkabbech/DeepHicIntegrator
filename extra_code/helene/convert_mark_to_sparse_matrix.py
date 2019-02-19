#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./run_convert_mark_to_sparse_matrix.py <info_file> <chr_sizes_file> <input_path>
                                               <output_path> [--resolution INT] [--threshold INT]
                                               [--chr_list STR] [--cpu INT]

    Arguments:
        <info_file>                         Path to the bam information file
                                            (bam_name\treads_number).
        <chr_sizes_file>                    Path to the chromosome sizes file
                                            (chr_name\tchr_size).
        <input_path>                        Path to the folder containing bam files.
        <output_path>                       Path to the output folder.

    Options:
        -r INT, --resolution INT            Resolution value. [default: 25000]
        -m INT, --threshold INT             Minimum reads threshold value. [default: 2500]
        -l STR, --chr_list STR              List of chromosome names separeted by "_".
                                            [default: all]
        -c INT, --cpu INT                   Number of cpus to use for parallelisation.
                                            By default using all available (0).
                                            [default: 0]
        -h, --help                          Show this.
"""

# Third-party modules
# import subprocess
# from multiprocessing import Pool, cpu_count
from multiprocessing import cpu_count
# from functools import partial
from datetime import datetime
# from tqdm import tqdm
import os
from docopt import docopt
from schema import Schema, And, Use, SchemaError
from pysam import Samfile

def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        '<info_file>': Use(open),
        '<chr_sizes_file>': And(Use(open)),
        '<input_path>': os.path.exists,
        '--resolution': And(Use(int), lambda n: n >= 0,
                            error='--resolution=INT should be a positive integer'),
        '--threshold': And(Use(int), lambda n: n >= 0,
                           error='--threshold=INT should be a positive integer'),
        '--chr_list': And(Use(str)),
        '--cpu': And(Use(int), lambda n: 0 <= n <= cpu_count(),
                     error='--cpus=NUM should be integer 1 <= N <= ' + str(cpu_count())),
        # We skip the output path check because it may does not exist yet.
        object: object
        })
    try:
        schema.validate(ARGS)
    except SchemaError as err:
        exit(err)

def create_dictionary(filename):
    """
        Create a dictionary with key = str and value = int.

        Args:
            filname(str): Path to the file containing the informations
                          (format: column_1\tcolumn_2).
        Returns:
            dict: The created dictionary.
    """
    info_dict = {}
    with open(filename, "r") as file:
        for line in file:
            contents = line.split("\t")
            info_dict[contents[0]] = int(contents[1])
    return info_dict

def fetch_total_reads(bamfile, region):
    """
        fsdfsd

        Args:
            bamfile (file): fs
            region (list): sdf
        Returns:
            int: dsf
    """
    count = 0
    for _ in bamfile.fetch(region[0], region[1], region[2]):
        count += 1
    return count

if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()
    BAM_READS_DICT = create_dictionary(ARGS['<info_file>'])
    CHR_SIZES_DICT = create_dictionary(ARGS['<chr_sizes_file>'])
    RESOLUTION = ARGS['--resolution']
    THRESHOLD = ARGS['--threshold']
    CHR_LIST = ["chr"+str(k) for k in range(1, 23)] + ["chrX"] if ARGS["--chr_list"] == "all"\
                else ARGS["--chr_list"].split("_")
    NB_PROC = cpu_count() if int(ARGS["--cpu"]) == "all" else int(ARGS["--cpu"])
    # Bam Loop
    for BAM_NAME, TOTAL_READS in BAM_READS_DICT.items():
        BAM_FILE = Samfile(ARGS['<input_path>'] + "/" + BAM_NAME + ".bam", "rb")
        OUTPUT_LOCATION = ARGS['<output_path>'] + "/" + BAM_NAME + "/"
        os.makedirs(OUTPUT_LOCATION, exist_ok=True)

        # Chromosome Loop
        for CHROMOSOME in CHR_LIST:
            OUTPUT_FILE = open(OUTPUT_LOCATION + CHROMOSOME + "_" + BAM_NAME + ".txt", "w")

            # Iterating on genome
            for START_1 in range(0, CHR_SIZES_DICT[CHROMOSOME], RESOLUTION):

                # Fetching locations
                END_1 = START_1 + RESOLUTION
                if END_1 > CHR_SIZES_DICT[CHROMOSOME]:
                    END_1 = CHR_SIZES_DICT[CHROMOSOME]

                # Fetching reads
                TOTAL_READS_1 = fetch_total_reads(BAM_FILE, [CHROMOSOME, START_1, END_1])
                if (TOTAL_READS_1 == 0 or TOTAL_READS_1 < THRESHOLD):
                    TOTAL_READS_1 = 1

                # Iterating on genome
                for START_2 in range(START_1, CHR_SIZES_DICT[CHROMOSOME], RESOLUTION):

                    # Fetching locations
                    END_2 = START_2 + RESOLUTION
                    if END_2 > CHR_SIZES_DICT[CHROMOSOME]:
                        END_2 = CHR_SIZES_DICT[CHROMOSOME]

                    # Fetching reads
                    TOTAL_READS_2 = fetch_total_reads(BAM_FILE, [CHROMOSOME, START_2, END_2])
                    if TOTAL_READS_2 == 0:
                        TOTAL_READS_2 = 1
                    if (TOTAL_READS_1 < THRESHOLD\
                        and TOTAL_READS_2 < THRESHOLD):
                        continue

                    # Writing to file
                    TOTAL = (TOTAL_READS_1 * TOTAL_READS_2) / TOTAL_READS
                    OUTPUT_FILE.write("{}\t{}\t{}\t{}\n".\
                                     format(CHROMOSOME, str(START_1), str(START_2), str(TOTAL)))
            # Closing files
            OUTPUT_FILE.close()
        BAM_FILE.close()

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
