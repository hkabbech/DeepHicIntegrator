#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./convert_mark_to_sparse_matrix.py <info_file> <chr_sizes_file> <input_path>
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
                                            By default using all available.
                                            [default: 0]
        -h, --help                          Show this.
"""

# Third-party modules
from multiprocessing import Pool, cpu_count
from functools import partial
from datetime import datetime
import os
from tqdm import tqdm
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
            filname (str): Path to the file containing the informations
                          (format: column_1\tcolumn_2)
        Returns:
            dict: The created dictionary
    """
    info_dict = {}
    with open(filename, "r") as file:
        for line in file:
            contents = line.split("\t")
            info_dict[contents[0]] = int(contents[1])
    return info_dict

def fetch_total_reads(bamfile, region):
    """
        Fetch the total number of reads in a bam file.

        Args:
            bamfile (file): Bam file
            region (list): Region to fetch in the bam file [chromosome, start, end]
        Returns:
            int: Total number of reads
    """
    count = 0
    for _ in bamfile.fetch(region[0], region[1], region[2]):
        count += 1
    return count

def process(job):
    """
        Create a sparse matrix file from a given bam file.

        Args:
            job (dict): Dictionary containing variables for the surrent job :
                        chromosome, bam_filename, output_filename and total_reads
    """
    bam_file = Samfile(job["bam_filename"], "rb")
    output_file = open(job["output_filename"], "w")


    # Iterating on genome
    for start_1 in range(0, CHR_SIZES_DICT[job["chromosome"]], RESOLUTION):

        # Fetching locations
        end_1 = start_1 + RESOLUTION
        if end_1 > CHR_SIZES_DICT[job["chromosome"]]:
            end_1 = CHR_SIZES_DICT[job["chromosome"]]

        # Fetching reads 1
        total_reads_1 = fetch_total_reads(bam_file, [job["chromosome"], start_1, end_1])
        if (total_reads_1 == 0 or total_reads_1 < THRESHOLD):
            total_reads_1 = 1

        # Iterating on genome
        for start_2 in range(start_1, CHR_SIZES_DICT[job["chromosome"]], RESOLUTION):

            # Fetching locations
            end_2 = start_2 + RESOLUTION
            if end_2 > CHR_SIZES_DICT[job["chromosome"]]:
                end_2 = CHR_SIZES_DICT[job["chromosome"]]

            # Fetching reads 2
            total_reads_2 = fetch_total_reads(bam_file, [job["chromosome"], start_2, end_2])
            if total_reads_2 == 0:
                total_reads_2 = 1
            if (total_reads_1 < THRESHOLD and total_reads_2 < THRESHOLD):
                continue
            else:
                total = (total_reads_1 * total_reads_2) / job["total_reads"]
                output_file.write("{}\t{}\t{}\t{}\n".\
                                 format(job["chromosome"], str(start_1), str(start_2), str(total)))
    # Closing files
    bam_file.close()
    output_file.close()


if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()
    BAM_READS_DICT = create_dictionary(ARGS['<info_file>'])
    CHR_SIZES_DICT = create_dictionary(ARGS['<chr_sizes_file>'])
    RESOLUTION = int(ARGS['--resolution'])
    THRESHOLD = int(ARGS['--threshold'])
    CHR_LIST = ["chr"+str(k) for k in range(1, 23)] + ["chrX"] if ARGS["--chr_list"] == "all"\
               else ARGS["--chr_list"].split("_")
    NB_PROC = cpu_count() if int(ARGS["--cpu"]) == 0 else int(ARGS["--cpu"])


    ## Multiprocessing
    ##################

    # Creating list for the multiprocessing
    # Each element of the list corresponds to a dictionary containing specific variables for a job
    LIST_FOR_MULTIPROCESSING = []
    # Loop crossing all bam name
    for BAM_NAME, TOTAL_READS in BAM_READS_DICT.items():
        OUTPUT_LOCATION = ARGS['<output_path>'] + "/" + BAM_NAME + "/"
        os.makedirs(OUTPUT_LOCATION, exist_ok=True)
        # Loop crossing all chromosome name
        for CHROMOSOME in CHR_LIST:
            LIST_FOR_MULTIPROCESSING.append({"chromosome": CHROMOSOME,
                                             "bam_filename": ARGS['<input_path>'] + "/" + BAM_NAME + ".bam",
                                             "output_filename": OUTPUT_LOCATION + CHROMOSOME + "_" + BAM_NAME + ".txt",
                                             "total_reads": TOTAL_READS})
    # Run the multiprocess
    with Pool(processes=NB_PROC) as pool:
        FUNC = partial(process)
        # Progress bar
        print("\n\n" + str(cpu_count()) + " cpus detected, using " + str(NB_PROC))
        print("Processing on output files (Sparse matrices) ...\n")
        tqdm(pool.map(FUNC, LIST_FOR_MULTIPROCESSING), total=len(LIST_FOR_MULTIPROCESSING))

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
