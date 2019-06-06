#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
            ./convert_1d_2d_histone_modification.py <signal_bam_file> [--chromosome STR]
                                                                      [--resolution INT]
                                                                      [--chrom_size FILE]
                                                                      [--region_bam_file FILE]
                                                                      [--signal_bam_file FILE]
                                                                      [--output FILE]

    Arguments:
        <signal_bam_file>                   Path of the bam file containing the signals.

    Options:
        -c STR, --chromosome STR            [default: chr20]
        -r INT, --resolution INT            Resolution value. [default: 25000]
        -s FILE, --chrom_size FILE          [default: data/chrom.sizes.hg19.filter]
        -b FILE, --region_bam_file FILE     [default: data/histone_modification/all_peaks.bam]
        -o FILE, --output FILE              [default: results/sparse_matrix]
        -h, --help                          Show this.
"""

# Third-party modules
from datetime import datetime
from schema import Schema, And, Use, SchemaError
from docopt import docopt
import numpy as np
from pysam import Samfile
import pandas as pd
from scipy.sparse import coo_matrix


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        '<signal_bam_file>': Use(open),
        '--chromosome': And(Use(str)),
        '--resolution': And(Use(int), lambda n: n >= 0,
                            error='--resolution=INT should be a positive integer'),
        '--chrom_size': And(Use(open)),
        '--region_bam_file': And(Use(open)),
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
                           format: column_1 \t column_2
        Returns:
            dict: The created dictionary
    """
    info_dict = {}
    with open(filename, "r") as file:
        for line in file:
            contents = line.split("\t")
            info_dict[contents[0]] = int(contents[1])
    return info_dict

def count_total_reads(bam_file, region):
    """
        Count the number of reads in the region.

        Args:
            bam_file():
            region(dict):

        Returns:
            int: The number of reads in the region
    """
    total = 0
    for _ in bam_file.fetch(region['chrom'], region['start'], region['end']):
        total += 1
    return total

def fetch_regions_in_iterval(bam_file, region):
    """
        Fetch all regions in the intervalle.

        Args:
            bam_file():
            region(dict):

        Returns:
            list: List of all regions
    """
    # Initializing region_list that will contain the regions fetched
    region_list = []
    # Fetching regions
    for read in bam_file.fetch(region['chrom'], region['start'], region['end']):
        pos_1 = max(region['start'], read.reference_start)
        pos_2 = min(region['end'], read.reference_end)
        if pos_1 < pos_2:
            new_region = {'chrom': region['chrom'], 'start': pos_1, 'end': pos_2}
            region_list.append(new_region)
    return region_list

def write_mark_matrix_sparse(chromosome, size, resolution, args, output, pseudo_count=1):
    """
        TO DO

        Args:
            chromosome(str):
            size(int):
            resolution(int):
            args():
            output(str):
            pseudo_count(int):
    """
    signal_bam_file = Samfile(args['<signal_bam_file>'], "rb")
    region_bam_file = Samfile(args['--region_bam_file'], "rb")

    lines = ""
    # Iterating on genome
    for start_1 in range(0, size+resolution, resolution):
        # Initializing read_list_1
        read_list_1 = []

        # Calculating start/end locations
        end_1 = min(start_1+resolution, size)
        region_1 = {'chrom': chromosome, 'start': start_1, 'end': end_1}

        # Fetching the regions that fall within [start_1, end_1]
        region_list_1 = fetch_regions_in_iterval(region_bam_file, region_1)

        # Fetching reads
        for region in region_list_1:
            try:
                region_size = region['end'] - region['start']
                total_reads_1 = count_total_reads(signal_bam_file, region) / region_size
            except Exception:
                continue
            if not np.isfinite(total_reads_1):
                continue
            read_list_1.append(total_reads_1)
        total_sum_reads_1 = sum(read_list_1) + pseudo_count

        # Iterating on genome
        for start_2 in range(start_1, size+resolution, resolution):

            # Initializing read_list_2
            read_list_2 = []

            # Calculating start/end locations
            end_2 = min(start_2+resolution, size)
            region_2 = {'chrom': chromosome, 'start': start_2, 'end': end_2}

            # Fetching the regions that fall within [start_2, end_2]
            region_list_2 = fetch_regions_in_iterval(region_bam_file, region_2)

            # Fetching reads
            for region in region_list_2:
                try:
                    region_size = region['end'] - region['start']
                    total_reads_2 = count_total_reads(signal_bam_file, region) / region_size
                except Exception:
                    continue
                if not np.isfinite(total_reads_2):
                    continue
                read_list_2.append(total_reads_2)
            total_sum_reads_2 = sum(read_list_2) + pseudo_count

            total_reads = total_sum_reads_1 * total_sum_reads_2
            # Upper triangle
            lines = lines+"\t".join([chromosome, str(start_1), str(start_2), str(total_reads)])+"\n"

    with open(output+'.bed', 'w') as file:
        file.write(lines)

    # Closing files
    signal_bam_file.close()
    region_bam_file.close()

def create_sparse_numpy(resolution, output):
    """
        TO DO

        Args:
            resolution(int):
            output(str):
    """
    mark_df = pd.read_csv(output+'.bed', sep='\t', header=None)
    mark_df.columns = ['chr', 'base_1', 'base_2', 'value']

    data = mark_df.value
    # base_1 and base_2 columns must be converted to index by dividing by the resolution number
    # This step is necesary for the creation of the sparse matrix with scipy.sparse
    row = ((mark_df.base_1 / resolution)).astype(int)
    col = ((mark_df.base_2 / resolution)).astype(int)
    # Creation of the sparse matrix and conversion into a numpy array
    size = int(max(max(mark_df['base_2']), max(mark_df['base_1'])) / resolution)
    matrix = coo_matrix((data, (row, col)), shape=(size+1, size+1)).toarray()
    matrix = np.float32(matrix)
    # matrix = np.log(matrix+1)
    # Rescaling of the values in range 0-1 (min-max scaling method)
    matrix = (matrix - matrix.min()) / (matrix.max() - matrix.min())
    np.save(output+'.npy', matrix)


if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################

    ARGS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()

    CHROMOSOME = ARGS['--chromosome']
    RESOLUTION = int(ARGS['--resolution'])
    CHR_SIZE_DICT = create_dictionary(ARGS['--chrom_size'])
    SIZE = CHR_SIZE_DICT[CHROMOSOME]
    OUTPUT = ARGS['--output']

    # Writing distribution of signals
    write_mark_matrix_sparse(CHROMOSOME, SIZE, RESOLUTION, ARGS, OUTPUT)
    # create_sparse_numpy(RESOLUTION, OUTPUT)

    TIME = datetime.now() - START_TIME
    HOURS = TIME.seconds // 3600
    MINUTES = TIME.seconds // 60 % 60
    print('\nTotal runtime: {}h {}min'.format(HOURS, MINUTES))
