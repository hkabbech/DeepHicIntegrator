
# Import
import os
import gc
import sys
import math
sys.path = ["/home/egusmao/.local/lib/python2.7/site-packages"] + sys.path
import numpy as np
from scipy import stats
from pysam import Samfile

###################################################################################################
# INPUT
###################################################################################################

# Input
chromosome = sys.argv[1]
resolution = int(sys.argv[2])
minimum_reads_threshold = float(sys.argv[3])
chromSizesFileName = sys.argv[4]
signalBamFileName = sys.argv[5]
outputLocation = sys.argv[6]

# Initialization
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# FUNCTIONS
###################################################################################################

def get_chromosome_sizes_dictionary(chromosome_sizes_file_name):

  # Read chromosome sizes
  chrom_sizes_dict = dict()
  chrom_sizes_file = open(chromosome_sizes_file_name,"rU")
  for line in chrom_sizes_file:
    ll = line.strip().split("\t")
    chrom_sizes_dict[ll[0]] = int(ll[1])
  chrom_sizes_file.close()
  chromosome_list = sorted(chrom_sizes_dict.keys())

  # Returning objects
  return chromosome_list, chrom_sizes_dict

def fetch_total_reads_bam(bam_file, region):

  # Reding bam file
  returnN = 0
  for read in bam_file.fetch(region[0], region[1], region[2]): returnN += 1

  # Returning objects
  return returnN

def distribution_of_signal_dictionary(chromosome, resolution, minimum_reads_threshold, chromosome_list, chrom_sizes_dict, signal_bam_file):

  # Resulting list
  signal_distribution_dict = dict()
  signal_distribution_list = []

  # Iterating on chromosomes
  for chrom in chromosome_list:

    if(chrom != chromosome): continue

    # Iterating on genome
    for pos in range(0, chrom_sizes_dict[chrom]+resolution, resolution):

      # Fetching locations
      start = pos
      end = min(pos+resolution, chrom_sizes_dict[chrom])
      if(end <= start): continue

      # Fetching reads
      region = [chrom, start, end]
      try: total_reads = fetch_total_reads_bam(signal_bam_file, region)
      except Exception: continue

      # Writing to dictionary
      if(total_reads <= minimum_reads_threshold or math.isnan(total_reads) or not np.isfinite(total_reads)): continue
      signal_distribution_dict[":".join([str(e) for e in region])] = total_reads
      signal_distribution_list.append(total_reads)

  # Returning objects
  return signal_distribution_list, signal_distribution_dict

def write_mark_dictionary(chromosome, resolution, minimum_reads_threshold, chromosome_sizes_file_name, signal_bam_file_name, output_location):

  # Create output location
  command = "mkdir -p "+output_location
  os.system(command)

  # Reading chromosome sizes
  chromosome_list, chrom_sizes_dict = get_chromosome_sizes_dictionary(chromosome_sizes_file_name)

  # Opening files
  signal_bam_file = Samfile(signal_bam_file_name, "rb")

  # Reading distribution of signals
  signal_distribution_list, signal_distribution_dict = distribution_of_signal_dictionary(chromosome, resolution, minimum_reads_threshold, chromosome_list, chrom_sizes_dict, signal_bam_file)

  # Iterating on chromosomes
  for chrom in chromosome_list:

    if(chrom != chromosome): continue

    # Opening files
    output_matrix_file_name = output_location + chrom + ".txt"
    output_matrix_file = open(output_matrix_file_name, "w", buffering=1)

    # Iterating on genome
    for pos in range(0, chrom_sizes_dict[chrom]+resolution, resolution):

      # Fetching locations
      start = pos
      end = min(pos+resolution, chrom_sizes_dict[chrom])
      if(end <= start): continue

      # Fetching reads
      region = [chrom, start, end]
      try: total_reads = signal_distribution_dict[":".join([str(e) for e in region])]
      except Exception: continue

      # Calculating percentile
      percentile = round(stats.percentileofscore(signal_distribution_list, total_reads),2)

      # Writing to dictionary
      output_matrix_file.write("\t".join([str(e) for e in region + [total_reads, percentile]])+"\n")
      output_matrix_file.flush()

      start = None
      end = None
      region = None
      total_reads = None
      percentile = None
      gc.collect()

    # Closing files
    output_matrix_file.close()
    output_matrix_file = None
    gc.collect()

  # Closing files
  signal_bam_file.close()

  # Free memory from unused objects
  signal_bam_file = None
  signal_distribution = None
  chromosome_list = None
  chrom_sizes_dict = None
  gc.collect()

###################################################################################################
# EXECUTION
###################################################################################################

# Writing mark dictionary
write_mark_dictionary(chromosome, resolution, minimum_reads_threshold, chromSizesFileName, signalBamFileName, outputLocation)


