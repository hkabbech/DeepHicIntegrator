
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
totalCount = float(sys.argv[4])
chromSizesFileName = sys.argv[5]
signalBamFileName = sys.argv[6]
outputLocation = sys.argv[7]

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

def writing_mark_matrix_sparse(chromosome, resolution, minimum_reads_threshold, total_count, chromosome_list, chrom_sizes_dict, signal_bam_file, output_location):

  # Resulting dictionary
  signal_distribution_dict = dict()

  # Opening files
  output_matrix_file_name = output_location + chromosome + ".txt"
  output_matrix_file = open(output_matrix_file_name, "w")

  # Iterating on genome
  for pos1 in range(0, chrom_sizes_dict[chromosome]+resolution, resolution):

    # Fetching locations
    start1 = pos1
    end1 = min(pos1+resolution, chrom_sizes_dict[chromosome])
    if(end1 <= start1): continue

    # Fetching reads
    region1 = [chromosome, start1, end1]
    try: total_reads1 = fetch_total_reads_bam(signal_bam_file, region1)
    except Exception: continue
    if(math.isnan(total_reads1) or not np.isfinite(total_reads1)): continue
    if(total_reads1 == 0 or total_reads1 < minimum_reads_threshold): total_reads1 = 1

    # Iterating on genome
    for pos2 in range(start1, chrom_sizes_dict[chromosome]+resolution, resolution):

      # Fetching locations
      start2 = pos2
      end2 = min(pos2+resolution, chrom_sizes_dict[chromosome])
      if(end2 <= start2): continue

      # Fetching reads
      region2 = [chromosome, start2, end2]
      try: total_reads2 = fetch_total_reads_bam(signal_bam_file, region2)
      except Exception: continue
      if(math.isnan(total_reads2) or not np.isfinite(total_reads2)): continue
      if(total_reads2 == 0): total_reads2 = 1
      if(total_reads1 < minimum_reads_threshold and total_reads2 < minimum_reads_threshold): continue

      # Writing to file
      total_reads = (total_reads1 * total_reads2) / total_count
      output_matrix_file.write("\t".join([chromosome, str(start1), str(start2), str(total_reads)])+"\n")
      output_matrix_file.write("\t".join([chromosome, str(start2), str(start1), str(total_reads)])+"\n")

  # Closing files
  output_matrix_file.close()
  output_matrix_file = None
  gc.collect()

  # Returning objects
  return 0

def write_mark_dictionary(chromosome, resolution, minimum_reads_threshold, total_count, chromosome_sizes_file_name, signal_bam_file_name, output_location):

  # Create output location
  command = "mkdir -p "+output_location
  os.system(command)

  # Reading chromosome sizes
  chromosome_list, chrom_sizes_dict = get_chromosome_sizes_dictionary(chromosome_sizes_file_name)

  # Opening files
  signal_bam_file = Samfile(signal_bam_file_name, "rb")

  # Writing distribution of signals
  writing_mark_matrix_sparse(chromosome, resolution, minimum_reads_threshold, total_count, chromosome_list, chrom_sizes_dict, signal_bam_file, output_location)

  # Closing files
  signal_bam_file.close()

  # Free memory from unused objects
  chromosome_list = None
  signal_bam_file = None
  chrom_sizes_dict = None
  gc.collect()

###################################################################################################
# EXECUTION
###################################################################################################

# Writing mark dictionary
write_mark_dictionary(chromosome, resolution, minimum_reads_threshold, totalCount, chromSizesFileName, signalBamFileName, outputLocation)


