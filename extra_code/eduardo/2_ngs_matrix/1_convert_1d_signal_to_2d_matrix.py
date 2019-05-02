
# Import
import os
import gc
import sys
import math
import numpy as np
from pysam import Samfile


# Input
chromosome = sys.argv[1]
resolution = int(sys.argv[2])
chromosome_sizes_file_name = sys.argv[3]
region_file_name = sys.argv[4]
signal_bam_file_name = sys.argv[5]
output_matrix_file_name = sys.argv[6]

###################################################################################################
# AUXILIARY FUNCTIONS
###################################################################################################

def get_chromosome_sizes_dictionary(chromosome_sizes_file_name):

  # Read chromosome sizes
  chrom_sizes_dict = dict()
  chrom_sizes_file = open(chromosome_sizes_file_name,"r")
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
  try:
    for read in bam_file.fetch(region[0], region[1], region[2]): returnN += 1
  except Exception: pass

  # Returning objects
  return returnN

def fetch_regions_in_iterval(bam_file, region):

  # Initializing region_list that will contain the regions fetched
  region_list = []

  # Fetching regions
  try:
    for read in bam_file.fetch(region[0], region[1], region[2]):
      pos1 = max(region[1], read.reference_start)
      pos2 = min(region[2], read.reference_end)
      if(pos1 < pos2): region_list.append([region[0],pos1,pos2])
  except Exception: pass

  # Returning objects
  return region_list

def writing_mark_matrix_sparse(chromosome, resolution, chromosome_list, chrom_sizes_dict, region_file, signal_bam_file, output_matrix_file_name):

  # Opening output file
  output_matrix_file = open(output_matrix_file_name, "w")

  # Iterating on genome
  for pos1 in range(0, chrom_sizes_dict[chromosome]+resolution, resolution):
    # Initializing read_list1
    read_list1 = []

    # Calculating start/end locations
    start1 = pos1
    end1 = min(pos1+resolution, chrom_sizes_dict[chromosome])
    if(end1 <= start1): continue
    region1 = [chromosome, start1, end1]

    # Fetching the regions that fall within [start1, end1]
    region_list1 = fetch_regions_in_iterval(region_file, region1)

    # Fetching reads
    for region in region_list1:
      try: total_reads1 = fetch_total_reads_bam(signal_bam_file, region)
      except Exception: continue
      if(math.isnan(total_reads1) or not np.isfinite(total_reads1)): continue
      read_list1.append(total_reads1)

    # Iterating on genome
    for pos2 in range(start1, chrom_sizes_dict[chromosome]+resolution, resolution):

      # Initializing read_list2
      read_list2 = []

      # Calculating start/end locations
      start2 = pos2
      end2 = min(pos2+resolution, chrom_sizes_dict[chromosome])
      if(end2 <= start2): continue
      region2 = [chromosome, start2, end2]

      # Fetching the regions that fall within [start1, end1]
      region_list2 = fetch_regions_in_iterval(region_file, region2)

      # Fetching reads
      for region in region_list2:
        try: total_reads2 = fetch_total_reads_bam(signal_bam_file, region)
        except Exception: continue
        if(math.isnan(total_reads2) or not np.isfinite(total_reads2)): continue
        read_list2.append(total_reads2)

      # Writing to file
      if read_list1 or read_list2:
        total_reads = np.average(read_list1 + read_list2)
      else:
        continue
      if total_reads == 0 or math.isnan(total_reads) or not np.isfinite(total_reads):
        continue
      output_matrix_file.write("\t".join([chromosome, str(start1), str(start2), str(total_reads)])+"\n") # Upper triangle
      output_matrix_file.write("\t".join([chromosome, str(start2), str(start1), str(total_reads)])+"\n") # Lower triangle

  # Closing files
  output_matrix_file.close()
  output_matrix_file = None
  gc.collect()

  # Returning objects
  return 0

###################################################################################################
# MAIN FUNCTION
###################################################################################################

def convert_1d_to_2d(chromosome, resolution, chromosome_sizes_file_name, region_file_name, signal_bam_file_name, output_matrix_file_name):

  # Reading chromosome sizes
  chromosome_list, chrom_sizes_dict = get_chromosome_sizes_dictionary(chromosome_sizes_file_name)

  # Opening bam files
  region_file = Samfile(region_file_name, "rb")
  signal_bam_file = Samfile(signal_bam_file_name, "rb")

  # Writing distribution of signals
  writing_mark_matrix_sparse(chromosome, resolution, chromosome_list, chrom_sizes_dict, region_file, signal_bam_file, output_matrix_file_name)

  # Closing files
  region_file.close()
  signal_bam_file.close()

  # Free memory from unused objects
  chromosome_list = None
  signal_bam_file = None
  chrom_sizes_dict = None
  gc.collect()


convert_1d_to_2d(chromosome, resolution, chromosome_sizes_file_name, region_file_name, signal_bam_file_name, output_matrix_file_name)