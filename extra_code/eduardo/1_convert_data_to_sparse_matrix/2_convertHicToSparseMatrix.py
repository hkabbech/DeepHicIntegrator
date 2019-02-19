
# Import
import os
import gc
import sys
import math
sys.path = ["/home/egusmao/.local/lib/python2.7/site-packages"] + sys.path
import numpy as np

###################################################################################################
# INPUT
###################################################################################################

# Input
juicerCommand = sys.argv[1] # "juicertools"
kindOfMatrix = sys.argv[2] # "observed"
kindOfNormalization = sys.argv[3] # "KR"
unitOfResolution = sys.argv[4] # "BP"
resolution = sys.argv[5] # "25000"
chromSizesFileName = sys.argv[6] # "/home/egg/Projects/Papantonis_Intrinsic/Tests/pe_scan/chrom.sizes.hg19.filter"
inputHicFileName = sys.argv[7] # "/home/egg/Projects/Papantonis_Intrinsic/Tests/pe_scan/data/GSE63525_NHEK_combined_30.hic"
outputLocation = sys.argv[8] # "/home/egg/Projects/Papantonis_Intrinsic/Tests/pe_scan/data/NHEK/"

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

def write_hic_dictionary(juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chromosome_list, chrom_sizes_dict, input_hic_file_name, output_location):

  # Creating individual sparse chromosome matrices
  for chrom in chromosome_list:

    # Initialization
    chrWoChr = chrom.split("chr")[-1]
    region = ":".join([chrWoChr,"1",str(chrom_sizes_dict[chrom])])

    # Creating sparse matrix
    tempFileName = output_location + chrom + "_temp.txt"
    command = " ".join([juicer_command, "dump", kind_of_matrix, kind_of_normalization, input_hic_file_name,
                      region, region, unit_of_resolution, resolution, tempFileName])
    os.system(command)

    # Writing entries
    tempFile = open(tempFileName, "rU")
    outputFileName = output_location + chrom + ".txt"
    outFile = open(outputFileName, "w", buffering=1)
    for line in tempFile:
      ll = line.strip().split("\t")
      value = float(ll[2])
      if(math.isnan(value) or not np.isfinite(value)): continue
      outFile.write("\t".join([chrom]+ll[:2]+[str(value)])+"\n")
      outFile.flush()
    tempFile.close()
    outFile.close()

    # Removing temporary file
    command = "rm "+tempFileName
    os.system(command)

    # Cleaning objects
    region = None
    tempFile = None
    outFile = None
    gc.collect()    

def create_hic_dictionary(juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chromosome_sizes_file_name, input_hic_file_name, output_location):

  # Create output location
  command = "mkdir -p "+output_location
  os.system(command)

  # Reading chromosome sizes
  chromosome_list, chrom_sizes_dict = get_chromosome_sizes_dictionary(chromosome_sizes_file_name)

  # Creating individual sparse chromosome matrices by chromosome
  write_hic_dictionary(juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chromosome_list, chrom_sizes_dict, input_hic_file_name, output_location)

  # Cleaning structures
  chromosome_list = None
  chrom_sizes_dict = None
  gc.collect()

###################################################################################################
# EXECUTION
###################################################################################################

# Writing mark dictionary
create_hic_dictionary(juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution, resolution, chromSizesFileName, inputHicFileName, outputLocation)


