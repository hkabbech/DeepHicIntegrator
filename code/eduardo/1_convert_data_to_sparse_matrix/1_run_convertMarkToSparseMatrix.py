
# Import
import os
import sys

# Folder List
chromSizesFile = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/projects/ag-papan/eduardo/Papantonis_Integrative/code/eduardo/1_convert_data_to_sparse_matrix/"
il = "/projects/ag-papan/eduardo/Papantonis_Integrative/data/histone_modification/"
ol = "/projects/ag-papan/eduardo/Papantonis_Integrative/results/eduardo/1_convert_data_to_sparse_matrix/1_regulatory_matrices/"
bamFileNameList = ["HUVEC_H3K4me1", "HUVEC_H3K4me3", "HUVEC_H3K27ac", "HUVEC_H3K27me3"]
totalReadsList = ["32422952", "13940804", "16902517", "16516826"]

# Opening input matrix file
inFileName = fl + "1_mts.txt"
inFile = open(inFileName,"w")

# Bam Loop
for i in range(0,len(bamFileNameList)):

  # Parameters
  bamFileName = bamFileNameList[i]
  counts = totalReadsList[i]

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Input
    chromosome = chrom
    resolution = "25000"
    minimum_reads_threshold = "25000"
    totalCounts = counts
    chromSizesFileName = chromSizesFile
    signalBamFileName = il + bamFileName + ".bam"
    outputLocation = ol + bamFileName + "/"

    # Write to input matrix
    inFile.write(" ".join([chromosome, resolution, minimum_reads_threshold, totalCounts, chromSizesFileName, signalBamFileName, outputLocation])+"\n")
    counter += 1

# Closing input matrix file
inFile.close()


