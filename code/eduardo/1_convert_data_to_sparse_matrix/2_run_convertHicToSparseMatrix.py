
# Import
import os
import sys

# Hic List
juicerCommand = "juicertools"
kindMatrix = "observed"
kindNormalization = "KR"
unitResolution = "BP"
chromSizesFile = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/projects/ag-papan/eduardo/Papantonis_Integrative/code/eduardo/1_convert_data_to_sparse_matrix/"
il = "/projects/ag-papan/eduardo/Papantonis_Integrative/data/rao_juicer_hic/"
ol = "/projects/ag-papan/eduardo/Papantonis_Integrative/results/eduardo/1_convert_data_to_sparse_matrix/2_hic_matrices/"
hicList = ["GSE63525_HUVEC"]

# Opening input matrix file
inFileName = fl + "2_chs.txt"
inFile = open(inFileName,"w")

# Hic Loop
for hicName in range(0,len(hicList)):

  # Resolution List
  resList = ["25000"]
  resLabel = ["25K"]

  # Resolution Loop
  for j in range(0,len(resList)):

    # Parameters
    outName = hicName + "_" + resLabel[j]

    # Input
    juicerCommand = juicerCommand
    kindOfMatrix = kindMatrix
    kindOfNormalization = kindNormalization
    unitOfResolution = unitResolution
    resolution = resList[j]
    chromSizesFileName = chromSizesFile
    inputHicFileName = il + hicName + "_combined_30.hic"
    outputLocation = ol + outName + "/"

    # Write to input matrix
    inFile.write("\n".join([juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution, resolution, chromSizesFileName, inputHicFileName, outputLocation]))

# Closing input matrix file
inFile.close()


