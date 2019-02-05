
# Import
import os
import sys

# Folder List
counter = 1
juicerCommand = "juicertools"
kindMatrix = "observed"
kindNormalization = "KR"
unitResolution = "BP"
chromSizesFile = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Code/12_DirichletV1/input/"
itl = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Results/10_Process_All_HiC_Data/0_Input_Tables/"
il = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Results/10_Process_All_HiC_Data/2_Merged_Hic/"
ol = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Results/12_DirichletV1/2_Hic_Sparse_Matrices/"
folderList = ["2_mg/hic/fasta/", "1_4dn/hic/fasta/"]

# Folder List
for fd in folderList:

  # Hic List
  hicFileNameList = []
  tableFileName = itl + "_".join(fd.split("/")[:2]) + ".txt"
  tableFile = open(tableFileName, "rU")
  for line in tableFile: hicFileNameList.append(line.strip().split("#")[1])
  tableFile.close()

  # Hic Loop
  for i in range(0,len(hicFileNameList)):

    # Resolution List
    resList = ["5000", "10000", "25000", "50000"]
    resLabel = ["5K", "10K", "25K", "50K"]

    # Resolution Loop
    for j in range(0,len(resList)):

      # Parameters
      outName = fd.split("/")[0] + "/" + hicFileNameList[i] + "_" + resLabel[j]

      # Input
      juicerCommand = juicerCommand
      kindOfMatrix = kindMatrix
      kindOfNormalization = kindNormalization
      unitOfResolution = unitResolution
      resolution = resList[j]
      chromSizesFileName = chromSizesFile
      inputHicFileName = il + hicFileNameList[i] + "/inter_30.hic"
      outputLocation = ol + outName + "/"

      # Creating files
      inFileName = fl + str(counter) + "_chs.txt"
      inFile = open(inFileName,"w")
      inFile.write("\n".join([juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution, resolution, chromSizesFileName, inputHicFileName, outputLocation]))
      inFile.close()
      counter += 1


