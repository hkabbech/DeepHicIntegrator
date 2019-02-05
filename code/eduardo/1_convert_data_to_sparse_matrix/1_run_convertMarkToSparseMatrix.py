
# Import
import os
import sys

# Folder List
counter = 1
flag = False
chromSizesFile = "/home/egusmao/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Code/12_DirichletV1/input2/"
itl = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Results/9_Process_All_NGS_Data/0_Input_Tables/"
il = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Results/9_Process_All_NGS_Data/2_Merged_Bam_Files/"
ol = "/projects/ag-papan/eduardo/Papantonis_Intrinsic/Results/12_DirichletV1/1_Signal_Sparse_Matrices/"
folderList = ["1_4dn_atac", "2_encode_chipseq", "2_encode_dnaseseq", "3_roadmap_fetal_heart_chipseq", "3_roadmap_fetal_heart_dnaseseq", "6_ang2016_atac", "6_ang2016_chipseq", "7_banovich2018_atac", "8_sakabe2018_chipseq"]

# Folder Loop
for fd in folderList:

  # Bam List
  bamFileNameList = []
  tableFileName = itl + fd + ".txt"
  tableFile = open(tableFileName, "rU")
  for line in tableFile: bamFileNameList.append(line.strip().split("#")[1])
  tableFile.close()

  # Bam Loop
  for bamFileName in bamFileNameList:

    # Chromosome List
    chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

    # Chromosome Loop
    for chrom in chrList:

      if(counter > 1400 or flag == True):
        flag = True
        # Input
        chromosome = chrom
        resolution = "1000"
        minimum_reads_threshold = "10"
        chromSizesFileName = chromSizesFile
        signalBamFileName = il + fd + "/" + bamFileName + ".bam"
        outputLocation = ol + fd + "/" + bamFileName + "/"

        # Creating files
        inFileName = fl + str(counter-1400) + "_mts.txt"
        inFile = open(inFileName,"w")
        inFile.write("\n".join([chromosome, resolution, minimum_reads_threshold, chromSizesFileName, signalBamFileName, outputLocation]))
        inFile.close()
      counter += 1


