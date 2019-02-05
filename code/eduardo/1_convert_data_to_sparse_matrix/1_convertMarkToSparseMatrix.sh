#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=MTS
#SBATCH --output=MTS.%A_%a.out
#SBATCH --error=MTS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=10:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-1383

# Commands
# sbatch 1_convertMarkToSparseMatrix.sh
# squeue -u egusmao
# scancel 10831898, 10833884

# Modules
module add python/2.7.5-V2
module add intel/17.0_gnu_5.1.0
module add bowtie2/2.2.9
module add samtools/1.6
module add R/3.3.3_intel_mkl
module add cuda/8.0
export PATH=$PATH:"/home/egusmao/.local/bin"
export PATH=$PATH:"/projects/ag-papan/install/FastQC/"
export PATH=$PATH:"/projects/ag-papan/install/cutadapt-1.15/inst/bin/"
export PATH=$PATH:"/projects/ag-papan/install/TrimGalore-0.4.3/"
export PATH=$PATH:"/projects/ag-papan/install/sratoolkit.2.8.2-1-ubuntu64/bin/"
export PATH=$PATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/bin/"
export PATH=$PATH:"/projects/ag-papan/install/bedtools-2.25.0/bin/"
export PATH=$PATH:"/projects/ag-papan/install/juicertools-1.7.6/bin/"
export PATH=$PATH:"/projects/ag-papan/install/chromhmm-1.15/bin/"
export PATH=$PATH:"/home/egusmao/software/ibm-cbc-genomic-tools-master/bin/"
export PYTHONPATH=$PYTHONPATH:"/home/egusmao/.local/lib/python2.7/site-packages"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/cutadapt-1.15/inst/lib/python2.7/site-packages/"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/lib/python2.7/site-packages/"

# Input
inputFileName="/projects/ag-papan/eduardo/Papantonis_Intrinsic/Code/12_DirichletV1/input2/"${SLURM_ARRAY_TASK_ID}"_mts.txt"
chromosome=`sed '1q;d' $inputFileName`
resolution=`sed '2q;d' $inputFileName`
minimum_reads_threshold=`sed '3q;d' $inputFileName`
chromSizesFileName=`sed '4q;d' $inputFileName`
signalBamFileName=`sed '5q;d' $inputFileName`
outputLocation=`sed '6q;d' $inputFileName`

# Creating matrix
python /projects/ag-papan/eduardo/Papantonis_Intrinsic/Code/12_DirichletV1/1_convertMarkToSparseMatrix.py $chromosome $resolution $minimum_reads_threshold $chromSizesFileName $signalBamFileName $outputLocation


