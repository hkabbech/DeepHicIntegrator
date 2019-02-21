#! /bin/sh

# Sparse matrix for bam files
./1.1_convert_bam_to_sparse_matrix.py ../../data/toy_example/1_binaries/bam\
                                      ../../data/toy_example/1_binaries/bam_total-reads.txt\
                                      ../../data/toy_example/1_binaries/chrom.sizes.hg19.filter\
                                      ../../data/toy_example/2_sparse_matrices/ngs_data


