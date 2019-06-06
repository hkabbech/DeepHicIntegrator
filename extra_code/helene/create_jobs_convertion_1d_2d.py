#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os

RESOLUTION = [100000, 50000, 25000, 10000, 5000]
CHROM = ['chr20', 'chr10', 'chr1']
HISTONE_MODIFICATION = ['HUVEC_h3k27ac', 'HUVEC_h3k27me3', 'HUVEC_h3k4me1', 'HUVEC_h3k4me3']

os.makedirs('jobs/', exist_ok=True)

for res in RESOLUTION:
    for chrom in CHROM
        for hist in HISTONE_MODIFICATION:
            test = str(chrom)+'_'+str(res)+'_'+str(hist)
            with open('jobs/job_'+test+'.sh', 'w') as filout:
                filout.write('#!/bin/bash\n')
                filout.write('#BSUB -q mpi\n')
                filout.write('#BSUB -o results/'+test+'_out-%J.txt\n')
                filout.write('#BSUB -e results/'+test+'_err-%J.txt\n\n')
                filout.write('./convert_1d_2d_histone_modification.py data/histone_modification/{}.bam -r {} -c {} -o results/{}\n'.format(hist, str(res), chrom, test))
