#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os

N = [20, 40, 60, 80, 120]
BATCH_SIZE = [16, 32, 64, 128, 256]


os.makedirs('jobs/', exist_ok=True)

for n in N:
    for b in BATCH_SIZE:
        test = '04_29_2019_3Conv_1MaxPool_1Upsampl_N'+str(n)+'_B'+str(b)
        result_path = 'results/'+test+'/'
        os.makedirs(result_path, exist_ok=True)
        with open('jobs/job_'+test+'.sh', 'w') as filout:
            filout.write('#!/bin/bash\n')
            filout.write('#BSUB -q gpu\n')
            filout.write('#BSUB -W 10:00\n')
            filout.write('#BSUB -o '+result_path+'out-%J.txt\n')
            filout.write('#BSUB -e '+result_path+'err-%J.txt\n\n')

            filout.write('module load cuda90/toolkit/9.0.176\n')
            filout.write('module load cuda90/blas/9.0.176\n')
            filout.write('module load cudnn/90v7.3.1\n\n')

            filout.write('./deep_hic_integrator data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic -o {} -b {} -n {}\n'.format(result_path, b, n))
