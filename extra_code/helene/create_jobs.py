#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os

os.makedirs('jobs/', exist_ok=True)

for i in [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000]:
	os.makedirs('results/resolution_'+str(i), exist_ok=True)
	with open('jobs/job_resolution_'+str(i)+'.sh', 'w') as filout:
		filout.write('#!/bin/bash\n')
		filout.write('#BSUB -q gpu\n')
		filout.write('#BSUB -W 10:00\n')
		filout.write('#BSUB -R "ngpus=3"\n')
		filout.write('#BSUB -o results/resolution_'+str(i)+'/out-%J.txt\n')
		filout.write('#BSUB -e results/resolution_'+str(i)+'/err-%J.txt\n\n')

		filout.write('module load cuda90/toolkit/9.0.176\n')
		filout.write('module load cuda90/blas/9.0.176\n')
		filout.write('module load cudnn/90v7.3.1\n\n')

		filout.write('./deep_hic_integrator data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic -o results/resolution_'+str(i)+'/ --test 20 --train 1 --resolution '+str(i)+'\n')
