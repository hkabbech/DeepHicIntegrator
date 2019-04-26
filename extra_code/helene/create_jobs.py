#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os

EPOCHS = 200
RESOLUTION = 25000
N = [40, 60, 80]
BATCH_SIZE = [16, 32, 64, 128, 256, 512]

os.makedirs('jobs/', exist_ok=True)

for i in [25000, 10000, 5000]:
	res_path = 'results/interpol_resolution_'+str(i)+'/'
	os.makedirs(res_path, exist_ok=True)
	with open('jobs/job_resolution_'+str(i)+'.sh', 'w') as filout:
		filout.write('#!/bin/bash\n')
		filout.write('#BSUB -q gpu\n')
		filout.write('#BSUB -W 10:00\n')
		filout.write('#BSUB -R "ngpus=1"\n')
		filout.write('#BSUB -o '+res_path+'out-%J.txt\n')
		filout.write('#BSUB -e '+res_path+'err-%J.txt\n\n')

		filout.write('module load cuda90/toolkit/9.0.176\n')
		filout.write('module load cuda90/blas/9.0.176\n')
		filout.write('module load cudnn/90v7.3.1\n\n')

		filout.write('./deep_hic_integrator data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic data/1_binaries/hic/GSE63525_HUVEC_combined_30.hic -o '+res_path+' --test 20 --train 1 --resolution '+str(i)+'\n')
