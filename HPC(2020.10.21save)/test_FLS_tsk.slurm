#!/bin/bash
#BATCH -J tsk_FLS
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH -o test_FLS_HPC_tsk1.out
#SBATCH -e test_FLS_HPC_tsk1.err

matlab -nodesktop -nosplash -nodisplay -r "test_parallel_TSK_HPC"

