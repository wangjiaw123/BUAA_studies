#!/bin/bash
#BATCH -J test_FLS
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 27
#SBATCH -t 72:00:00
#SBATCH -o test_FLS_HPC_MC_tsk.out
#SBATCH -e test_FLS_HPC_MC_tsk.err

matlab -nodesktop -nosplash -nodisplay -r "test_parallel_TSK_MC_HPC"

