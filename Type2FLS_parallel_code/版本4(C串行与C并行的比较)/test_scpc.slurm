#!/bin/bash
#BATCH -J test_scpc
#SBATCH -p cpu-low
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH -o test_scpc.out
#SBATCH -e test_scpc.err

matlab -nodesktop -nosplash -nodisplay -r "test_parallel_HPC_7_scpc"

