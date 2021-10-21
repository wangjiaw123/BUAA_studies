#!/bin/bash
#SBATCH -J test_parallel_C
#SBATCH -p cpu-normal
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH -o test_parallel_C.out
#SBATCH -e test_parallel_C.err

matlab -nodesktop -nosplash -nodisplay -r "test_parallel_C_HPC6"

