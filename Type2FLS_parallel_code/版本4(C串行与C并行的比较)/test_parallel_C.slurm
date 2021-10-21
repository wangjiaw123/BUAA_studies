#!/bin/bash
#SBATCH -J test_parallel_C
#SBATCH -p cpu-normal
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 10:00:00
#SBATCH -o test_parallel_C1.out
#SBATCH -e test_parallel_C1.err

matlab -nodesktop -nosplash -nodisplay -r "test_parallel_C_HPC6"

