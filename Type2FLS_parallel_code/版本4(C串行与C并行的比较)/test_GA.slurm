#!/bin/bash
#BATCH -J test_GA
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH -o test_GA.out
#SBATCH -e test_GA.err

matlab -nodesktop -nosplash -nodisplay -r "Intelligent_Algorithm_mainFunction_GA"

