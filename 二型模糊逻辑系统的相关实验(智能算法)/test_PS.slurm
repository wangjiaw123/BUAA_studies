#!/bin/bash
#BATCH -J tsk_PS
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH -o test_PS.out
#SBATCH -e test_PS.err

matlab -nodesktop -nosplash -nodisplay -r "tIntelligent_Algorithm_mainFunction_PS"

