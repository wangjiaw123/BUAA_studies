#!/bin/bash
#BATCH -J test_serial_HPC
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 72:00:00
#SBATCH -o test_PS.out
#SBATCH -e test_PS.err

matlab -nodesktop -nosplash -nodisplay -r "test_serial_HPC"

