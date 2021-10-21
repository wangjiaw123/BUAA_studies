#!/bin/bash
#SBATCH -J fig1
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 36
#SBATCH -t 4000:00
#SBATCH -o fig1.out
#SBATCH -e fig1.err

matlab -nodesktop -nosplash -nodisplay -r "fig"

