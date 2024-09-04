#!/bin/bash
#SBATCH -p cl2
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --ntasks 1
#SBATCH --job-name KATZ

module purge
module load matlab
module list

matlab -batch "Test_Katz_Forcing" 
