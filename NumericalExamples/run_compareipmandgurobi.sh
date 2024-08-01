#!/bin/bash
#SBATCH -p gpu
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=256
#SBATCH --ntasks 1
#SBATCH --nodelist=gpu02
#SBATCH --job-name CIPMGUR

module purge
module load matlab
module list

matlab -batch "compare_ipm_gurobi" 
