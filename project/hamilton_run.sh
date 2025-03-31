#!/bin/bash

# REQUEST RESOURCES:
# Allocate 1 core
#SBATCH -c 90
#SBATCH --mem=5G
#SBATCH --time=2-0:50:0
#SBATCH --gres=tmp:5G
# Run in the 'shared' queue
#SBATCH -p shared

# Commands to be run:
module load python/3.12.5
bash report.sh
