#!/bin/bash

# REQUEST RESOURCES:
# Allocate 90 cores (approximately, the number of network snapshots for which to build a topology)
#SBATCH -c 90
#SBATCH --mem=5G
# Maximum run time of 2 days
#SBATCH --time=2-0:50:0
#SBATCH --gres=tmp:5G
# Run in the 'shared' queue
#SBATCH -p shared

# Change executable permissions of files
chmod --recursive u+x ../*

# Commands to be run:
module load python/3.12.5
bash ./../install_dependencies.sh
bash report.sh

# References:
# https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/software/
# developers/python/
