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

module load python/3.9.9

echo "Loaded python module..."

bash ./../install_dependencies.sh

echo "Installed dependencies..."

chmod --recursive u+x ./../project/*

echo "Modified permissions..."

bash report.sh

echo "DONE"

# References:
# https://www.durham.ac.uk/research/institutes-and-centres/advanced-research-computing/hamilton-supercomputer/software/
# developers/python/
# https://stackoverflow.com/questions/51390968/python-ssl-certificate-verify-error
# https://stackoverflow.com/questions/76321508/importerror-cannot-import-name-csr-array-from-scipy-sparse
# https://stackoverflow.com/questions/73328773/how-do-i-make-list-of-files-executable-at-once-in-git-bash
# https://stackoverflow.com/questions/10319652/check-if-a-file-is-executable