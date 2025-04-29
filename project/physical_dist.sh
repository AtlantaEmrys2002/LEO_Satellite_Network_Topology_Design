export PYTHONPATH="$PWD"

# Used to evaluate topologies that were constructed according to physical distance and demonstrate that measure can be
# used as a separate evaluation tool (not just with cost function coefficient optimisation)

# GENERATION

# Kuiper
python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise False --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --weights 0 1 0

# Starlink
python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise False --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --weights 0 1 0

# Telesat
python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise False --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --weights 0 1 0

# EVALUATION
python analysis/measure.py --constellation_name Kuiper-630 --num_satellites 1156 --num_snapshots 97 --topology_file novel/aco/kuiper-630 --topology_type dynamic
