export PYTHONPATH="$PWD"
# USED TO EVALUATE FAULT TOLERANCE

# GENERATE TOPOLOGIES (HIGH GAMMA COEFFICIENT)

# Kuiper
#python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise False --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --weights 0 1 0.5

## Starlink
#python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise False --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --weights 0 1 0.5
#
## Telesat
#python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise False --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --weights 0 1 0.5

# EVALUATE FAULT TOLERANCE
python analysis/satellite_failures.py --num_snapshots 97 --topology_location novel/primal/kuiper-630 --num_sat 1156 --constellation_name Kuiper-630 --prob_failure 0.2

#python analysis/satellite_failures.py --num_snapshots 90 --topology_location novel/primal/starlink-550 --num_sat 1584 --constellation_name Starlink-550
#
#python analysis/satellite_failures.py --num_snapshots 105 --topology_location novel/primal/telesat-1015 --num_sat 351 --constellation_name Telesat-1015
