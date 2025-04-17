# RECREATES REPORT
# Runs commands that recreate (approximately, as some functions are stochastic) the results of the final report. By
# running commands separately, checkpoints are effectively built into the code. Topologies are built for each satellite
# network (Kuiper, Starlink, and Telesat)

echo "REPORT RECREATION"
echo ""
echo "Please note that you may have to brew install proj to ensure pip install geopandas."
echo "Starting report recreation ..."

# BENCHMARK TOPOLOGIES #

echo "4 ISLS/SATELLITE"

echo "Generating benchmark topologies..."
echo "Please note that +Grid topologies are generated with the assumption that each satellite has 4 active ISL terminals."
echo ""

## +Grid Topology Design
#
## Kuiper
#
#echo -n "Kuiper +Grid Topology: "
#python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology plus-grid --isl_terminals 4 --snapshot_interval 60
#
## Starlink
#
#echo -n "Starlink +Grid Topology: "
#python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology plus-grid --isl_terminals 4 --snapshot_interval 60
#
## Telesat
#
#echo -n "Telesat +Grid Topology: "
#python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology plus-grid --isl_terminals 4 --snapshot_interval 60

# xGrid Topology Design

# Kuiper

#echo -n "Kuiper xGrid Topology: "
#python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology x-grid --isl_terminals 4 --snapshot_interval 60
#
## Starlink
#
#echo -n "Starlink xGrid Topology: "
#python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology x-grid --isl_terminals 4 --snapshot_interval 60
#
## Telesat
#
#echo -n "Telesat xGrid Topology: "
#python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology x-grid --isl_terminals 4 --snapshot_interval 60

## Minimum Delay Topology Design Algorithm

# Kuiper

#echo -n "Kuiper MDTD Topology: "
#python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology mdtd --isl_terminals 4 --snapshot_interval 60
#
## Starlink
#
#echo -n "Starlink MDTD Topology: "
#python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology mdtd --isl_terminals 4 --snapshot_interval 60
#
## Telesat
#
#echo -n "Telesat MDTD Topology: "
#python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology mdtd --isl_terminals 4 --snapshot_interval 60

## PROPOSED ALGORITHM TOPOLOGIES #
#
## Unlike the benchmark topologies (which are generated and evaluated), the proposed algorithm also undergoes cost
## function optimisation, as well as topology builds and evaluation (according to given metrics)
#
#echo "Generating topologies with novel algorithm..."
#echo ""
#
## Kuiper
#
## Random Search Optimisation
#
#echo -n "Kuiper Novel Algorithm Topology (ACO with Random Search): "
#python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method random

##echo -n "Kuiper Novel Algorithm Topology (GA with Random Search): "
##python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method random
#
#echo -n "Kuiper Novel Algorithm Topology (Primal with Random Search): "
#python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method random

## Evol Strategy Search Optimisation
#
#echo -n "Kuiper Novel Algorithm Topology (ACO with Evolutionary Strategy Search): "
#python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method evolutionary

##echo -n "Kuiper Novel Algorithm Topology (GA with Evolutionary Strategy Search): "
##python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method evolutionary
##
echo -n "Kuiper Novel Algorithm Topology (Primal with Evolutionary Strategy Search): "
python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method evolutionary

## Starlink
#
#echo -n "Starlink Novel Algorithm Topology: "
#
## Random Search Optimisation
#
#echo -n "Starlink Novel Algorithm Topology (ACO with Random Search): "
#python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method random
#
##echo -n "Starlink Novel Algorithm Topology (GA with Random Search): "
##python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method random
##
##echo -n "Starlink Novel Algorithm Topology (Primal with Random Search): "
##python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method random
#
## Evolutionary Strategy Search Optimisation
#
##echo -n "Starlink Novel Algorithm Topology (ACO with Evolutionary Strategy Search): "
##python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method evolutionary
#
##echo -n "Starlink Novel Algorithm Topology (GA with Evolutionary Strategy Search): "
##python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method evolutionary
##
echo -n "Starlink Novel Algorithm Topology (Primal with Evolutionary Strategy Search): "
python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method evolutionary

## Telesat
#
#echo -n "Telesat Novel Algorithm Topology: "
#
## Random Search Optimisation
#
#echo -n "Telesat Novel Algorithm Topology (ACO with Random Search): "
#python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method random
##
##echo -n "Telesat Novel Algorithm Topology (GA with Random Search): "
##python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method random
##
##echo -n "Telesat Novel Algorithm Topology (Primal with Random Search): "
##python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method random
#
## Evolutionary Strategy Search Optimisation
#
##echo -n "Telesat Novel Algorithm Topology (ACO with Evolutionary Strategy Search): "
##python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method evolutionary
#
##echo -n "Telesat Novel Algorithm Topology (GA with Evolutionary Strategy Search): "
##python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method evolutionary
##
echo -n "Telesat Novel Algorithm Topology (Primal with Evolutionary Strategy Search): "
python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method evolutionary

#echo""
#echo "Report Recreation Finished"
#
#### MACHINE LEARNING ###
#
## Machine Learning Search Optimisation
#
##python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method machine-learning
##python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method machine-learning
##python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method machine-learning
#
## Machine Learning Search Optimisation
#
##echo -n "Kuiper Novel Algorithm Topology (ACO with Evolutionary Strategy Search): "
##python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method machine-learning
##
##echo -n "Kuiper Novel Algorithm Topology (GA with Evolutionary Strategy Search): "
##python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method machine-learning
##
##echo -n "Kuiper Novel Algorithm Topology (Primal with Evolutionary Strategy Search): "
##python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method machine-learning
#
## Machine Learning Search Optimisation
#
##python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst aco --optimisation_method machine-learning
##python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst ga --optimisation_method machine-learning
##python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology novel --isl_terminals 4 --snapshot_interval 60 --dcmst primal --optimisation_method machine-learning
#
