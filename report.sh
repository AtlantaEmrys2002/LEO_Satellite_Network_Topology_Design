# RECREATES REPORT
# Runs commands that recreate (approximately, as some functions are stochastic) the results of the final report. By
# running commands separately, checkpoints are effectively built into the code. Topologies are built for each satellite
# network (Kuiper, Starlink, and Telesat)

# NB TIME EVERYTHING!!!!!!!!!!!!

echo "REPORT RECREATION"
echo ""
echo "Starting report recreation ..."

# BENCHMARK TOPOLOGIES #

echo "Generating benchmark topologies..."
echo "Please note that +Grid topologies are generated with the assumption that each satellite has 4 active ISL terminals."
echo ""

# +Grid Topology Design

# Kuiper

echo -n "Kuiper +Grid Topology: "
python __main__.py --tles kuiper-constellation_tles.txt.tmp --constellation Kuiper-630 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise True --topology plus-grid --isl_terminals 4 --snapshot_interval 60

# Starlink

echo -n "Starlink +Grid Topology: "
python __main__.py --tles starlink-constellation_tles.txt.tmp --constellation Starlink-550 --m 72 --n 22 --i 53 --rev 15.9 --multi False --optimise True --topology plus-grid --isl_terminals 4 --snapshot_interval 60

# Telesat

echo -n "Telesat +Grid Topology: "
python __main__.py --tles telesat-constellation_tles.txt.tmp --constellation Telesat-1015 --m 27 --n 13 --i 98.98 --rev 13.66 --multi False --optimise True --topology plus-grid --isl_terminals 4 --snapshot_interval 60

echo ""

# Minimum Delay Topology Design Algorithm

# Kuiper

echo -n "Kuiper MDTD Topology: "
# BUILD!


# Starlink

echo -n "Starlink MDTD Topology: "
# BUILD !


# Telesat

echo -n "Telesat MDTD Topology: "
# BUILD!


# PROPOSED ALGORITHM TOPOLOGIES #

# Unlike the benchmark topologies (which are generated and evaluated), the proposed algorithm also undergoes cost
# function optimisation, as well as topology builds and evaluation (according to given metrics)

echo "Generating topologies with novel algorithm..."
echo ""

echo -n "Kuiper Novel Algorithm Topology: "
# BUILD!


# Starlink

echo -n "Starlink Novel Algorithm Topology: "
# BUILD !


# Telesat

echo -n "Telesat Novel Algorithm Topology: "
# BUILD!


# SUMMARY OF RESULTS #
# Write python for this !!!!!!!!

echo""
echo "Report Recreation Finished"
