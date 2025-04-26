# Toy example used in project demonstration exam - Iridium is relatively small compared with constellations such as Kuiper,
# Starlink, and Lightspeed (Telesat), with only 66 satellites distributed across 11 near-polar orbits.

echo "DEMONSTRATION EXAMPLE"
echo ""
echo "Please note that you may have to brew install proj to ensure pip install geopandas."
echo "Starting report recreation ..."

echo -n "Iridium +Grid Topology: "
python __main__.py --tles iridium-constellation_tles.txt.tmp --constellation Iridium-780 --m 6 --n 11 --i 86.4 --rev 14.8 --multi False --optimise False --topology plus-grid --isl_terminals 4 --snapshot_interval 60 --weights 1 1 1

echo -n "Visualising..."
python analysis/visualisation.py --tles iridium-constellation_tles.txt.tmp --location ./plus_grid/iridium-780 --constellation_name Iridium-780 --snapshot_interval 60 --method plus_grid --topology-type static --num_snapshots 96

echo -n "Iridium MDTD Topology: "
python __main__.py --tles iridium-constellation_tles.txt.tmp --constellation Iridium-780 --m 6 --n 11 --i 86.4 --rev 14.8 --multi False --optimise False --topology mdtd --isl_terminals 4 --snapshot_interval 60 --weights 1 1 1

echo -n "Visualising..."
python analysis/visualisation.py --tles iridium-constellation_tles.txt.tmp --location ./Results/mdtd/iridium-780 --constellation_name Iridium-780 --snapshot_interval 60 --method mdtd --topology-type dynamic --num_snapshots 96

echo -n "Iridium Novel Algorithm Topology (Primal DCMST): "
python __main__.py --tles iridium-constellation_tles.txt.tmp --constellation Iridium-780 --m 6 --n 11 --i 86.4 --rev 14.8 --multi False --optimise False --topology novel --dcmst "primal" --isl_terminals 4 --snapshot_interval 60 --weights 1 1 1

echo -n "Visualising..."
python analysis/visualisation.py --tles iridium-constellation_tles.txt.tmp --location ./novel/primal/iridium-780 --constellation_name Iridium-780 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 96

echo""
echo "EXAMPLE TOPOLOGIES CONSTRUCTED"