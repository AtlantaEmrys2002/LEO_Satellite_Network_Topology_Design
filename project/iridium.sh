echo "DEMONSTRATION EXAMPLE"
echo ""
echo "Please note that you may have to brew install proj to ensure pip install geopandas."
echo "Starting report recreation ..."

echo -n "Iridium +Grid Topology: "
python __main__.py --tles iridium-constellation_tles.txt.tmp --constellation Iridium-780 --m 34 --n 34 --i 51.9 --rev 14.8 --multi False --optimise False --topology plus-grid --isl_terminals 4 --snapshot_interval 60 --weights 1 1 1

