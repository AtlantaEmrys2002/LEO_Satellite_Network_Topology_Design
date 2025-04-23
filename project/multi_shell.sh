# This demonstrates that the multi-shell feature works and a topology can be generated for a multi-shell constellation
# (here the constellation is a combination of Iridium and Kuiper, as they are close in orbital altitude)

echo "MULTI-SHELL EXAMPLE"
echo ""
echo "Please note that you may have to brew install proj to ensure pip install geopandas."
echo "Starting report recreation ..."

echo -n "Multishell Novel Algorithm Topology (Primal DCMST): "
python __main__.py --tles multishell-constellation_tles.txt.tmp --constellation Multishell --m 6 34 --n 11 34 --i 86.4 51.9 --rev 14.8 14.8 --multi True --optimise False --topology novel --dcmst "primal" --isl_terminals 4 --snapshot_interval 60 --weights 1 1 1

echo""
echo "EXAMPLE TOPOLOGIES CONSTRUCTED"