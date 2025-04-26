# Run visualisation code
echo "TOPOLOGY VISUALISATION"

echo "Visualising +Grid Topologies: "

python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/plus_grid/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method plus_grid --topology-type static --num_snapshots 97

python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/plus_grid/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method plus_grid --topology-type static --num_snapshots 90

python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/plus_grid/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method plus_grid --topology-type static --num_snapshots 105

echo "Visualising MDTD Topologies:"

python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/mdtd/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method mdtd --topology-type dynamic --num_snapshots 97

python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/mdtd/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method mdtd --topology-type dynamic --num_snapshots 90

python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/mdtd/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method mdtd --topology-type dynamic --num_snapshots 105

echo "Visualising Novel Topologies (Primal):"

python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/novel/random/primal/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 97
python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/primal/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 97

python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/novel/random/primal/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 90
python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/primal/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 90

python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/novel/random/primal/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 105
python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/primal/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 105

echo "Visualising Novel Topologies (Genetic):"

python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/novel/random/ga/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 97
python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/ga/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 97

python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/novel/random/ga/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 90
python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/ga/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 90

python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/novel/random/ga/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 105
python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/ga/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 105

echo "Visualising Novel Topologies (Ant-Based):"

python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/novel/random/aco/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 96
python analysis/visualisation.py --tles kuiper-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/aco/kuiper-630 --constellation_name Kuiper-630 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 96

python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/novel/random/aco/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 89
python analysis/visualisation.py --tles starlink-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/aco/starlink-550 --constellation_name Starlink-550 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 89

python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/novel/random/aco/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 104
python analysis/visualisation.py --tles telesat-constellation_tles.txt.tmp --location ./Results/novel/evolutionary/aco/telesat-1015 --constellation_name Telesat-1015 --snapshot_interval 60 --method novel --topology-type dynamic --num_snapshots 104

echo "Visualisation build completed."