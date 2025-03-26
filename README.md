# Low Earth Orbit Satellite Network Topology Design

This Python package was created by Millicent Riordan for their third-year project *Designing High-Performance
Inter-Satellite Link Topologies for Satellite Megaconstellations* at Durham University. 

## Set-Up

### System

1. Python version 3.11+
2. To ensure full operation, recommended operating systems include recent versions of Linux (e.g. Ubuntu 18+, similar to the Hypatia simulation software) or macOS 14+ (Sonoma onwards).

### Python Dependencies

Dependencies for this package include the Python libraries: `astropy`, `certifi`, `ephem`, `networkx`, `numpy`, `nx-parallel`, `pytest`, `scipy`, `sgp4`, `skyfield`, and `sphinx`. 

To install all dependencies automatically, run ```bash install_dependencies.sh```if you wish to install on a macOS system or ```bash install_dependencies_linux.sh``` if you wish to install on a Linux-based system.

## Modules and Functionality

This package involves numerous functions that are grouped as follows:

NEEDS DESCRIPTIONS INCLUDED

- `analysis`
- `cost_function_optimisation_algorithms`
- `data_handling`
- `dcmst_construction_algorithms`
- `metrics`
- `satellite_network_attribute_functions`
- `satellite_topology_construction_algorithms`
- `tests`

## Tutorial and Documentation

Documentation of each function within this Python package can be found [here](). 

## Recreating Report

To recreate the results in the project's final report, please run ```bash report.sh```. Please note that this scripts 
takes significant execution time. However, a large proportion of this code is automatically parallelised -  if you have 
access to a significant number of CPU cores, this code will run significantly faster.
