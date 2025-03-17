#TODO
# Can add in random failures (like original paper)
# CALCULATE propagation latency using classic speed = distance/time with speed being c (speed of light) and
# distance being straight-line distance between satellites
# NEED 2 FUNCtiONS - one with random failures and one without random failures (see which
# set of params works best with random failures - works well with solar flares idea)
# State in report precision of signal_speed (speed of light used)
# If update alg to use cone of visibility from Load Balancing Paper and Energy Efficiency - three equations to calculate energy
# consumed by network - include function to calculate total energy consumed by topology (possibly just maintenance energy,
# as already considering on/off energy as looking at link churn). Energy efficiency also tied to distance.
# CHECK CSR ARRAY GENERATED CORRECTLY - 0 SHOULD BE NO DIST NOT 0KM

# Import Libraries
from astropy.constants import c
from scipy.sparse import csr_array
from scipy.sparse.csgraph import dijkstra
from scipy.sparse.csgraph import floyd_warshall
import csv
import random
from random import sample
import numpy as np
import time

# Seed random generators to ensure reproducible results
np.random.seed(42)
random.seed(42)

### Link Churn Calculation ###
# Takes all snapshots of topology and finds the total number of link changes over the course of 1 orbit - need all snapshots!
# def link_churn(snapshots):
#     # Convert to numpy arrays
#     snapshots = [np.array(x) for x in snapshots]
#
#     # Initialise total link changes to 0
#     total = 0
#
#     # CHECK IT WORKS WITH 0 and 1 (INSTEAD OF T AND F)
#     for x in range(len(snapshots)):
#         if x != len(snapshots) - 1:
#             total += np.sum(np.logical_xor(snapshots[x], snapshots[x + 1]))
#
#     return total

### Latency and Hop Count Calculations ###

# Function to calculate propagation delay (in seconds) between satellite pair
# Set default signal speed to speed of light. Uses time = propagation speed * total distance of ISL path between
# satellites
def propagation_latency(total_distance):
    return total_distance / c.to('km/s').value


def dijkstras_algorithm(distance_matrix, source, destination):

    # Calculate path from source to destination using Dijkstra's Shortest Path Algorithm
    graph = csr_array(distance_matrix)
    dist, prev, _ = dijkstra(csgraph=graph, directed=False, indices=source, return_predecessors=True, min_only=True)

    # If no path exists between nodes, throw ValueError
    if dist[destination] == -9999:
        raise ValueError("No path exists between given satellite pair.")

    # Source == destination
    if prev[destination] == -9999:
        raise ValueError("No path exists between given satellite pair.")

    # Calculate hop count from source to destination - hop count is always >= 1 as cannot calculate path from vertex to
    # itself
    hop_count = 1
    current_node = prev[destination]
    while current_node != source:
        if prev[current_node] == -9999:
            raise ValueError("No path exists between given satellite pair.")
        current_node = prev[current_node]
        hop_count += 1

    return dist[destination], hop_count


# Calculates the max propagation delay, mean propagation delay, and average hop count for a given satellite network
# topology
# num_random_pairs states the number of random satellite pairs to calculate measurements for (default is 1000 based on
# pLEO paper) - could change
def latency_hop_count_calculation(topology_matrix, distance_matrix, num_random_pairs=10):

    # Calculate the number of satellites in the network
    num_satellites = len(distance_matrix[0])

    # Generate random satellite pairs
    satellite_pairs = [sample(range(num_satellites),2) for _ in range(num_random_pairs)]

    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included
    distance_matrix = np.where(topology_matrix == 1, distance_matrix, 0)

    # Calculate distance and hop count between each pair of vertices
    path_distances = [dijkstras_algorithm(distance_matrix, satellite_pairs[i][0], satellite_pairs[i][1]) for i in range(num_random_pairs)]

    # Calculate propagation latencies between each pair of satellites
    prop_latencies = [propagation_latency(x[0]) for x in path_distances]

    # Calculate mean and max propagation latency for topology
    max_latency = max(prop_latencies)
    mean_latency = sum(prop_latencies) / num_random_pairs

    # Calculate average hop count for topology
    average_hop_count = sum([x[1] for x in path_distances]) / num_random_pairs

    return [max_latency, mean_latency, average_hop_count]


def latency_hop_count_floyd(topology_matrix, distance_matrix):
    # Calculate the number of satellites in the network
    num_satellites = len(distance_matrix[0])

    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included
    distance_matrix = np.where(topology_matrix == 1, distance_matrix, 0)

    # Calculate path from source to destination using Dijkstra's Shortest Path Algorithm
    graph = csr_array(distance_matrix)
    dist, prev = floyd_warshall(csgraph=graph, directed=False, return_predecessors=True)

    print(dist)


# Main Function - Monte Carlo Simulation. Arguments are topology matrix, distance matrix - calculates set of random
# pairs (e.g. a thousand random satellite pairs) - stores then calls Dijkstra to find the shortest path between each set
# of pairs (storing the sequence of satellites visited along path, e.g. 70, 31, 54). Find distances of paths and stores
# in corresponding distance matrix. Then, calculates propagation delay for paths (using distance matrix and seq from
# total distances. Experiments with different weights for cost function. Theoretical Monte Carlo simulations are similar
# to and inspired by those used in paper titled 'Network Performance of pLEO Topologies in a High-Inclination Walker
# Delta Satellite Constellation' (see report for full reference).
def experiment(num_param_sets, num_snapshots, num_snapshots_to_sample, num_satellites, constellation_name, num_random_pairs):

    # Temporarily Store Results
    results = []

    # Randomly sample sets of parameters (where alpha, beta, and gamma can be random variables)
    parameter_sets = np.random.rand(num_param_sets, 3)

    # Randomly sample snapshots to analyse for each set of parameters
    snapshot_ids = np.random.randint(0, num_snapshots, (num_param_sets, num_snapshots_to_sample))

    # Run experiments with given parameters
    for k in range(num_param_sets):

        # TEMPORARILY COMMENTED OUT
        # Build topologies with given parameters
        # topology_build.main("starlink-constellation_tles.txt.tmp", "Starlink-550", 72, 22, 53, 15.19, snapshot_ids[k], parameter_sets[k])

        # Fetch topologies generated by algorithm
        for j in snapshot_ids[k]:

            # Read in topology built for given snapshot
            topology_isls = np.loadtxt("./isl_topologies/" + constellation_name + "_isls_" + str(j) + ".txt").astype(int).T

            # Create topology matrix for ISLs
            topology_matrix = np.zeros((num_satellites, num_satellites))
            topology_matrix[topology_isls[0], topology_isls[1]] = 1

            max_latency, mean_latency, average_hop_count = latency_hop_count_calculation(topology_matrix, np.load("./"+constellation_name+"/distance_matrices/dist_matrix_"+str(j)+".npy"), num_random_pairs)

            results.append(dict(alpha=parameter_sets[k][0], beta=parameter_sets[k][1], gamma=parameter_sets[k][2], max_latency=max_latency, mean_latency=mean_latency, average_hop_count=average_hop_count))

            print('hi')

            start = time.time()

            latency_hop_count_floyd(topology_matrix, np.load("./"+constellation_name+"/distance_matrices/dist_matrix_"+str(j)+".npy"))

            print(time.time() - start)


            print('ji')

    # Write Results to CSV Format - this code was adapted from documentation - https://docs.python.org/3/library/csv.html#csv.DictWriter
    with open('results.csv', 'w', newline='') as csvfile:
        fieldnames = ['alpha', 'beta', 'gamma', 'max_latency', 'mean_latency', 'average_hop_count']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)


# At the moment, we are taking a snapshot every minute (therefore, num_snapshots is 94). Number of snapshots that are
# analysed per parameter set is currently 2
experiment(10, 94, 1, 1584, "Starlink-550", 1000)


# References:
# Astropy Documentation - https://docs.astropy.org/en/latest/index_user_docs.html
# CSV File - https://docs.python.org/3/library/csv.html#csv.DictWriter
# Dijkstra's Algorithm Pseudocode - https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
# Floyd-Warshall - https://medium.com/100-days-of-algorithms/day-65-floyd-warshall-2d10a6d6c49d
# If/Else List Comprehension Order - https://stackoverflow.com/questions/4406389/if-else-in-a-list-comprehension
# IBM What is Monte Carlo Simulation? - https://www.ibm.com/think/topics/monte-carlo-simulation
# Monte Carlo Method - https://en.wikipedia.org/wiki/Monte_Carlo_method#
# Multi-Objective Optimisation - https://en.wikipedia.org/wiki/Multi-objective_optimization
# Numpy Documentation - https://numpy.org/doc/2.1/user/index.html
# Numpy-List Conversion - https://stackoverflow.com/questions/10346336/list-of-lists-into-numpy-array
# Numpy Type Conversion - https://stackoverflow.com/questions/12648624/python-converting-an-numpy-array-data-type-from-int64-to-int
# # Python Documentation - https://docs.python.org/3/reference/simple_stmts.html#raise
# Raising Errors and Termination - https://stackoverflow.com/questions/65152176/does-raise-alone-make-the-program-terminate-or-not
# Scipy Documentation - https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.dijkstra.html#scipy.sparse.csgraph.dijkstra
# Writing Dictionaries to Files - https://pythonspot.com/save-a-dictionary-to-a-file/