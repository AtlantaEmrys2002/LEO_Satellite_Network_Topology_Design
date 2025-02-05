#TODO
# Based on theoretical recommendations by Network Performance of pLEO Topologies in a
# High-Inclination Walker Delta Satellite Constellation
# Involves using Monte Carlo methods to evaluate different topology designs
# Satellite links already in adjacency matrix (need to create function in heuristic_topology_
# design_algorithm_isls.py that allows you to specify format (file or adjacency matrix)
# Need to figure out latency matrix. Use Shortest Path (Dijkstra) to find shortest paths between all pairs of
# vertices (need to state that that is our routing algorithm)
# Can add in random failures (like original paper)
# CALCULATE propagation latency using classic speed = distance/time with speed being c (speed of light) and
# distance being straight-line distance between satellites
# NEED 2 FUNCtiONS - one with random failures and one without random failures (see which
# set of params works best with random failures - works well with solar flares idea)
# Say assume no queueing delay at actual satellites
# NEED TO DO FOR DIFFERENT SNAPSHOTS OF NETWORK
# Put everything in main function (make sure all tidy and neat)
# State in report precision of signal_speed (speed of light used)
# NEED TO CHECK IF DISTANCE MATRIX IN KM OR M - SO CAN MATCH SPEED OF LIGHT UNITS
# See if need to calculate number of nodes
# Adapt Dijkstra to use priority queue?
# might need case where - raise error if no path between two satellites in Dijkstra
# NEED LOTS OF ERROR HANDLING IN BOTH FILES
# Do i need to calculate the number of satellites or jsut pass to function - cleaner to calculate length once (so no confusion between
# function calls)
# RANDOMLY SELECT SNAPSHOTS SO NOT ALL FOR SAME NETWORK OR PERFORM FOR MULTIPLE RANDOM SNAPSHOTS (EACH RUN, NETWORK IS STATIC)
# If update alg to use cone of visibility from Load Balancing Paper and Energy Efficiency - three equations to calculate energy
# consumed by network - include function to calculate total energy consumed by topology (possibly just maintenance energy,
# as already considering on/off energy as looking at link churn). Energy efficiency also tied to distance.

# CHANGE SO MATRIX BUILT FROM ISL FILE


# Import Libraries
from astropy.constants import c
from random import sample
import numpy as np
import heuristic_topology_design_algorithm_isls as topology_build # NEED TO CREATE TEST FUNCTION THAT RETURNS MATRICES RATHER THAN FILE

### Link Churn Calculation ###
# Takes all snapshots of topology and finds the total number of link changes over the course of 1 orbit - need all snapshots!
def link_churn(snapshots):
    # Convert to numpy arrays
    snapshots = [np.array(x) for x in snapshots]

    # Initialise total link changes to 0
    total = 0

    # CHECK IT WORKS WITH 0 and 1 (INSTEAD OF T AND F)
    for x in range(len(snapshots)):
        if x != len(snapshots) - 1:
            total += np.sum(np.logical_xor(snapshots[x], snapshots[x + 1]))

    return total

### Latency and Hop Count Calculations ###

# Function to calculate propagation delay (in seconds) between satellite pair
# Set default signal speed to speed of light. Uses time = propagation speed * total distance of ISL path between
# satellites
def propagation_latency(total_distance):
    return total_distance * c.value

# Implementation of Dijkstra's Algorithm - used to calculate the length of the shortest path between two satellites
def dijkstras_algorithm(distance_matrix, source, destination, num_nodes):
    # Stores unexplored vertices
    q = list(range(num_nodes))
    # Stores distance to vertex v from source
    dist = [-1 for _ in range(num_nodes)]
    # Stores previous node in the shortest path from source u to vertex v
    prev = [-1 for _ in range(num_nodes)]

    # Initialise source
    dist[source] = 0
    q.remove(source)
    current_node = source

    # Find Shortest Path
    while len(q) != 0:
        # NEED TO CHECK ALL DISTANCES MUST BE GREATER THAN 0 IF ISL EXISTS BETWEEN THEM
        current_node = min([distance_matrix[current_node].index(x) for x in distance_matrix[current_node] if x not in q and x > 0])

        if current_node == destination:

            # Calculate hop count of path between two satellites - count must be at least one (as no satellite ever
            # establishes ISL with itself)
            hop_count = 1
            new_node = prev[current_node]
            while new_node != source:
                new_node = prev[new_node]
                hop_count += 1
            return [dist[current_node], hop_count]

        q.remove(current_node)

        for v in range(num_nodes):
            if v not in q:
                continue
            else:
                alt = dist[current_node] + distance_matrix[current_node][v]
                if alt < dist[v] or dist[v] == -1:
                    dist[v] = alt
                    prev[v] = current_node

    # If destination node has been found
    if (dist[current_node] != -1) and (current_node == destination):
        hop_count = 1
        new_node = prev[current_node]
        while new_node != source:
            new_node = prev[new_node]
            hop_count += 1
        return [dist[current_node], hop_count]
    else:
        raise ValueError("No path exists between given satellite pair.")


# Calculates the max propagation delay, mean propagation delay, and average hop count for a given satellite network
# topology
# num_random_pairs states the number of random satellite pairs to calculate measurements for (default is 1000 based on
# pLEO paper) - could change
def latency_hop_count_calculation(topology_matrix, distance_matrix, num_random_pairs=1000):

    # Calculate the number of satellites in the network
    num_satellites = len(distance_matrix[0])

    # Generate random satellite pairs
    satellite_pairs = [sample(range(num_satellites),2) for _ in range(num_random_pairs)]

    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included - USE NUMPY MASK WHEN CHANGE TO NUMPY
    for i in range(num_satellites):
        for j in range(num_satellites):
            if topology_matrix[i][j] == 0:
                distance_matrix[i][j] = 0

    # Calculate distance and hop count between each pair of vertices
    path_distances = [dijkstras_algorithm(distance_matrix, satellite_pairs[i][0], satellite_pairs[i][1], num_satellites - 1) for i in range(num_random_pairs)]

    # Calculate propagation latencies between each pair of satellites
    prop_latencies = [propagation_latency(x[0]) for x in path_distances]

    # Calculate mean and max propagation latency for topology
    max_latency = max(prop_latencies)
    mean_latency = sum(prop_latencies) / num_random_pairs

    # Calculate average hop count for topology
    average_hop_count = sum([x[1] for x in path_distances])

    return max_latency, mean_latency, average_hop_count


# Run with different weights for cost function
def experiment():
    alpha_values = [x/10 for x in range(11)]
    beta_values = [x / 10 for x in range(11)]
    gamma_values = [x / 10 for x in range(11)]

    # Store Results
    results = []

    for x in alpha_values:
        for y in beta_values:
            for z in gamma_values:
                # RUN TEST FUNCTION THAT RETURNS ALL THE MATRICES AND DOESN'T SAVE TO FILE - MAY WANT TO SAVE COPY OF TOPOLOGY THOUGH IF MAKE PROBABILISTIC




# Main Function - arguments are topology matrix, distance matrix - calculates set of random pairs (e.g. a thousand random
# satellite pairs) - stores then calls Dijkstra to find the shortest path between each set of pairs (storing the sequence of
# satellites visited along path, e.g. 70, 31, 54). Find distances of paths and stores in corresponding distance matrix.
# Then, calculates propagation delay for paths (using distance matrix and seq from total distances.


# References:
# Dijkstra's Algorithm Pseudocode - https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
# If/Else List Comprehension Order - https://stackoverflow.com/questions/4406389/if-else-in-a-list-comprehension
# Python Documentation - https://docs.python.org/3/reference/simple_stmts.html#raise
# Astropy Documentation - https://docs.astropy.org/en/latest/index_user_docs.html
# Numpy Documentation - https://numpy.org/doc/2.1/user/index.html
# Raising Errors and Termination - https://stackoverflow.com/questions/65152176/does-raise-alone-make-the-program-terminate-or-not
# Monte Carlo Method - https://en.wikipedia.org/wiki/Monte_Carlo_method#
# IBM What is Monte Carlo Simulation? - https://www.ibm.com/think/topics/monte-carlo-simulation
# Numpy-List Conversion - https://stackoverflow.com/questions/10346336/list-of-lists-into-numpy-array