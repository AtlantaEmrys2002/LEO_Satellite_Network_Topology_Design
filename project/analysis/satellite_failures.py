import numpy as np
from measure import read_isl_file
import networkx as nx
import random


def satellite_failure(num_snapshots: int, topology_file_location: str, num_satellites: int, constellation_name: str,
                      prob_failure: float):

    if prob_failure > 1 or prob_failure < 0:
        raise ValueError("Probabilities must be in [0, 1].")

    av_disjoint_paths = 0
    av_disconnected = 0

    # For each topology over the course of one orbital period
    for k in range(num_snapshots):

        # Read in current topology
        topology = read_isl_file(topology_file_location + "/isls_" + str(k) + ".txt", num_satellites)

        # Read in sunlight matrix
        sunlight = np.load("./" + constellation_name + "/sunlight_matrices/sunlight_matrix_" + str(k) + ".npy")

        # Multiply sunlight matrix by probability of satellite failure
        sunlight *= prob_failure

        # Calculate whether each satellite is in sunlight or not (i.e. columns of all ones
        sunlit_satellites = np.array([k for k in range(num_satellites) if np.all(sunlight[k] == 1)])

        failed_satellites = []

        # Randomly select satellites to fail all at once due to solar flares - MONTE CARLO
        for i in range(len(sunlit_satellites)):
            rand_val = random.random()
            if rand_val < prob_failure:
                failed_satellites.append(sunlit_satellites[i])

        # Set all columns and rows to 0 where satellite failed to 0 (no connections)
        for v in failed_satellites:
            topology[:, v] = 0
            topology[v, :] = 0

        # Calculate the number of disjoint paths between all satellite pairs that contain none of the satellites that
        # will fail
        G = nx.from_numpy_array(topology)

        # Remove failed satellites
        for f in failed_satellites:
            G.remove_node(f)

        # Check if subgraph is connected
        connected = nx.is_connected(G)

        # Counts number of disjoint paths between all remaining satellite pairs
        path_counter = 0

        for i in range(num_satellites):
            for j in range(i+1, num_satellites):
                path_counter += len(list(nx.node_disjoint_paths(G, i, j)))

        path_counter /= ((pow(num_satellites, 2) - num_satellites)/2)

        if connected is False:
            av_disconnected += 1

        av_disjoint_paths += path_counter

    return av_disjoint_paths/num_snapshots, av_disconnected
