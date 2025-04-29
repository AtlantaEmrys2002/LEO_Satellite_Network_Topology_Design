import argparse
import numpy as np
from measure import read_isl_file
import networkx as nx
import random


def satellite_failure(num_snapshots: int, topology_file_location: str, num_satellites: int, constellation_name: str,
                      prob_failure: float, topology_type: str = "dynamic"):

    if prob_failure > 1 or prob_failure < 0:
        raise ValueError("Probabilities must be in [0, 1].")

    av_disjoint_paths = 0
    av_disconnected = 0

    total_failed_satellites = 0

    # For each topology over the course of one orbital period
    for k in range(num_snapshots):

        if topology_type == "static":
            topology = read_isl_file(topology_file_location + "/isls_0.txt", num_satellites)
        else:
            # Read in current topology
            topology = read_isl_file(topology_file_location + "/isls_" + str(k) + ".txt", num_satellites)

        # Read in sunlight matrix
        sunlight = np.load("./" + constellation_name + "/sunlight_matrices/sunlight_matrix_" + str(k) + ".npy")

        # # Multiply sunlight matrix by probability of satellite failure
        # sunlight *= prob_failure

        # Calculate whether each satellite is in sunlight or not (i.e. columns of all ones)
        # sunlit_satellites = np.array([k for k in range(num_satellites) if np.all(sunlight[k] == 1)])
        sunlit_satellites = np.array([m for m in range(num_satellites) if np.all(np.isclose(sunlight[m], 1.0))])

        failed_satellites = []

        # Randomly select satellites to fail all at once due to solar flares - MONTE CARLO
        for i in range(len(sunlit_satellites)):
            rand_val = random.random()
            if rand_val < prob_failure:
                failed_satellites.append(sunlit_satellites[i])

        total_failed_satellites += len(failed_satellites)

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

        # for i in range(num_satellites):
        #     for j in range(i+1, num_satellites):
        #         try:
        #             path_counter += len(list(nx.node_disjoint_paths(G, i, j)))
        #         except nx.NetworkXNoPath:
        #             print("No path between satellites" + str(i) + " and " + str(j) + ".")

        # path_counter /= ((pow(num_satellites, 2) - num_satellites)/2)

        if connected is False:
            av_disconnected += 1

        av_disjoint_paths += path_counter

    print('\n')
    print("Total Satellites Failures: " + str(total_failed_satellites))
    print("Average Disjoint Paths: " + str(av_disjoint_paths/num_snapshots))
    print("Number of Disconnected Component Occurrences: " + str(av_disconnected))
    print('\n')

    # return av_disjoint_paths/num_snapshots, av_disconnected


if __name__ == "__main__":

    print("Calculating failures...")

    # Parse inputs to module
    parser = argparse.ArgumentParser()

    parser.add_argument("--num_snapshots", type=int, help="number of snapshots for which a topology was "
                                                          "constructed", required=True)

    parser.add_argument("--topology_location", type=str, help="location of files containing topology",
                        required=True)

    parser.add_argument("--num_sat", type=int, help="number of satellites within the constellation",
                        required=True)

    parser.add_argument("--constellation_name", type=str, help="name of constellation for which topology"
                                                               "was constructed", required=True)

    parser.add_argument("--prob_failure", type=float, help="probability of failure due to solar flare",
                        required=True)

    parser.add_argument("--topology_type", type=str, help="type of topology - 'dynamic' or 'static'. "
                                                          "Default is dynamic.")

    args = parser.parse_args()

    num_snapshots = args.num_snapshots
    topology_location = args.topology_location
    num_sat = args.num_sat
    constellation_name = args.constellation_name
    prob_failure = args.prob_failure
    topology_type = "dynamic"

    if args.topology_type:
        topology_type = args.topology_type

    satellite_failure(num_snapshots, topology_location, num_sat, constellation_name, prob_failure,
                      topology_type=topology_type)
