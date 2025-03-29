# Libraries
import numpy as np
import networkx as nx

# Networkx Warning Suppression - warning already acknowledged and recommendation followed
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def hop_count(topology_matrix: np.ndarray, distance_matrix: np.ndarray, num_satellites: int) -> float:
    """
    Calculate the average hop count for a network topology based on calculations of hop count between all satellite
    pairs in the network.
    :param topology_matrix: an adjacency matrix representing a satellite network topology
    :param distance_matrix: an adjacency matrix where each element represents the distance between a satellite pair
    :param num_satellites: the number of satellites within the network
    :return: the average hop count for the given satellite network
    """
    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included
    distance_matrix = np.where(topology_matrix == 1, distance_matrix, 0)

    # Calculate path from source to destination using Dijkstra's Shortest Path Algorithm
    graph = nx.from_numpy_array(distance_matrix)
    result = dict(nx.shortest_path(graph))

    average_hop_count = 0

    # Calculate hop counts between all pairs of satellites
    for i in range(num_satellites):
        for j in range(i+1, num_satellites):
            average_hop_count += (len(result[i][j]) - 1)

    # Calculate average
    average_hop_count /= ((pow(num_satellites, 2) - num_satellites) / 2)

    return average_hop_count
