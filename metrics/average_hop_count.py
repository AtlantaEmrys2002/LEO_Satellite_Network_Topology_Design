# Libraries
import numpy as np
import networkx as nx

def hop_count(topology_matrix, distance_matrix, num_satellites):
    """
    Calculate the average hop count for a network topology based on calculations of hop count between all satellite pairs in the network.
    :param topology_matrix:
    :param distance_matrix:
    :param num_satellites:
    :return:
    """
    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included
    distance_matrix = np.where(topology_matrix == 1, distance_matrix, 0)

    # Calculate path from source to destination using Dijkstra's Shortest Path Algorithm
    graph = nx.from_numpy_array(distance_matrix)

    print("HI")

    # result = nx.shortest_path(graph)
    path = dict(nx.all_pairs_shortest_path(graph))

    average_hop_count = 0

    # Calculate hop counts between all pairs of satellites
    for i in range(num_satellites):
        for j in range(i+1, num_satellites):
            # average_hop_count += (len(result[i][j]) - 1)
            average_hop_count += (len(path[i][j]) - 1)

    # Calculate average
    average_hop_count /= ((pow(num_satellites, 2) - num_satellites) / 2)

    return average_hop_count