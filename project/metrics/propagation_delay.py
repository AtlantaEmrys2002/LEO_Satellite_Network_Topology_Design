# Libraries
from astropy.constants import c
import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import dijkstra


def propagation_delay(topology_matrix: np.ndarray, distance_matrix: np.ndarray, num_satellites: int) \
        -> tuple[float, float]:
    """
    Calculates the max propagation delay and mean propagation delay for a satellite network for all paths between all
    satellite pairs within the satellite network.
    :param topology_matrix: an adjacency matrix representing a satellite network topology
    :param distance_matrix: an adjacency matrix where each element represents the distance between a satellite pair
    :param num_satellites: the number of satellites in the network
    :return:
    """
    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included
    distance_matrix = np.where(topology_matrix == 1, distance_matrix, 0)

    # Calculate path from source to destination using Dijkstra's Shortest Path Algorithm
    graph = csr_array(distance_matrix)
    dist = np.array(dijkstra(csgraph=graph, directed=False))

    # Check if any path between nodes where satellite a cannot reach satellite b
    if len(dist[dist == -9999.]) != 0:
        raise ValueError("Cannot calculate propagation delay - a path does not exist between 2 satellites in the "
                         "network.")

    # Undirected graph - only consider unique paths between satellites (i.e. do not consider path from 2 to 1 if 1 to 2
    # has already been considered)
    for i in range(num_satellites):
        for j in range(i):
            dist[i, j] = 0

    # Calculate latencies using velocity = distance/time formula
    dist /= c.to('km/s').value

    # Calculate propagation delay
    max_pd = np.max(dist)

    # Calculate mean propagation delay - don't consider paths between satellites and themselves. Additionally, this is
    # an undirected graph
    mean_pd = np.sum(dist) / ((num_satellites ** 2 - num_satellites)/2)

    return max_pd, mean_pd
