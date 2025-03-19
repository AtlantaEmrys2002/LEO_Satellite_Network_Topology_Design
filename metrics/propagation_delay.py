# Libraries
from astropy.constants import c
import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import floyd_warshall

def propagation_delay(topology_matrix, distance_matrix, num_satellites):
    """
    Calculates the max propagation delay and mean propagation delay for a satellite network for all paths between all satellite pairs within the satellite network.
    :param topology_matrix:
    :param distance_matrix:
    :param num_satellites:
    :return:
    """
    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included
    distance_matrix = np.where(topology_matrix == 1, distance_matrix, 0)

    # Calculate path from source to destination using Dijkstra's Shortest Path Algorithm
    graph = csr_array(distance_matrix)
    dist = floyd_warshall(csgraph=graph, directed=False)

    # Check if any path between nodes where satellite a cannot reach satellite b
    if len(dist[dist == -9999.]) != 0:
        raise ValueError("Cannot calculate propagation delay - a path does not exist between 2 satellites in the network.")

    # Calculate latencies using velocity = distance/time formula
    dist /= c.to('km/s').value

    # Calculate propagation delay
    max_pd = np.max(dist)

    # Calculate mean propagation delay - don't consider paths between satellites and themselves. Additionally, this is an undirected graph
    mean_pd = np.sum(dist) / ((num_satellites ** 2 - num_satellites)/2)

    return max_pd, mean_pd



