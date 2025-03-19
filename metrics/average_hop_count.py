import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import floyd_warshall

def hop_count(topology_matrix, distance_matrix, num_satellites):

    # Calculate distance matrix for dijkstra function, such that only distances between satellites with an active ISL
    # are included
    distance_matrix = np.where(topology_matrix == 1, distance_matrix, 0)

    # Calculate path from source to destination using Dijkstra's Shortest Path Algorithm
    graph = csr_array(distance_matrix)
    _, prev = floyd_warshall(csgraph=graph, directed=False, return_predecessors=True)

    print("NEEDS IMPLEMENTING")


    # return average_hop_count