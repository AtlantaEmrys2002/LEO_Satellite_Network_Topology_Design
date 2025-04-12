# Libraries
from collections import deque
import copy
import numpy as np
import random
from scipy.cluster.hierarchy import DisjointSet


def update_edge_pheromone(edges: np.ndarray, eta: float, minPhm: float, maxPhm: float) -> np.ndarray:
    """
    Used ACO DCMST construction algorithm to update pheromones of each potential edge in graph.
    :param edges: list of edges in graph with their respective values of edge cost, initial pheromone deposit, current
     pheromone deposit, and number of ants that have traversed the edge
    :param eta: hyperparameter used to determine amount to update edge pheromone
    :param minPhm: minimum pheromone that can be on any edge
    :param maxPhm: maximum pheromone that can be on any edge
    :return: returns list of edges where the pheromone level of each edge has been updated according to set of
     conditions
    """
    const = 1 - eta

    # phm = (1 - eta) * phm + nVisited * initPhm
    edges[:, 4] = const * edges[:, 4] + edges[:, 5] * edges[:, 3]

    # set nVisited to 0
    edges[:, 5] = 0

    for e in edges:
        # If pheromones are greater than maximum amount possible on edge, or less than minimum amount possible on edge
        if e[4] > maxPhm:
            e[4] = maxPhm - e[3]
        if e[4] < minPhm:
            e[4] = minPhm - e[3]

    return edges


def min_max(weights: np.ndarray) -> np.ndarray:
    """
    Min-max scales array of values to range [0, 1].
    :param weights: array of values to be min-max scaled (in this case, pheromone levels for each edge in graph)
    :return: min-max scaled values of weights parameter.
    """
    return (weights - np.min(weights)) / (np.max(weights) - np.min(weights))


# def initialise_ants_and_edges(cost_matrix: np.ndarray, num_sat: int) -> tuple[list[Ant], np.ndarray, float, float]:
def initialise_ants_and_edges(cost_matrix: np.ndarray, num_sat: int) -> tuple[list, np.ndarray, float, float]:
    """
    Initialises ants and edges (they are associated with other values, such as pheromones) for ACO DCMST construction
    algorithm.
    :param cost_matrix: costs assigned to each edge within the graph that represents the satellite network
    :param num_sat: the number of satellites within the network
    :return: the initialised ants and edges (in the correct format), as well as the maximum pheromone level and minimum
     pheromone level that may be assigned to any given edge
    """
    # CREATE ANTS #

    # Initialise an ant at each vertex (cardinality will always be the number of satellites in the network) with a tabu
    # list of all vertices ant has recently visited
    ants = [[v, deque([])] for v in range(num_sat)]

    # CREATE EDGES #

    max_cost = np.max(cost_matrix)
    min_cost = np.min(cost_matrix[cost_matrix > 0])

    # Select all possible edges (remove duplicates, as undirected)
    graph_edges = np.argwhere(cost_matrix > 0)
    graph_edges = np.unique(np.sort(graph_edges), axis=0)

    # Each edge in array consists of u, v (the two nodes incident to edge), the edge cost, the initial pheromone deposit
    # on the edge, the current pheromone deposit on the edge (same as initial when initialised), and the number of ants
    # that have traversed the edge
    edges = np.array([[u[0], u[1], cost_matrix[u[0], u[1]], (max_cost - cost_matrix[u[0], u[1]]) + (max_cost - min_cost)
                       / 3, (max_cost - cost_matrix[u[0], u[1]]) + (max_cost - min_cost) / 3, 0] for u in graph_edges])

    # Calculate min and max pheromone each edge can have
    maxPhm = 1000 * (max_cost - min_cost) + (max_cost - min_cost) / 3
    minPhm = (max_cost - min_cost) / 3

    return ants, edges, maxPhm, minPhm


def modified_kruskal(edges: np.ndarray, constraints: np.ndarray, nCandidates: int, num_sat: int) -> (
        tuple)[list[np.ndarray], dict]:
    """
    Modified version of Kruskal's algorithm used to construct a DCMST - used in the ACO algorithm to create a DCMST
    based on pheromones laid by ants and edge costs.
    :param edges: edges within the graph
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param nCandidates: the number of candidate edges to evaluate at a time
    :param num_sat: the number of satellites within the network
    :return: a degree-constrained minimum spanning tree of the graph and a dict of the indices of incident edges for
    each vertex
    """

    # Initialise spanning tree and degrees of each vertex
    T_n = []
    degrees = np.array([0 for _ in range(num_sat)])

    # Sort edges in decreasing order according to pheromone level
    edges = edges[np.argsort(edges[:, 4])]

    # Select top nCandidates edges from E
    if len(edges) > nCandidates:
        candidate_edges = copy.deepcopy(edges[:nCandidates])
    else:
        candidate_edges = copy.deepcopy(edges)

    # Sort edges according to increasing edge cost
    candidate_edges = candidate_edges[np.argsort(candidate_edges[:, 2])]

    # Initialise Disjoint set data structure of all nodes
    disjoint_set = DisjointSet(list(range(num_sat)))

    level_of_candidates = 1

    num_edges_tree = num_sat - 1

    while len(T_n) < num_edges_tree:

        candidate = candidate_edges[0]

        candidate_edges = candidate_edges[1:]

        u = int(candidate[0])
        v = int(candidate[1])

        # Check if degree constraint violated by adding edge to tree
        if degrees[u] < constraints[u] and degrees[v] < constraints[v]:

            if disjoint_set.connected(u, v) is False:
                disjoint_set.merge(u, v)
                T_n.append(candidate)
                degrees[u] += 1
                degrees[v] += 1

        if len(candidate_edges) == 0:
            # Select next best candidates (if more edges are needed)
            if len(edges) <= nCandidates:
                raise ValueError("DCMST cannot be constructed, as not enough edges.")
            else:
                # Calculate new slice indices
                lowest = level_of_candidates * nCandidates
                level_of_candidates += 1
                highest = level_of_candidates * nCandidates
                if highest > len(edges) and lowest > len(edges):
                    raise ValueError("DCMST cannot be constructed, as not enough edges.")
                elif highest > len(edges):
                    candidate_edges = copy.deepcopy(edges[lowest:])
                else:
                    candidate_edges = copy.deepcopy(edges[lowest:highest])

                # Sort C in order of increasing edge cost
                candidate_edges = candidate_edges[np.argsort(candidate_edges[:, 2])]

    in_edges = incident_edges(edges, num_sat)

    return T_n, in_edges


def move_ant(a: list, edges: np.ndarray, in_edges: dict) -> list:
    """
    Moves ant through graph, effectively performing the "exploration phase" of the ant. Moves according to set number of
    rules - see comments and original paper for more details.
    :param a: list of ants
    :param edges: edges within the graph
    :param in_edges: dictionary of vertices and the indices of the edges incident to them
    :return: edges within the graph
    """
    nAttempts = 0
    moved = False

    # Move between five vertices or until ant cannot move (which over occurs first).
    while moved is False and nAttempts < 5:

        v_1 = a[0]

        # Select edges incident to v_1
        potential_edge_indices = in_edges[str(int(v_1))]

        # print(np.array_equal(potential_edge_indices, potential_edge_indices_2))

        potential_edges = edges[potential_edge_indices]

        # Select an edge adjacent to the vertex at which ant is located randomly (but proportional to the
        # pheromone on that edge). If no more edges to explore, increase nAttempts

        if len(potential_edges) != 0:

            random_edge = random.choices(list(range(len(potential_edges))), min_max(potential_edges[:, 4]).
                                         tolist())[0]

            # Find neighbouring vertex
            if potential_edges[random_edge, 0] == v_1:
                v_2 = potential_edges[random_edge, 1]
            else:
                v_2 = potential_edges[random_edge, 0]

            # If ant has not visited that vertex recently, move ant to that vertex
            if v_2 not in a[1]:
                a[1].append(v_2)
                a[0] = v_2

                edges[potential_edge_indices[random_edge], 5] += 1
                moved = True
            else:
                nAttempts += 1
        else:
            nAttempts += 1
            break

    return a


def incident_edges(edges: np.ndarray, num_sat: int) -> dict:
    """
    Creates a dictionary, such that the key is the satellite/vertex and the values are the indices (within the edges
    array of all edges incident to the vertex).
    :param edges: edges within the graph
    :param num_sat: the number of satellites within the network
    :return: dictionary of vertices and the indices of their incident edges
    """

    # Initialise dictionary
    incident = dict()

    view = edges.T

    for v in range(num_sat):

        # Find all edges where one vertex is equal to currently analysed vertex
        first = view[0] == v
        second = view[1] == v

        # Store indices of relevant edges
        potential_edge_indices = np.concatenate((np.flatnonzero(first), np.flatnonzero(second)))

        # Add vertex and indices as key-value pair to dictionary
        incident.update({str(v): potential_edge_indices})

    return incident


def move_ants(ants: list, edges: np.ndarray, max_steps: int, update_period: float, eta: float, minPhm: float,
              maxPhm: float, in_edges: dict) -> tuple[list, np.ndarray]:
    """
    Moves edges according to given constraints within the graph. Each ant's position is updated accordingly, as is the
    pheromone level of each edge. This is often deemed the 'exploration phase' of the algorithm
    :param ants: the ants used to explore the graph
    :param edges: edges within the graph
    :param max_steps: maximum number of edges each ant can traverse
    :param update_period: determines when to update edge pheromones (don't update every time an ant moves)
    :param eta: hyperparameter used to determine amount to update edge pheromone
    :param minPhm: minimum pheromone that can be on any edge
    :param maxPhm: maximum pheromone that can be on any edge
    :param in_edges: dictionary of vertices and the indices of the edges incident to them
    :return: returns updated ants (according to position and nodes recently visited) and edges (according to pheromone
     levels)
    """

    # The ants explore for a given number of steps
    for s in range(max_steps):

        # Update pheromones of edges after a given number of steps - helps reduce execution time.
        if s % update_period == 0:
            edges = update_edge_pheromone(edges, eta, minPhm, maxPhm)

        # Move each ant
        ants = [move_ant(a, edges, in_edges) for a in ants]

    return ants, edges


def solution_fitness(tree) -> float:
    """
    Calculates the sum of all edge costs within a given DCMST - used by the ACO DCMST construction algorithm.
    :param tree: a proposed DCMST of the given graph
    :return: the fitness of a given proposed DCMST.
    """
    return np.sum(np.asarray(tree)[:, 2])


def ant_colony(cost_matrix, constraints, num_sat: int, max_iterations: int = 100,
               max_iterations_without_improvement: int = 25, max_steps: int = 12, eta: float = 0.5,
               gamma: float = 1.5, eta_change: float = 0.95, gamma_change: float = 1.05, R: int = 30) \
        -> tuple[np.ndarray, np.ndarray]:
    """
    Algorithm that constructs a DCMST using an ant-based algorithm. Adapted from T. N. Bui, X. Deng, and C. M. Zrncic,
    “An Improved Ant-Based Algorithm for the Degree-Constrained Minimum Spanning Tree Problem,” IEEE Trans. On Evol.
    Computation, vol. 16, no. 2, pp. 266-278, Apr. 2012., doi: 10.1109/TEVC.2011.2125971. Default values for function
    recommended by this original paper.
    :param cost_matrix: an adjacency matrix, such that element cost_matrix[i][j] represents the cost of the graph edge
     ij
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param num_sat: the number of satellites within the network
    :param max_iterations: the maximum number of iterations performed by the function. Changed paper default from 10000.
    :param max_iterations_without_improvement: the maximum number of iterations performed by the function without any
     improvement before the current best topology is returned
    :param max_steps: the maximum number of exploration steps performed by each ant per iteration. Changed paper default
     from 75.
    :param eta: hyperparameter used to determine amount to update edge pheromone
    :param gamma: hyperparameter used to "enhance" the edge pheromones of the edges in the best spanning tree found so
     far
    :param eta_change: hyperparameter that determines the amount to adjust eta by each iteration
    :param gamma_change: hyperparameter that determines the amount to adjust gamma by each iteration
    :param R: used to prevent algorithm getting stuck in local optima. Changed from 100 to 30
    :return: a DCMST and the degree of each vertex within the tree
    """

    # Initialise counters
    i = 1
    i_best = 0
    i_restart = 0

    # Initialise candidate set cardinality (value recommended by original paper)
    candidate_set_cardinality = 5 * num_sat

    # Initialise update period (value recommended by original paper)
    update_period = max_steps / 3

    # Initialise |V| ants and edges with corresponding pheromone levels
    ants, edges, maxPhm, minPhm = initialise_ants_and_edges(cost_matrix, num_sat)

    # Construct spanning tree
    best_spanning_tree, in_edges = modified_kruskal(edges, constraints, candidate_set_cardinality, num_sat)

    # Calculate fitness of spanning tree
    best_fitness = solution_fitness(best_spanning_tree)

    # Continues iteratively improving DCMST until either the maximum number of iterations exceeded or no improvements
    # have been made in a set number of iterations

    while i < max_iterations and (i - i_best) < max_iterations_without_improvement:

        # ANT EXPLORATION #
        ants, edges = move_ants(ants, edges, max_steps, update_period, eta, minPhm, maxPhm, in_edges)

        # Construct new spanning tree based on exploration
        T, in_edges = modified_kruskal(edges, constraints, candidate_set_cardinality, num_sat)

        # Compare the fitness of current best solution to new one
        current_solution_fitness = solution_fitness(T)

        # Update spanning tree if better DCMST has been found
        if current_solution_fitness < best_fitness:
            best_spanning_tree = T
            best_fitness = current_solution_fitness
            i_best = i

        # ENHANCE #

        # Enhance edges (lay pheromones) - update pheromones of all edges in the best solution found so far

        # Find indices of all edges in best spanning tree and update their pheromone level accordingly

        for edge in best_spanning_tree:
            # Find edge in main set of edges and update accordingly
            index_of_edge = np.where((edges[:, 0] == edge[0]) & (edges[:, 1] == edge[1]))[0]
            edges[index_of_edge, 4] *= gamma

        # RESTART #

        # This prevents algorithm getting stuck in local optimum - restart algorithm (i.e. adjust edge pheromones)
        if i - max(i_best, i_restart) > R:
            i_restart = i
            for edge in best_spanning_tree:
                # Find edge index
                index_of_edge = np.where((edges[:, 0] == edge[0]) & (edges[:, 1] == edge[1]))[0]
                # Restart
                edges[index_of_edge, 4] *= random.uniform(0.1, 0.3)

        i += 1

        # RESET ANTS #

        # In next iteration, reset ants - approximately half (coin toss) of all ants remain in their current locations
        for a in ants:
            coin_flip = random.random()
            if coin_flip < 0.5:
                a[0] = random.randint(0, num_sat - 1)
            a[1].clear()

        # Update gamma and eta
        gamma *= gamma_change
        eta *= eta_change

    best_spanning_tree_adjacency = np.zeros((num_sat, num_sat))

    for edge in best_spanning_tree:

        best_spanning_tree_adjacency[int(edge[0]), int(edge[1])] = 1
        best_spanning_tree_adjacency[int(edge[1]), int(edge[0])] = 1

    return best_spanning_tree_adjacency, np.sum(best_spanning_tree_adjacency, axis=1).astype(np.int32)

# References
# Accelerating Masking - https://stackoverflow.com/questions/64109948/how-to-accelerate-numpy-array-masking
# Array of Objects - https://stackoverflow.com/questions/4877624/numpy-array-of-objects
# Column to Zero - https://stackoverflow.com/questions/37251992/set-a-column-in-numpy-array-to-zero
# Converting Tuples - https://stackoverflow.com/questions/47389447/how-convert-a-list-of-tuples-to-a-numpy-array-of-
# tuples
# Cycle Detection - https://en.wikipedia.org/wiki/Cycle_(graph_theory)
# Deep Copies - https://stackoverflow.com/questions/37593013/deep-copy-of-a-np-array-of-np-array
# Dictionaries - https://stackoverflow.com/questions/8424942/creating-a-new-dictionary-in-python
# Disjoint-set Data Structure - https://en.wikipedia.org/wiki/Disjoint-set_data_structure
# Extracting Rows - https://stackoverflow.com/questions/61025485/numpy-array-how-to-extract-whole-rows-based-on-values-
# in-a-column
# Faster Indexing - https://stackoverflow.com/questions/74482961/faster-numpy-array-indexing-when-using-condition-numpy-
# where
# Incident Edges - https://stackoverflow.com/questions/72134948/how-to-find-the-list-of-edges-incident-to-a-particular-
# vertex
# Incident Edges Networkx - https://stackoverflow.com/questions/33078907/get-all-edges-linked-to-a-given-node-in-a-
# networkx-graph
# Kruskal Algorithms - https://en.wikipedia.org/wiki/Kruskal%27s_algorithm#:~:text=Kruskal's%20algorithm%20finds%20a%20
# minimum,will%20not%20form%20a%20cycle.
# List Operations - https://stackoverflow.com/questions/2703310/list-comprehension-map-and-numpy-vectorize-performance
# Map Efficiency - https://stackoverflow.com/questions/35215161/most-efficient-way-to-map-function-over-numpy-array
# Matching Rows - https://stackoverflow.com/questions/25823608/find-matching-rows-in-2-dimensional-numpy-array
# Min-Max Normalisation - https://stackoverflow.com/questions/48178884/min-max-normalisation-of-a-numpy-array
# Numpy Equality - https://stackoverflow.com/questions/56449784/is-np-array-num-comparison-very-slow-can-multiprocessing
# -be-used-to-accelera
# Numpy to Tuples - https://stackoverflow.com/questions/10016352/convert-numpy-array-to-tuple
# Numpy Where - https://www.reddit.com/r/learnpython/comments/yo8az5/faster_then_npwhere/
# Row Search - https://stackoverflow.com/questions/14766194/testing-whether-a-numpy-array-contains-a-given-row
# Searching Numpy Arrays - https://stackoverflow.com/questions/3030480/how-do-i-select-elements-of-an-array-given-
# condition
# Sorting Objects - https://stackoverflow.com/questions/403421/how-do-i-sort-a-list-of-objects-based-on-an-attribute-of-
# the-objects
# Sorting on Multiple Columns - https://stackoverflow.com/questions/29352511/numpy-sort-ndarray-on-multiple-columns
# Specific Columns - https://stackoverflow.com/questions/8386675/extracting-specific-columns-in-numpy-array
# Speed-Up Numpy - https://stackoverflow.com/questions/77990840/how-to-speed-up-loops-of-numpy-array-calculations
