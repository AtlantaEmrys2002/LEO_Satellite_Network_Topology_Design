# Libraries
from collections import deque
import copy
import numpy as np
import random
from scipy.cluster.hierarchy import DisjointSet
import time


class Ant:
    def __init__(self, location: int):
        """
        Used to initialise a member of the Ant class in the ACO DCMST construction algorithm.
        :param location: node at which ant is currently situated
        """
        self.location = location
        self.tabu_list = deque([])

    def reset(self, num_sat: int):
        """
        Used by ACO DCMST construction algorithm to "reset" an ant for the next iteration of the algorithm.
        :param num_sat: the number of satellites within the network
        """
        coin_flip = random.random()
        if coin_flip < 0.5:
            self.location = random.randint(0, num_sat - 1)
        self.tabu_list.clear()


class Edge:
    def __init__(self, max_cost: float, min_cost: float, edge_cost: float, u: int, v: int):
        """
        Used to initialise a member of the Edge class in the ACO DCMST construction algorithm.
        :param max_cost: the maximum cost on any edge in the network
        :param min_cost: the minimum cost of any edge in the network
        :param edge_cost: the cost of the edge
        :param u: one of the nodes incident to edge
        :param v: one of the nodes incident to edge
        """
        self.u = u
        self.v = v
        self.edge_cost = edge_cost
        self.initPhm = (max_cost - edge_cost) + (max_cost - min_cost) / 3
        self.phm = self.initPhm
        self.nVisited = 0

    def update_pheromone(self, eta: float, minPhm: float, maxPhm: float):
        """
        Used by Edge class in ACO DCMST construction algorithm to update pheromones of the edge.
        :param eta: used to determine amount to update edge pheromone
        :param minPhm: minimum pheromone that can be on any edge
        :param maxPhm: maximum pheromone that can be on any edge
        """
        self.phm = (1 - eta) * self.phm + self.nVisited * self.initPhm
        self.nVisited = 0
        if self.phm > maxPhm:
            self.phm = maxPhm - self.initPhm
        if self.phm < minPhm:
            self.phm = minPhm + self.initPhm

    def enhance(self, gamma: float):
        """
        Used by ACO DCMST construction algorithm to update pheromones of edges found in the best DCMST found (so far) by
        the algorithm.
        :param gamma: used to determine amount of pheromone with which to update edges in best DCMST discovered
        """
        self.phm = self.phm * gamma

    def restart(self):
        """
        Used by ACO DCMST construction algorithm to escape local optimum.
        """
        self.phm = self.phm * random.uniform(0.1, 0.3)


def update_edge_pheromone(edges: np.ndarray, eta: float, minPhm: float, maxPhm: float) -> np.ndarray:
    """
    Used ACO DCMST construction algorithm to update pheromones of each potential edge in graph.
    :param edges: list of edges in graph with their respective values of edge cost, initial pheromone deposit, current
     pheromone deposit, and number of ants that have traversed the edge
    :param eta: used to determine amount to update edge pheromone
    :param minPhm: minimum pheromone that can be on any edge
    :param maxPhm: maximum pheromone that can be on any edge
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


def initialise_ants_and_edges(cost_matrix: np.ndarray, num_sat: int) -> tuple[list[Ant], np.ndarray, float, float]:
    """
    Initialises ants and edges (they are associated with other values, such as pheromones) for ACO DCMST construction
    algorithm.
    :param cost_matrix: costs assigned to each edge within the graph that represents the satellite network
    :param num_sat: the number of satellites within the network
    :return:
    """
    # CREATE ANTS #

    # Initialise an ant at each vertex (cardinality will always be the number of satellites in the network)
    ants = [Ant(v) for v in range(num_sat)]

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


def modified_kruskal(edges: np.ndarray, constraints: np.ndarray, nCandidates: int, num_sat: int) -> list[np.ndarray]:
    """
    Modified version of Kruskal's algorithm used to construct a DCMST - used in the ACO algorithm to create a DCMST
    based on pheromones laid by ants and edge costs.
    :param edges: edges within the graph
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param nCandidates: the number of candidate edges to evaluate at a time
    :param num_sat: the number of satellites within the network
    :return:
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

    return T_n


def move_ants(ants: list[Ant], edges: list[Edge], max_steps: int, update_period: float, eta: float, minPhm: float,
              maxPhm: float):
    """
    Moves edges according to given constraints within the graph. Each ant's position is updated accordingly, as is the
    pheromone level of each edge. This is often deemed the 'exploration phase' of the algorithm
    :param ants: the ants used to explore the graph
    :param edges: edges within the graph
    :param max_steps: maximum number of edges each ant can traverse
    :param update_period:
    :param eta:
    :param minPhm: minimum pheromone that can be on any edge
    :param maxPhm: maximum pheromone that can be on any edge
    :return:
    """
    # The ants explore for a given number of steps
    for s in range(max_steps):
        # Update pheromones of edges after a given number of steps - helps reduce execution time.
        if s % update_period == 0:
            edges = update_edge_pheromone(edges, eta, minPhm, maxPhm)

            # GOT TO ABOVE - NEEDS TESTING

            # for e in edges:
            #     e.update_pheromone(eta, minPhm, maxPhm)

        # Move each ant
        for a in ants:

            nAttempts = 0
            moved = False

            # Move between five vertices or until ant cannot move (which over occurs first).
            while moved is False and nAttempts < 5:

                v_1 = a.location

                # Select random edge from edges incident to v_1 with probability proportional to pheromone on edge
                potential_edges = [[edges[edge], edge] for edge in range(len(edges)) if edges[edge].u == v_1 or
                                   edges[edge].v == v_1]

                # Select an edge adjacent to the vertex at which ant is located randomly (but proportional to the
                # pheromone on that edge)
                random_edge = random.choices([e[1] for e in potential_edges], [e[0].phm for e in potential_edges])[0]

                # Find neighbouring vertex
                if edges[random_edge].u == v_1:
                    v_2 = edges[random_edge].v
                else:
                    v_2 = edges[random_edge].u

                # If ant has not visited that vertex recently, move ant to that vertex
                if v_2 not in a.tabu_list:
                    a.tabu_list.append(v_2)
                    a.location = v_2
                    edges[random_edge].nVisited += 1
                    moved = True
                else:
                    nAttempts += 1

    return ants, edges


# def two_edge_replacement(T, num_sat):
#     T_n = T
#     nTries = 0
#     while nTries < num_sat / 2:
#         # Select random edge in T_n
#         e_1 = random.randint(0, len(T_n))
#         e_b = e_1
#         c_b = 0
#         # for e in T
#
#
# def one_edge_replacement():
#     pass


def solution_fitness(tree) -> float:
    """
    Calculates the sum of all edge costs within a given DCMST - used by the ACO DCMST construction algorithm.
    :param tree:
    :return:
    """
    return np.sum(tree[:, 2])


def ant_colony(cost_matrix, constraints, num_sat: int, max_iterations: int = 10000,
               max_iterations_without_improvement: int = 2500, max_steps: int = 75, eta: float = 0.5,
               gamma: float = 1.5, eta_change: float = 0.95, gamma_change: float = 1.05, R: int = 100):
    """
    Algorithm that constructs a DCMST using an ant-based algorithm. Adapted from T. N. Bui, X. Deng, and C. M. Zrncic,
     “An Improved Ant-Based Algorithm for the Degree-Constrained Minimum Spanning Tree Problem,” IEEE Trans. On Evol.
     Computation, vol. 16, no. 2, pp. 266-278, Apr. 2012., doi: 10.1109/TEVC.2011.2125971. Default values for function
     recommended by this original paper.
    :param cost_matrix:
    :param constraints:
    :param num_sat:
    :param max_iterations:
    :param max_iterations_without_improvement:
    :param max_steps:
    :param eta:
    :param gamma:
    :param eta_change:
    :param gamma_change:
    :param R:
    :return:
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
    best_spanning_tree = modified_kruskal(edges, constraints, candidate_set_cardinality, num_sat)

    # Calculate fitness of spanning tree
    best_fitness = solution_fitness(best_spanning_tree)

    # Continues iteratively improving DCMST until either the maximum number of iterations exceeded or no improvements
    # have been made in a set number of iterations
    while i < max_iterations and (i - i_best) < max_iterations_without_improvement:

        # ANT EXPLORATION #

        ants, edges = move_ants(ants, edges, max_steps, update_period, eta, minPhm, maxPhm)

        # Construct new spanning tree based on exploration
        # T = construct_spanning_tree(edges, constraints, candidate_set_cardinality, num_sat)
        T = modified_kruskal(edges, constraints, candidate_set_cardinality, num_sat)

        # REMOVED AS SIGNIFICANTLY INCREASES RUN TIME - IMPLEMENT!

        # LOCAL OPTIMISATION #
        # T = two_edge_replacement(T)
        # T = one_edge_replacement()

        # Compare the fitness of current best solution to new one
        current_solution_fitness = solution_fitness(T)

        # Update spanning tree if better DCMST has been found
        if current_solution_fitness < best_fitness:
            best_spanning_tree = T
            best_fitness = current_solution_fitness
            i_best = i

        # ENHANCE #

        # Enhance edges (lay pheromones) - update pheromones of all edges in the best solution found so far
        for edge in best_spanning_tree:
            # Find edge in main set of edges and update accordingly
            for e in range(len(edges)):
                if edges[e].u == edge.u and edges[e].v == edge.v:
                    edges[e].enhance(gamma)
                    break

        # RESTART #

        # This prevents algorithm getting stuck in local optimum - restart algorithm (i.e. adjust edge pheromones)
        if i - max(i_best, i_restart) > R:
            i_restart = i
            for edge in best_spanning_tree:
                for e in range(len(edges)):
                    if edges[e].u == edge.u and edges[e].v == edge.v:
                        edges[e].restart()
                        break

        i += 1

        # RESET ANTS #

        # In next iteration, reset ants - approximately half (coin toss) of all ants remain in their current locations
        for a in ants:
            a.reset(num_sat)

        # Update gamma and eta
        gamma *= gamma_change
        eta *= eta_change

    edges = [[e.u, e.v] for e in best_spanning_tree]

    best_spanning_tree_adjacency = np.zeros((num_sat, num_sat))

    # for edge in edges:
    for edge in best_spanning_tree:
        best_spanning_tree_adjacency[edge[0], edge[1]] = 1
        best_spanning_tree_adjacency[edge[1], edge[0]] = 1

    return best_spanning_tree_adjacency, np.sum(best_spanning_tree_adjacency, axis=1).astype(np.int32)

# References
# Array of Objects - https://stackoverflow.com/questions/4877624/numpy-array-of-objects
# Cycle Detection - https://en.wikipedia.org/wiki/Cycle_(graph_theory)
# Deep Copies - https://stackoverflow.com/questions/37593013/deep-copy-of-a-np-array-of-np-array
# Kruskal Algorithms - https://en.wikipedia.org/wiki/Kruskal%27s_algorithm#:~:text=Kruskal's%20algorithm%20finds%20a%20
# minimum,will%20not%20form%20a%20cycle.
# Searching Numpy Arrays - https://stackoverflow.com/questions/3030480/how-do-i-select-elements-of-an-array-given-
# condition
# Sorting Objects - https://stackoverflow.com/questions/403421/how-do-i-sort-a-list-of-objects-based-on-an-attribute-of-
# the-objects
