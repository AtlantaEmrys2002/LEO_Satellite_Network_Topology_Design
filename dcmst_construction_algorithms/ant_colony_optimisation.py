# Libraries
from collections import deque
import copy
import networkx as nx
import numpy as np
import random


class Ant:
    def __init__(self, location):
        """
        Used to initialise a member of the Ant class in the ACO DCMST construction algorithm.
        :param location:
        """
        self.location = location
        self.tabu_list = deque([])

    def reset(self, num_sat):
        """
        Used by ACO DCMST construction algorithm to "reset" an ant for the next iteration of the algorithm.
        :param num_sat:
        """
        coin_flip = random.random()
        if coin_flip < 0.5:
            self.location = random.randint(0, num_sat - 1)
        self.tabu_list.clear()


class Edge:
    def __init__(self, max_cost, min_cost, edge_cost, u, v):
        """
        Used to initialise a member of the Edge class in the ACO DCMST construction algorithm.
        :param max_cost:
        :param min_cost:
        :param edge_cost:
        :param u:
        :param v:
        """
        self.u = u
        self.v = v
        self.edge_cost = edge_cost
        self.initPhm = (max_cost - edge_cost) + (max_cost - min_cost) / 3
        self.phm = self.initPhm
        self.nVisited = 0

    def update_pheromone(self, eta, minPhm, maxPhm):
        """
        Used by Edge class in ACO DCMST construction algorithm to update pheromones of the edge.
        :param eta:
        :param minPhm:
        :param maxPhm:
        """
        self.phm = (1 - eta) * self.phm + self.nVisited * self.initPhm
        self.nVisited = 0
        if self.phm > maxPhm:
            self.phm = maxPhm - self.initPhm
        if self.phm < minPhm:
            self.phm = minPhm + self.initPhm

    def enhance(self, gamma):
        """
        Used by ACO DCMST construction algorithm to update pheromones of edges found in the best DCMST found (so far) by
        the algorithm.
        :param gamma:
        """
        self.phm = self.phm * gamma

    def restart(self):
        """
        Used by ACO DCMST construction algorithm to escape local optimum.
        """
        self.phm = self.phm * random.uniform(0.1, 0.3)


def initialise_ants_and_edges(cost_matrix, num_sat):
    """
    Initialises ants and edges (they are associated with other values, such as pheromones) for ACO DCMST construction
    algorithm.
    :param cost_matrix:
    :param num_sat:
    :return:
    """
    # CREATE ANTS #

    # Initialise an ant at each vertex (cardinality will always be the number of satellites in the network)
    ants = [Ant(v) for v in range(num_sat)]

    # CREATE EDGES #

    max_cost = np.max(cost_matrix)
    min_cost = np.min(cost_matrix[cost_matrix > 0])

    graph_edges = np.argwhere(cost_matrix > 0)

    graph_edges = np.unique(np.sort(graph_edges), axis=0)

    edges = [Edge(max_cost, min_cost, cost_matrix[u[0], u[1]], u[0], u[1]) for u in graph_edges]

    # Calculate min and max pheromone each edge can have
    maxPhm = 1000 * (max_cost - min_cost) + (max_cost - min_cost) / 3
    minPhm = (max_cost - min_cost) / 3

    return ants, edges, maxPhm, minPhm


def construct_spanning_tree(edges, constraints, nCandidates, num_sat):
    """
    Modified version of Kruskal's algorithm used to construct a DCMST - used in the ACO algorithm to create a DCMST
    based on pheromones laid by ants and edge costs.
    :param edges:
    :param constraints:
    :param nCandidates:
    :param num_sat:
    :return:
    """
    # Initialise spanning tree and degrees of each vertex
    T_n = []
    degrees = np.array([0 for _ in range(num_sat)])

    # Sort edges in decreasing order according to pheromone level
    edges.sort(key=lambda x: x.phm, reverse=True)

    # Select top nCandidates edges from E
    if len(edges) > nCandidates:
        candidate_edges = copy.deepcopy(edges[:nCandidates])
    else:
        candidate_edges = copy.deepcopy(edges)

    # Sort edges according to increasing edge cost
    candidate_edges.sort(key=lambda x: x.edge_cost)

    G = nx.Graph()

    level_of_candidates = 1

    while len(T_n) < num_sat - 1:

        candidate = candidate_edges.pop(0)
        # Check if degree constraint violated by adding edge to tree
        if degrees[candidate.u] < constraints[candidate.u] and degrees[candidate.v] < constraints[candidate.v]:

            G.add_edge(candidate.u, candidate.v)

            # Check if cycle created (i.e. no longer a tree if edge added)
            if len(list(nx.simple_cycles(G))) != 0:
                G.remove_edge(candidate.u, candidate.v)
            else:
                T_n.append(candidate)
                degrees[candidate.u] += 1
                degrees[candidate.v] += 1

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
                candidate_edges.sort(key=lambda x: x.edge_cost)

    return T_n







# OPTIMISE BY CHANGING SPANNING TREE ALGORITHM - USE MODIFIED PRIM'S INSTEAD OF KRUSKAL







def move_ants(ants, edges, max_steps, update_period, eta, minPhm, maxPhm):

    # The ants explore for a given number of steps
    for s in range(max_steps):
        # Update pheromones of edges after a given number of steps - helps reduce execution time.
        if s % update_period == 0:
            for e in edges:
                e.update_pheromone(eta, minPhm, maxPhm)

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


def solution_fitness(tree):
    """
    Calculates the sum of all edge costs within a given DCMST - used by the ACO DCMST construction algorithm.
    :param tree:
    :return:
    """
    return sum([edge.edge_cost for edge in tree])


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
    best_spanning_tree = construct_spanning_tree(edges, constraints, candidate_set_cardinality, num_sat)

    # Calculate fitness of spanning tree
    best_fitness = solution_fitness(best_spanning_tree)

    # Continues iteratively improving DCMST until either the maximum number of iterations exceeded or no improvements
    # have been made in a set number of iterations
    while i < max_iterations and (i - i_best) < max_iterations_without_improvement:

        # ANT EXPLORATION #

        ants, edges = move_ants(ants, edges, max_steps, update_period, eta, minPhm, maxPhm)

        # Construct new spanning tree based on exploration
        T = construct_spanning_tree(edges, constraints, candidate_set_cardinality, num_sat)

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

    for edge in edges:
        best_spanning_tree_adjacency[edge[0], edge[1]] = 1
        best_spanning_tree_adjacency[edge[1], edge[0]] = 1

    return best_spanning_tree_adjacency, np.sum(best_spanning_tree_adjacency, axis=1).astype(np.int32)

# References
# Array of Objects - https://stackoverflow.com/questions/4877624/numpy-array-of-objects
# Cycle Detection - https://en.wikipedia.org/wiki/Cycle_(graph_theory)
# Deep Copies - https://stackoverflow.com/questions/37593013/deep-copy-of-a-np-array-of-np-array
# Searching Numpy Arrays - https://stackoverflow.com/questions/3030480/how-do-i-select-elements-of-an-array-given-
# condition
# Sorting Objects - https://stackoverflow.com/questions/403421/how-do-i-sort-a-list-of-objects-based-on-an-attribute-of-
# the-objects
