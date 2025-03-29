# Libraries
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_array


def subtree_builder(tree, deleted_edge):
    """
    Finds the two subtrees created by deleting a given edge from a given degree constrained minimum spanning tree.
    :param tree:
    :param deleted_edge:
    :return:
    """
    # Edge deleted from tree
    edge = deleted_edge

    # CREATE TWO SUBTREES CREATED BY EDGE REMOVAL #

    # Tree without edge (i.e. two subtrees with edge removed)
    tree[edge[0], edge[1]] = 0
    tree[edge[1], edge[0]] = 0

    graph = csr_array(tree)

    # Find connected components (i.e. 2 subtrees)
    n, labels = connected_components(graph, directed=False, return_labels=True, connection='strong')

    # If more than 2 subtrees, not connected
    if n > 2:
        raise ValueError("No path - DCMST does not exist.")

    # Find subtrees
    subtree_i = np.argwhere(np.asarray(labels) > 0).flatten()
    subtree_j = np.argwhere(labels == 0).flatten()

    # Restore tree
    tree[edge[0], edge[1]] = 1
    tree[edge[1], edge[0]] = 1

    return subtree_i, subtree_j


# This performs the edge exchange portion of the primal cut algorithm
def edge_exchange(cost_matrix, constraints, total_satellites, tree, degree):

    # List edges in the tree
    tree_edges = np.argwhere(tree > 0)

    # Sort edges and remove duplicates (undirected edges)
    tree_edges = np.unique(np.sort(tree_edges), axis=0)

    # Evaluate each edge - as DCMST, |E| = |V| - 1
    for m in range(total_satellites - 1):

        # Edge in tree
        edge = tree_edges[m]
        cost_of_edge = cost_matrix[edge[0], edge[1]]

        # Construct subtrees created if edge m is deleted
        subtree_i, subtree_j = subtree_builder(tree, tree_edges[m])

        # Look at all edges connecting subtree i to subtree j
        potential_better_edges = np.array(np.meshgrid(subtree_i, subtree_j)).T.reshape(-1, 2)

        # Sort (smaller node first) and remove duplicates (undirected graph)

        potential_better_edges = np.sort(potential_better_edges)

        # Select all edges with a cost greater than zero
        potential_better_edge_costs = cost_matrix[potential_better_edges.T[0], potential_better_edges.T[1]]
        potential_better_edges = np.vstack((potential_better_edge_costs, potential_better_edges.T[0],
                                            potential_better_edges.T[1])).T

        # Sort according to cost in increasing order
        potential_better_edges = potential_better_edges[potential_better_edges[:, 0].argsort()]

        # Remove all costs less than 0
        costs_less_than_zero = np.searchsorted(potential_better_edges.T[0], 0)
        sorted_costs = potential_better_edges[costs_less_than_zero:]

        # Find all edges with cost less than current edge's cost
        costs_less_than_current_cost = np.searchsorted(sorted_costs.T[0], cost_of_edge, side='left')

        # Find all edges with cost equal to current edge's cost
        costs_equal_to_current_cost = np.searchsorted(sorted_costs.T[0], cost_of_edge, side='right')

        sorted_costs_less_than_current = sorted_costs[:costs_less_than_current_cost]

        # Check this line correct
        sorted_costs_equal_to_current = sorted_costs[costs_less_than_current_cost:costs_equal_to_current_cost]

        new_edge = np.array([])

        # Find new edge with smaller cost with incident nodes not at maximum degree
        if sorted_costs_less_than_current.size != 0:
            for k in sorted_costs_less_than_current:
                k = k.astype(int)
                if (degree[k[1]] != constraints[k[1]]) and (degree[k[2]] != constraints[k[2]]):
                    new_edge = k[1:]
                    break

        # If no new edge found, check if an edge with incident nodes of a smaller degree can be found
        if ((new_edge.size == 0) and ((degree[edge[0]] == constraints[edge[0]]) or
                                      (degree[edge[1]] == constraints[edge[1]])) and
                (sorted_costs_equal_to_current.size != 0)):
            for k in sorted_costs_equal_to_current:
                k = k.astype(int)
                if (degree[k[1]] < constraints[k[1]]) and (degree[k[2]] < constraints[k[2]]):
                    new_edge = k[1:]
                    break

        # Update with new edge
        if new_edge.size != 0:

            # Update tree
            tree[edge[0], edge[1]] = 0
            tree[edge[1], edge[0]] = 0

            tree[new_edge[0], new_edge[1]] = 1
            tree[new_edge[1], new_edge[0]] = 1

            # Update degree
            degree[edge[0]] -= 1
            degree[edge[1]] -= 1

            degree[new_edge[0]] += 1
            degree[new_edge[1]] += 1

            # Update list of tree edges
            if new_edge[0] < new_edge[1]:
                tree_edges[m][0], tree_edges[m][1] = new_edge[0], new_edge[1]
            else:
                tree_edges[m][1], tree_edges[m][0] = new_edge[0], new_edge[1]

    return tree, degree


# Function constructs initial DCMST by greedily adding the shortest edges that connect vertices not currently within the
# tree to vertices already within the tree. Function returns tree and degree of each vertex in the tree.
def modified_prims_algorithm(cost_matrix, constraints, total_satellites: int, initial_node: int):
    """
    # Function constructs initial DCMST by greedily adding the shortest edges that connect vertices not currently
    within the tree to vertices already within the tree. However, degree constraints for each node are obeyed. I.e.
    this is the modified version of Prim's Minimum Spanning Tree Algorithm presented in the original paper on DCMST (
    see report for full reference). Function returns tree and degree of each vertex in the tree.

    :param cost_matrix:
    :param constraints:
    :param total_satellites:
    :param initial_node:
    :return:
    """
    # Initialise list of edges
    tree = np.zeros((total_satellites, total_satellites))

    # All the vertices within the tree - select random initial vertex. Set quicker to search. Initial node chosen
    # randomly.
    tree_vertices = {initial_node}

    # Stores the current degree of all satellites
    degree = np.zeros(total_satellites)

    # Create array of edges and their associated costs - take j in range (k+1, total_satellites) as the matrix is
    # symmetric and reduces search space
    sorted_costs = np.asarray([[cost_matrix[k, j], k, j] for k in range(total_satellites) for j in
                               range(k+1, total_satellites)])

    # Sort the costs in increasing order according to cost
    sorted_costs = sorted_costs[sorted_costs[:, 0].argsort()]

    # Ignore all costs < 0
    costs_less_than_zero = np.searchsorted(sorted_costs.T[0], 0)

    sorted_costs = sorted_costs[costs_less_than_zero:]

    # Only need the edges (not the costs) now that they are sorted
    sorted_costs = sorted_costs.T[1:].T.astype(int)

    # Number of potential edges that could be in tree
    potential_edges_num = len(sorted_costs)

    # While vertices not included in tree - i.e. while there does not exist a path between every pair of vertices
    for _ in range(total_satellites - 1):

        # Take the first edge (as sorted by cost) that meets the condition that it connects one vertex in tree to one
        # vertex not in tree and degree constraint of both vertices is not at maximum
        current_pos = 0

        # Find the edge with the lowest cost that satisfies conditions - i.e. edge that connects vertex in tree to
        # vertex not in tree where both vertices have a degree less than their assigned maximum
        while True:
            first = sorted_costs[current_pos, 0]
            second = sorted_costs[current_pos, 1]

            if ((((first in tree_vertices) and (second in tree_vertices)) or
                 ((first not in tree_vertices) and (second not in tree_vertices))) or
                    ((degree[second] == constraints[second]) or (degree[first] == constraints[first]))):
                current_pos += 1
                # If no way to construct DCMST, no way to construct topology
                if current_pos == potential_edges_num:
                    raise AttributeError("A DCMST cannot be constructed.")
            else:
                break

        # Update tree
        tree[first, second] = 1
        tree[second, first] = 1

        # Update list of vertices in tree
        if first not in tree_vertices:
            tree_vertices.add(first)
        else:
            tree_vertices.add(second)

        # # Update degree count
        degree[first] += 1
        degree[second] += 1

    return tree.astype(int), degree


def primal_algorithm(cost_matrix, constraints, total_satellites: int, initial_node):
    """
    Function builds DCMST according to Primal Algorithm presented in original paper on DCMSTs ('Degree-Constrained
    Minimum Spanning Tree' - Narula and Ho - see report for full citation). See original paper on DCMST and Prim's
    Algorithm (https://en.wikipedia.org/wiki/Prim%27s_algorithm). Tree holds a graphical representation of the network
    topology (1 where an ISL exists between satellites i and j, 0 otherwise). Current ISL number holds the degree of
    each node in the graph - i.e. the number of active ISLs each satellite.
    # possesses
    :param cost_matrix:
    :param constraints:
    :param total_satellites:
    :param initial_node:
    :return:
    """
    # Construct initial DCMST (Degree-Constrained Spanning Tree) using modified version of Prim's algorithm
    # (modified so number of edges incident to any given vertex cannot be greater than constraint (maximum degree)
    # of given vertex)

    tree, degree = modified_prims_algorithm(cost_matrix, constraints, total_satellites, initial_node)

    # Exchange edges if better edge found according to given conditions
    tree, degree = edge_exchange(cost_matrix, constraints, total_satellites, tree, degree)

    return tree, degree


# References
# Connected Components - https://www.geeksforgeeks.org/connected-components-in-an-undirected-graph/
# Connected Components in Large Graphs - https://stackoverflow.com/questions/11016256/connected-components-in-a-graph-
# with-100-million-nodes
# csgraph - https://www.geeksforgeeks.org/scipy-csgraph-compressed-sparse-graph/
# Edge Removal - https://www.geeksforgeeks.org/check-removing-given-edge-disconnects-given-graph/
# Prim's Algorithm - https://en.wikipedia.org/wiki/Prim%27s_algorithm
# Comparison of Networkx Functions - https://stackoverflow.com/questions/56726562/is-all-pairs-dijkstra-faster-than-
# multiple-dijkstra-path
