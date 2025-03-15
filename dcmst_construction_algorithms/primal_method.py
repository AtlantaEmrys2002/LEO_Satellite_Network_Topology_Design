# Libraries
import numpy as np

# Function constructs initial DCMST by greedily adding the shortest edges that connect vertices not currently within the
# tree to vertices already within the tree. Function returns tree and degree of each vertex in the tree.
def modified_prims_algorithm(cost_matrix, constraints, total_satellites:int, initial_node:int):
    """
    # Function constructs initial DCMST by greedily adding the shortest edges that connect vertices not currently within the tree to vertices already within the tree. However, degree constraints for each node are obeyed. I.e. this is the modified version of Prim's Minimum Spanning Tree Algorithm presented in the original paper on DCMSTs (see report for full reference). Function returns tree and degree of each vertex in the tree.

    :param cost_matrix:
    :param constraints:
    :param total_satellites:
    :param initial_node:
    :return:
    """
    # Initialise list of edges
    tree = np.zeros((total_satellites, total_satellites))

    # All the vertices within the tree - select random initial vertex. Set quicker to search. Initial node chosen randomly.
    tree_vertices = {initial_node}

    # Stores the current degree of all satellites
    degree = np.zeros(total_satellites)

    # Create array of edges and their associated costs - take j in range (k+1, total_satellites) as the matrix is
    # symmetric and reduces search space
    sorted_costs = np.asarray([[cost_matrix[k, j], k, j] for k in range(total_satellites) for j in range(k+1, total_satellites)])

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

            if (((first in tree_vertices) and (second in tree_vertices)) or ((first not in tree_vertices) and (second
                    not in tree_vertices))) or ((degree[second] == constraints[second]) or (degree[first] == constraints[first])):
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
