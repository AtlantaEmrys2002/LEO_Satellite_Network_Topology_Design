# Libraries
from collections import deque
import numpy as np
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_array

def subtree_builder(total_satellites, tree_edges, pos):
    """
    Find the nodes in subtree i and subtree j if subtree i and subtree j are created by the deletion of edge in tree at position m in list of edges in tree.
    :param total_satellites:
    :param tree_edges:
    :param pos:
    :return:
    """

    # Edge deleted from tree
    edge = tree_edges[pos]

    ### CREATE TWO SUBTREES CREATED BY EDGE REMOVAL ###

    # Tree without edge (i.e. two subtrees with edge removed)
    temp_tree_edges = np.delete(tree_edges, pos, axis=0)

    # Identify 2 subtrees created by deleting edge
    subtree_i = {edge[0]}
    subtree_j = {edge[1]}

    current_i = deque([edge[0]])
    current_j = deque([edge[1]])

    tmp_i = current_i.popleft()
    tmp_j = current_j.popleft()

    # Find edges that have a node in a subtree

    while len(subtree_i) + len(subtree_j) != total_satellites:

        potential_new_nodes_for_subtree_j = np.concatenate((temp_tree_edges[temp_tree_edges[:, 0] == tmp_j], temp_tree_edges[temp_tree_edges[:, 1] == tmp_j]), axis=None)

        for j in potential_new_nodes_for_subtree_j:
            current_j.append(j)
            subtree_j.add(j)

        potential_new_nodes_for_subtree_i = np.concatenate((temp_tree_edges[temp_tree_edges[:, 0] == tmp_i], temp_tree_edges[temp_tree_edges[:, 1] == tmp_i]), axis=None)
        for i in potential_new_nodes_for_subtree_i:
            current_i.append(i)
            subtree_i.add(i)

        if len(current_j) == 0:
            subtree_i = set(range(0, total_satellites)) - subtree_j
            break
        elif len(current_i) == 0:
            subtree_j = set(range(0, total_satellites)) - subtree_i
            break
        else:
            tmp_i = current_i.popleft()
            tmp_j = current_j.popleft()

    return subtree_i, subtree_j


# def subtree_builder_2(total_satellites, tree_edges, tree, pos):
#
#     # Edge deleted from tree
#     edge = tree_edges[pos]
#
#     print(edge)
#
#     ### CREATE TWO SUBTREES CREATED BY EDGE REMOVAL ###
#
#     # Tree without edge (i.e. two subtrees with edge removed)
#     tree[edge[0], edge[1]] = 0
#     tree[edge[1], edge[0]] = 0
#
#     subtree_i = {edge[0]}
#     subtree_j = {edge[1]}
#
#     graph = csr_array(tree)
#
#     path_exists_i = dijkstra(csgraph=graph, directed=False, indices=edge[0])
#
#     path_exists_j = dijkstra(csgraph=graph, directed=False, indices=edge[1])
#
#     print(path_exists_j)
#     print(path_exists_i)
#
#     path_exists_i[edge[0]] = 0
#     path_exists_j[edge[1]] = 0
#
#     for sat in range(total_satellites):
#         if path_exists_i[sat] < np.inf:
#             subtree_i.add(path_exists_i[sat])
#         elif path_exists_j[sat] < np.inf:
#             subtree_j.add(path_exists_j[sat])
#         else:
#             raise ValueError("No path - DCMST does not exist.")
#
#
#
#     # subtree_i = set(np.argwhere(path_exists < np.inf).flatten())
#     # subtree_j = set(range(0, total_satellites)) - subtree_i
#
#     # print(subtree_i)
#     # print(subtree_j)
#     return subtree_i, subtree_j







# This performs the edge exchange portion of the primal cut algorithm
# def edge_exchange(cost_matrix, constraints, total_satellites, tree, degree):
#
#     # List edges in the tree
#     tree_edges = np.argwhere(tree > 0)
#
#     # Sort edges and remove duplicates (undirected edges)
#     tree_edges = np.unique(np.sort(tree_edges), axis=0)
#
#     # Evaluate each edge - as DCMST, |E| = |V| - 1
#     for m in range(total_satellites - 1):
#
#         # # Edge in tree
#         edge = tree_edges[m]
#
#         cost_of_edge = cost_matrix[edge[0], edge[1]]
#         #
#         # ### CREATE TWO SUBTREES CREATED BY EDGE REMOVAL ###
#         #
#         # # Tree without edge (i.e. two subtrees with edge removed)
#         # temp_tree_edges = np.delete(tree_edges, m, axis=0)
#         #
#         # # Identify 2 subtrees created by deleting edge
#         # subtree_i = {edge[0]}
#         # subtree_j = {edge[1]}
#         #
#         # current_i = deque([edge[0]])
#         # current_j = deque([edge[1]])
#         #
#         # # While 2 subtrees do not contain all vertices
#         # while len(current_i) != 0 or len(current_j) != 0:
#         #
#         #     # If all remaining vertices in subtree j
#         #     if len(current_i) == 0:
#         #         subtree_j.update(set(range(total_satellites)) - subtree_i - subtree_j)
#         #         break
#         #     # If all remaining vertices in subtree i
#         #     elif len(current_j) == 0:
#         #         subtree_j.update(set(range(total_satellites)) - subtree_i - subtree_j)
#         #         break
#         #     else:
#         #
#         #         # Fetch current_i[0] and current_j[0] from their respective queues
#         #         current_i_first_val = current_i.popleft()
#         #         current_j_first_val = current_j.popleft()
#         #
#         #         # Select all edges where vertex endpoint of edge is connected to current vertex
#         #         next_i = np.append(temp_tree_edges[temp_tree_edges[:,0] == current_i_first_val],
#         #                            temp_tree_edges[temp_tree_edges[:, 1] == current_i_first_val], axis=0)
#         #
#         #         next_j = np.append(temp_tree_edges[temp_tree_edges[:, 0] == current_j_first_val],
#         #                            temp_tree_edges[temp_tree_edges[:, 1] == current_j_first_val], axis=0)
#         #
#         #         # Select all points not in subtree i or j
#         #         next_i_tmp = set(next_i.flatten()) - subtree_i
#         #         next_i = np.fromiter(next_i_tmp, int, len(next_i_tmp))
#         #
#         #         next_j_tmp = set(next_j.flatten()) - subtree_j
#         #         next_j = np.fromiter(next_j_tmp, int, len(next_j_tmp))
#         #
#         #         # Add unexplored vertices to queues and subtrees
#         #         current_i.extend(next_i)
#         #         current_j.extend(next_j)
#         #
#         #         subtree_i.update(next_i)
#         #         subtree_j.update(next_j)
#         #
#
#         # Construct subtrees created if edge m is deleted
#         subtree_i, subtree_j = subtree_builder(total_satellites, tree_edges, m)
#
#         # # Convert sets to numpy arrays
#         subtree_i = np.fromiter(subtree_i, int, len(subtree_i))
#         subtree_j = np.fromiter(subtree_j, int, len(subtree_j))
#
#         ### ANALYSE EDGE COSTS ###
#
#         # Look at all edges connecting subtree i to subtree j
#         potential_better_edges = np.array(np.meshgrid(subtree_i, subtree_j)).T.reshape(-1, 2)
#
#         # Sort (smaller node first) and remove duplicates (undirected graph)
#         potential_better_edges = np.unique(np.sort(potential_better_edges), axis=1)
#
#         # Create list of edges with their associated costs - only take costs and edges where costs >= 0
#
#         # Select all costs for relevant edges - stack them with corresponding edges
#         potential_better_edge_costs = cost_matrix[potential_better_edges.T[0], potential_better_edges.T[1]]
#         tmp = np.vstack((potential_better_edge_costs, potential_better_edges.T[0], potential_better_edges.T[1])).T
#
#         # Sort according to cost in increasing order
#         tmp = tmp[tmp[:, 0].argsort()]
#
#         # Remove all costs less than 0
#         costs_less_than_zero = np.searchsorted(tmp.T[0], 0)
#         sorted_costs = tmp[costs_less_than_zero:]
#
#         # CHANGED THIS
#         costs_less_than_current_cost = np.searchsorted(sorted_costs.T[0], cost_of_edge, side='right')
#         # costs_less_than_current_cost = np.searchsorted(sorted_costs.T[0], cost_of_edge, side='left')
#         sorted_costs = sorted_costs[:costs_less_than_current_cost]
#
#         # sorted_costs = np.setxor1d(sorted_costs[:, 0], tree_edges[:])
#
#         # print("CURRENT")
#         # print(tree_edges)
#         #
#         # print("POTENTIAL:")
#         # print(sorted_costs)
#
#         # Stores potential new edge
#         new_edge = np.array([])
#
#         # If edge exists with cost smaller than or equal to current edge's cost
#         if sorted_costs.size > 1:
#
#             # If there exists edge with smaller cost than current edge
#             pos_of_smaller_cost = np.searchsorted(sorted_costs.T[0], cost_of_edge)
#
#             if pos_of_smaller_cost != 0:
#
#                 # Iterate over all edges with smaller cost than current edge
#                 for x in sorted_costs[:pos_of_smaller_cost].T[1:].T.astype(int):
#                     if degree[x[0]] != constraints[x[0]] and degree[x[1]] != constraints[x[1]]:
#                         new_edge = x
#                         break
#
#             else:
#                 # If degrees of either vertex are at maximum, see if edge with equal cost that does not have max degree
#                 # for one or both vertices
#                 if degree[edge[0]] == constraints[edge[0]] or degree[edge[1]] == constraints[edge[1]]:
#                     for x in sorted_costs.T[1:].T.astype(int):
#                         if degree[x[0]] != constraints[x[0]] and degree[x[1]] != constraints[x[1]]:
#                             new_edge = x
#                             break
#
#             # If new (better) edge has been found, update tree and degree values
#             if new_edge.size > 0:
#
#                 # Update tree
#                 tree[edge[0], edge[1]] = 0
#                 tree[edge[1], edge[0]] = 0
#
#                 tree[new_edge[0], new_edge[1]] = 1
#                 tree[new_edge[1], new_edge[0]] = 1
#
#                 # Update degree
#                 degree[edge[0]] -= 1
#                 degree[edge[1]] -= 1
#
#                 degree[new_edge[0]] += 1
#                 degree[new_edge[1]] += 1
#
#                 # Update list of tree edges
#                 if new_edge[0] < new_edge[1]:
#                     tree_edges[m][0], tree_edges[m][1] = new_edge[0], new_edge[1]
#                 else:
#                     tree_edges[m][1], tree_edges[m][0] = new_edge[0], new_edge[1]
#
#     print(tree)
#
#     return tree, degree

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
        # subtree_i, subtree_j = subtree_builder(total_satellites, tree_edges, m)
        subtree_i, subtree_j = subtree_builder_2(total_satellites, tree_edges, tree, m)

        # # Convert sets to numpy arrays
        subtree_i = np.fromiter(subtree_i, int, len(subtree_i))
        subtree_j = np.fromiter(subtree_j, int, len(subtree_j))

        # Look at all edges connecting subtree i to subtree j
        potential_better_edges = np.array(np.meshgrid(subtree_i, subtree_j)).T.reshape(-1, 2)

        # Sort (smaller node first) and remove duplicates (undirected graph)
        potential_better_edges = np.unique(np.sort(potential_better_edges), axis=1)

        # Select all edges with a cost greater than zero
        potential_better_edge_costs = cost_matrix[potential_better_edges.T[0], potential_better_edges.T[1]]
        potential_better_edges = np.vstack((potential_better_edge_costs, potential_better_edges.T[0], potential_better_edges.T[1])).T

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
                if (degree[k[1]] != constraints[k[1]]) and (degree[k[2]] != constraints[k[2]]):
                    new_edge = k[1:]
                    break

        # If no new edge found, check if can find edge with incident nodes of a smaller degree
        if (new_edge.size == 0) and ((degree[edge[0]] == constraints[edge[0]]) or (degree[edge[1]] == constraints[edge[1]])) and (sorted_costs_equal_to_current.size != 0):
            for k in sorted_costs_equal_to_current:
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

    return tree










print(edge_exchange(np.asarray([[-1, 4, -1, -1, -1, -1, -1, 8, -1],
                                 [4, -1, 8, -1, -1, -1, -1, 11, -1],
                                 [-1, 8, -1, 7, -1, 4, -1, -1, 2],
                                 [-1, -1, 7, -1, 9, 14, -1, -1, -1],
                                 [-1, -1, -1, 9, -1, 10, -1, -1, -1],
                                 [-1, -1, 4, 14, 10, -1, 2, -1, -1],
                                 [-1, -1, -1, -1, -1, 2, -1, 1, 6],
                                 [8, 11, -1, -1, -1, -1, 1, -1, 7],
                                 [-1, -1, 2, -1, -1, -1, 6, 7, -1]]), np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), 9, np.array(
                 [[0, 1, 0, 0, 0, 0, 0, 0, 0],
                  [1, 0, 1, 0, 0, 0, 0, 0, 0],
                  [0, 1, 0, 1, 0, 0, 0, 0, 1],
                  [0, 0, 1, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 1, 0, 1, 0, 0, 0],
                  [0, 0, 0, 0, 1, 0, 1, 0, 0],
                  [0, 0, 0, 0, 0, 1, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 1, 0, 0],
                  [0, 0, 1, 0, 0, 0, 0, 0, 0]]), np.array([1, 2, 3, 1, 2, 3, 2, 1, 1])))


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

# References
# Prim's Algorithm - https://en.wikipedia.org/wiki/Prim%27s_algorithm
