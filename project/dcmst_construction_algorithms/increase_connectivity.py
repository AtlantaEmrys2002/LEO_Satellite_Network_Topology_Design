# Libraries
import numpy as np


def increase_connectivity(tree: np.ndarray, degree_constraints: np.ndarray, current_isl_number: np.ndarray,
                          cost_matrix: np.ndarray, total_satellites: int) -> np.ndarray:
    """
    Gradually and greedily adds ISLs to graph to increase the connectivity between satellites until no more ISLs can be
    established without breaking constraints.

    :param tree: a degree-constrained minimum spanning tree of the graph
    :param degree_constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param current_isl_number: list that describes the current degree of each satellite within the network
    :param cost_matrix: costs assigned to each edge within the graph that represents the satellite network
    :param total_satellites: the number of satellites within the network
    :return: a network that contains the tree spanning the network, as well as other edges added such that the degree
     constraint of each vertex is not broken
    """
    # Ensure all feasible edges (i.e. ISL can be established between satellites) not in tree are considered
    cost_matrix = np.where(tree == 0, cost_matrix, -1)

    edges = np.argwhere(cost_matrix > 0)

    # Sort edges and remove duplicates (undirected edges)
    edges = np.unique(np.sort(edges), axis=0)

    # Create array of potential edges and their associated costs and sort by cost
    costs = np.asarray([[cost_matrix[edge[0], edge[1]], edge[0], edge[1]] for edge in edges])
    sorted_costs = costs[costs[:, 0].argsort()]

    # Once sorted, costs are no longer needed
    sorted_costs = sorted_costs.T[1:].T.astype(int)

    for k in range(total_satellites):
        # If satellite has already met degree constraint, ignore and move to next satellite
        if current_isl_number[k] == degree_constraints[k]:
            continue
        else:

            # Select all edges incident to vertex k
            potential_edges = sorted_costs[sorted_costs[:, 0] == k, :]

            # If no more potential edges
            if potential_edges.size == 0:
                continue
            else:
                pos = 0
                potential_edges_num = len(potential_edges)
                while (pos < potential_edges_num) and (current_isl_number[k] < degree_constraints[k]):

                    current_potential_a, current_potential_b = potential_edges[pos]

                    # If both satellites don't have maximum degree (maximum number of ISLs established)
                    if (current_isl_number[current_potential_a] < degree_constraints[current_potential_a] and
                            current_isl_number[current_potential_b] < degree_constraints[current_potential_b]):

                        # Add edge (ISL)
                        tree[current_potential_a, current_potential_b] = 1
                        tree[current_potential_b, current_potential_a] = 1

                        # Update the current number of active ISLs each satellite has established
                        current_isl_number[current_potential_b] += 1
                        current_isl_number[current_potential_a] += 1

                    # Select next edge
                    pos += 1

    return tree
