# Libraries
import numpy as np


# Encodes tree as its Prufer Number - uses Cayley's representation (see DCMST comparison paper)
def prufer_encode(tree, num_sat):

    P = []

    # Find all edges in tree
    edges = np.argwhere(tree > 0)

    # Sort (smaller node first) and remove duplicates (undirected graph)
    edges = np.unique(np.sort(edges), axis=1)

    # All edges labelled as temporary (i.e. with 0)
    labelled_edges = np.asarray([[0, edge[0], edge[1]] for edge in edges])

    # While there are more than 2 remaining nodes with a temporary label
    while np.sum(labelled_edges.T[0]) < num_sat - 2:

        # Select all edges labelled with temporary
        tmp_edges = labelled_edges[labelled_edges[:, 0].argsort()]

        # Select all temporary edges
        tmp_edges_index = np.searchsorted(tmp_edges.T[0], 0)
        tmp_edges = tmp_edges[:tmp_edges_index]

        # Select all "leaf nodes in temporarily labelled arcs of tree T"
        node_degrees = np.bincount(tmp_edges.flatten())

        # Find node, such that least index in D1
        k = np.searchsorted(node_degrees, 1)

        # Find edge incident to k and add other vertex to which arc is incident to P
        arc_location = np.searchsorted(tmp_edges.T[1], k)
        j = tmp_edges[arc_location, 2]
        P.append(j)

        # Update labels of edges
        labelled_edges[arc_location, 0] = 1

        # Return tree encoded as Prufer number
        return np.asarray(P, dtype=np.int32)


def prufer_encode_2(tree, num_sat):

    # Will store Prufer encoding
    P = []

    # Find all edges in tree
    edges = np.argwhere(tree > 0)

    # Sort (smaller node first) and remove duplicates (undirected graph)
    edges = np.unique(np.sort(edges), axis=0)

    while len(np.unique(edges)) > 2:

        # Find leaf nodes (will only appear once)
        leaf = np.argwhere(np.bincount(edges.flatten()) == 1)[0][0]

        print(leaf)

        # Find edge that connects leaf to another node
        # pos = np.searchsorted(edges.flatten(), leaf)


        if pos % 2 == 0:
            P.append(edges[int(pos // 2)][1])
        else:
            P.append(edges[int((pos // 2) - 1)][0])

        edges = np.delete(edges, int(pos // 2), axis=0)

    P = np.asarray(P, dtype=np.int32)

    print(P)

    return P





# def prufer_decode(prufer_number, num_sat):
#
#     # Initialise tree (using adjacency matrix representation) as array of all 0s
#     tree = np.zeros((num_sat, num_sat))
#
#     # Initialise temporary number as prufer number passed to function
#     P_dash = prufer_number.tolist()
#
#     # S is set of nodes that have already been considered
#     S = set([])
#
#     while len(P_dash) > 0:
#
#         # k is the smallest node index not present in P_dash or S
#         k = min(set(range(num_sat)) - set(P_dash).union(S))
#
#         # Initialise arc between first element of P_dash and k in tree
#         tree[P_dash[0], k] = 1
#         tree[k, P_dash[0]] = 1
#
#         # Add k to set of considered nodes
#         S.add(k)
#
#         # Remove j from Prufer encoding
#         j = P_dash.pop(0)
#
#         if len(P_dash) == 0:
#
#             # Nodes not already in tree
#             potential_final_node = set(np.argwhere(np.sum(tree, axis=1) == 0).flatten())
#
#             if len(potential_final_node) != 0:
#
#                 # Connect j to the smallest node index in S union j that is not already in tree
#                 final_node = min(potential_final_node.union(S.add(j)))
#
#                 tree[final_node, j] = 1
#                 tree[j, final_node] = 1
#
#     print(tree)
#
#     # Return tree in adjacency matrix format and degree of each node
#     return tree.astype(np.int32), np.sum(tree, axis=1) == 0


def prufer_decode(prufer_number):

    # Tree will have 2 more nodes than length of list
    num_sat = len(prufer_number) + 2

    # Initialise tree (using adjacency matrix representation) as array of all 0s
    tree = np.zeros((num_sat, num_sat))

    # Stores degree of each vertex
    degree = np.ones(num_sat, dtype=np.int32)

    np.add.at(degree, prufer_number, 1)

    for i in prufer_number:

        j = np.nonzero(degree == 1)[0][0]

        tree[i, j] = 1
        tree[j, i] = 1

        degree[i] -= 1
        degree[j] -= 1

    # Final two nodes in tree (not already included)
    degrees_of_1 = np.argwhere(degree == 1)
    u, v = degrees_of_1[0], degrees_of_1[1]

    tree[u, v] = 1
    tree[v, u] = 1

    # Return tree in adjacency matrix format and degree of each node
    return tree.astype(np.int32), np.sum(tree, axis=1) == 0























# References
# Adding 1 to Specific Indices - https://stackoverflow.com/questions/66315038/add-1-to-numpy-array-from-a-list-of-
# indices
# Comparison of Algorithms for the Degree Constrained Minimum Spanning Tree Paper (see report for full reference)
# Counting Occurrences in numpy - https://stackoverflow.com/questions/28663856/how-do-i-count-the-occurrence-of-a-
# certain-item-in-an-ndarray
# First occurrences - https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-
# something-in-an-array
# Prufer Sequence - https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
