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
        node_degrees = np.bincount(tmp_edges)

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


def prufer_decode(prufer_number, num_sat):

    # Initialise tree (using adjacency matrix representation) as array of all 0s
    tree = np.zeros((num_sat, num_sat))

    # Initialise temporary number as prufer number passed to function
    P_dash = prufer_number

    # S is set of nodes that have already been considered
    S = set([])

    while len(P_dash) > 0:

        # k is the smallest node index not present in P_dash or S
        k = min(set(range(num_sat)) - prufer_number.union(S))

        # Initialise arc between first element of P_dash and k in tree
        tree[P_dash[0], k] = 1
        tree[k, P_dash[0]] = 1

        # Add k to set of considered nodes
        S.add(k)

        # Remove j from Prufer encoding
        j = P_dash.pop(0)

        if len(P_dash) == 0:

            # Nodes not already in tree
            potential_final_node = set(np.argwhere(np.sum(tree, axis=1) == 0).flatten())

            # Connect j to the smallest node index in S union j that is not already in tree
            final_node = min(potential_final_node.union(S.add(j)))

            tree[final_node, j] = 1
            tree[j, final_node] = 1

    # Return tree in adjacency matrix format and degree of each node
    return tree, np.sum(tree, axis=1) == 0
























# References
# Comparison of Algorithms for the Degree Constrained Minimum Spanning Tree Paper (see report for full reference)
# Counting Occurrences in numpy - https://stackoverflow.com/questions/28663856/how-do-i-count-the-occurrence-of-a-
# certain-item-in-an-ndarray
