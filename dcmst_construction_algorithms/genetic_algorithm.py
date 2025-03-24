# Libraries
import numpy as np


# Encodes tree as its Prufer Number - uses Cayley's representation (see DCMST comparison paper)
def prufer_encode(tree):

    # Will store Prufer encoding
    P = []

    # Find all edges in tree
    edges = np.argwhere(tree > 0)

    # Sort (smaller node first) and remove duplicates (undirected graph)
    edges = np.unique(np.sort(edges), axis=0)

    while len(np.unique(edges)) > 2:

        # Find leaf nodes (will only appear once)
        leaf = np.argwhere(np.bincount(edges.flatten()) == 1)[0, 0]

        # Find edge that connects leaf to another node
        leaf_location = np.argwhere(np.logical_or(edges[:, 1] == leaf, edges[:, 0] == leaf))[0, 0]

        # Add leaf's neighbour to encoding
        if leaf == edges[leaf_location, 0]:
            P.append(edges[leaf_location, 1])
        else:
            P.append(edges[leaf_location, 0])

        # Delete edge connecting leaf to neighbour from temporary edges
        edges = np.delete(edges, leaf_location, axis=0)

    return np.asarray(P, dtype=np.int32)


# Decodes Prufer number and returns tree encoded in Prufer number as adjacency matrix
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
    return tree.astype(np.int32), np.sum(tree, axis=1).astype(np.int32)


# References
# Adding 1 to Specific Indices - https://stackoverflow.com/questions/66315038/add-1-to-numpy-array-from-a-list-of-
# indices
# Comparison of Algorithms for the Degree Constrained Minimum Spanning Tree Paper (see report for full reference)
# Counting Occurrences in numpy - https://stackoverflow.com/questions/28663856/how-do-i-count-the-occurrence-of-a-
# certain-item-in-an-ndarray
# First occurrences - https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-
# something-in-an-array
# Prufer Sequence - https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
