# Libraries
import numpy as np
import os


# Write resulting topology for given snapshot to file
def write_topology_to_file(file_name: str, topology, method: str):
    """
    Writes a given snapshot's topology (list of ISLs) to file in correct format for integration with Hypatia software.

    :param file_name:
    :param topology:
    :param method:
    """
    # Create directory to store the topologies for each snapshot
    # if os.path.isdir("./isl_topologies") is False:
    if os.path.isdir("./" + method + "/isl_topologies") is False:

        # Create directory in which to store distance matrices
        try:
            # os.mkdir("./isl_topologies")
            os.mkdir("./" + method + "/isl_topologies")
        except OSError:
            print("Directory to store distance matrices could not be created.")

    # Select all ISLs within topology - lists edges in the graph
    tree_edges = np.argwhere(topology > 0)

    # Sort edges and remove duplicates (undirected edges)
    tree_edges = np.unique(np.sort(tree_edges), axis=0)

    # Save results to file compatible with Hypatia simulation software
    # np.savetxt("./isl_topologies/" + file_name, tree_edges, fmt='%i %i')
    np.savetxt("./" + method + "/isl_topologies/" + file_name, tree_edges, fmt='%i %i')
