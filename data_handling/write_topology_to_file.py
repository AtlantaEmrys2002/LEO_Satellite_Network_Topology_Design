# Libraries
import numpy as np
import os


# Write resulting topology for given snapshot to file
def write_topology_to_file(location: str, topology, snapshot_num: int):
    """
    Writes a given snapshot's topology (list of ISLs) to file in correct format for integration with Hypatia software.

    :param location:
    :param topology:
    :param snapshot_num:
    """

    # Create directory to store the topologies for each snapshot
    if os.path.isdir(location) is False:

        os.makedirs(location)

        # # Create directory in which to store topology
        # try:
        #     os.mkdir(location)
        # except OSError:
        #     raise OSError("Directory to store resulting topologies could not be created.")

    # Select all ISLs within topology - lists edges in the graph
    tree_edges = np.argwhere(topology > 0)

    # Sort edges and remove duplicates (undirected edges)
    tree_edges = np.unique(np.sort(tree_edges), axis=0)

    # Save results to file compatible with Hypatia simulation software
    np.savetxt(location + "/isls_" + str(snapshot_num) + ".txt", tree_edges, fmt='%i %i')
