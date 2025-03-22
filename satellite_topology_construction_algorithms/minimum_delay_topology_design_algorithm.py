# Libraries
import data_handling
import dcmst_construction_algorithms as topology_build
from metrics import propagation_delay
import numpy as np
import networkx as nx
from scipy.sparse import csr_array
from scipy.sparse.csgraph import dijkstra, reconstruct_path, depth_first_order
import time
import warnings

# NEED TO CHECK THIS IS CORRECT - MAKE SURE NX LINES WORKING - IF NOT, THEN USE SCIPY NUMBER OF COMPONENTS


def minimum_delay_topology_design_algorithm(constellation_name, num_snapshots, num_satellites, constraints, constant,
                                            method):

    # Stores previous snapshot's topology and average propagation delay
    # N.B. Cannot parallelise this function, as relies on results of previous topology calculations
    former_topology = 0
    previous_propagation_delay = 0

    # Generates topology for each snapshot
    for k in range(num_snapshots):

        # Define output filename
        output_filename = constellation_name + "_isls" + str(k) + ".txt"

        # Initialise new topology
        new_topology = np.zeros((num_satellites, num_satellites))

        # Read in distance and visibility matrix
        distance_matrix = np.load("./" + constellation_name + "/distance_matrices/dist_matrix_" + str(k) + ".npy")
        visibility_matrix = np.load("./" + constellation_name + "/visibility_matrices/visibility_matrix_" + str(k)
                                    + ".npy")

        # Ignore edges where satellites not visible to one another
        distance_matrix = np.where(visibility_matrix == 1, distance_matrix, -1)

        # Ensure no negative values
        distance_matrix = np.where(distance_matrix > 0, distance_matrix, 0)

        # Calculate the shortest path between all satellite pairs in the network
        graph = csr_array(distance_matrix)

        # For each pair of satellites, use Dijkstra to find the shortest path between the pair and add all edges on
        # path to topology
        for sat_i in range(num_satellites):
            _, result = dijkstra(graph, directed=False, indices=sat_i, return_predecessors=True)
            # Find all edges to add to topology
            edges = np.argwhere((reconstruct_path(csgraph=graph, predecessors=result, directed=False).
                                 todense()).astype(np.int32) > 0).T
            new_topology[edges[0], edges[1]] = 1
            new_topology[edges[1], edges[0]] = 1

        # Let user know if the union of shortest paths is all paths
        if len(np.argwhere(new_topology == 0)) == num_satellites:
            warnings.warn("This has included every edge in the network.")

        # Calculate degree of all vertices
        degree = [np.sum(new_topology[k]) for k in range(num_satellites)]

        # Attempt to delete edges to prevent violation of degree constraint
        for v in range(num_satellites):

            start_v = time.time()

            if degree[v] > constraints[v]:

                # Find all links where v is endpoint
                connected_satellites = np.argwhere(new_topology[v] > 0).flatten()

                # Find associated costs AND only keeps edges were u greater than v (otherwise those edges will have
                # already been considered)
                links = np.array([[distance_matrix[v, u], v, u] for u in connected_satellites if u > v])

                # Sort according to distance (cost) in decreasing order - CHECK
                links = np.flip(links[links[:, 0].argsort()], 0).astype(int)

                # In decreasing order, check if deleting component leads to disconnected graph
                for link in links:

                    new_topology[link[1], link[2]] = 0
                    new_topology[link[2], link[1]] = 0

                    # See if graph still connected and if it causes disconnection, add edge back in
                    if len(depth_first_order(csr_array(new_topology), i_start=link[1], directed=False,
                                             return_predecessors=False)) != num_satellites:
                        new_topology[link[1], link[2]] = 1
                        new_topology[link[2], link[1]] = 1
                    else:
                        # Decrease degree of relevant satellite vertices
                        degree[link[1]] -= 1
                        degree[link[2]] -= 1

            print(time.time() - start_v)

        # Increase connectivity of topology (adding in edges such that degree constraints not violated)
        new_topology = topology_build.increase_connectivity(new_topology, constraints, degree, distance_matrix,
                                                            num_satellites)

        # If first snapshot, this is the new topology
        if k == 0:
            former_topology = new_topology
            _, current_prop_delay = propagation_delay(new_topology, distance_matrix, num_satellites)
            previous_propagation_delay = current_prop_delay

            # Write topology to file
            data_handling.write_topology_to_file(output_filename, new_topology, method)

        else:
            _, current_prop_delay = propagation_delay(new_topology, distance_matrix, num_satellites)
            # Assuming that if links are switched, changes are done concurrently, so there is no need to take into
            # account link switch time (it will be constant no matter the number of link switches). If no link switches
            # required, topology is same and will remain the same
            if current_prop_delay < constant * previous_propagation_delay:
                former_topology = new_topology
                previous_propagation_delay = current_prop_delay
                data_handling.write_topology_to_file(output_filename, new_topology, method)
            else:
                data_handling.write_topology_to_file(output_filename, former_topology, method)
