# Libraries
import data_handling
import dcmst_construction_algorithms as topology_build
from metrics import propagation_delay
import numpy as np
import networkx as nx

def minimum_delay_topology_design_algorithm(constellation_name, num_snapshots, num_satellites, constraints, constant, method):

    # Stores previous snapshot's topology and average propagation delay
    former_topology = 0
    previous_propagation_delay = 0

    # Generates topology for each snapshot
    for k in range(num_snapshots):

        # Define output filename
        output_filename = constellation_name + "_isls" + str(k) + ".txt"

        # Initialise new topology
        new_topology = np.zeros((num_satellites, num_satellites))

        # Read in distance matrix
        distance_matrix = np.load("./" + constellation_name + "/distance_matrices/dist_matrix_" + str(k) + ".npy")

        # Calculate the shortest path between all satellite pairs in the network
        graph = nx.from_numpy_array(distance_matrix)
        result = nx.shortest_path(graph)

        # Find union of all shortest paths in network
        for i in range(num_satellites):
            for j in range(i+1, num_satellites):

                # Look at all edges on path between satellites i and j and add to new topology
                path = result[i][j]
                for p in range(1, len(path)):
                    new_topology[path[p], path[p-1]] = 1
                    new_topology[path[p - 1], path[p]] = 1

        # Calculate degree of all vertices
        degree = [np.sum(new_topology[k]) for k in range(num_satellites)]

        # Attempt to delete edges to prevent violation of degree constraint
        for v in range(num_satellites):
            if degree[v] > constraints[v]:
                # Find all links where v is endpoint
                connected_satellites = np.argwhere(new_topology[v] > 0)

                # Find associated costs
                links = np.array([[distance_matrix[v, u], v, u] for u in connected_satellites])

                # Sort according to distance (cost)
                links = links[links[:, 0].argsort()]

                # In decreasing order, check if deleting component leads to disconnected graph
                for l in links:
                    new_topology[l[1], l[2]] = 0
                    new_topology[l[2], l[1]] = 0

                    # Convert to format
                    deleted_edge_graph = nx.from_numpy_array(new_topology)

                    # If causes disconnection, add edge back in
                    if nx.is_connected(deleted_edge_graph) is False:
                        new_topology[l[1], l[2]] = 1
                        new_topology[l[2], l[1]] = 1
                    else:
                        # Decrease degree
                        degree[l[1]] -= 1
                        degree[l[0]] -= 1

        # Increase connectivity of topology (adding in edges such that degree constraints not violated)
        new_topology = topology_build.increase_connectivity(new_topology, constraints, degree, distance_matrix, num_satellites)

        # If first snapshot, this is the new topology
        if k == 0:
            former_topology = new_topology
            _, current_prop_delay = propagation_delay(new_topology, distance_matrix, num_satellites)
            previous_propagation_delay = current_prop_delay

            # Write topology to file
            data_handling.write_topology_to_file(output_filename, new_topology, method)

        else:
            _, current_prop_delay = propagation_delay(new_topology, distance_matrix, num_satellites)
            # Assuming that if links are switched, changes are done concurrently, so there is no need to take into account link switch time (it will be constant no matter the number of link switches). If no link switches required, topology is same and will remain the same
            if current_prop_delay < constant * previous_propagation_delay:
                former_topology = new_topology
                previous_propagation_delay = current_prop_delay
                data_handling.write_topology_to_file(output_filename, new_topology, method)

        print("HI")







