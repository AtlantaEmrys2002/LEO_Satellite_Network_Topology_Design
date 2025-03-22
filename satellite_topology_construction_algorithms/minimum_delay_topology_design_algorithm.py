# Libraries
import data_handling
import dcmst_construction_algorithms as topology_build
import itertools
from metrics import propagation_delay
import numpy as np
import networkx as nx
from scipy.sparse import csr_array
from scipy.sparse.csgraph import dijkstra, reconstruct_path

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

        # Read in distance matrix
        distance_matrix = np.load("./" + constellation_name + "/distance_matrices/dist_matrix_" + str(k) + ".npy")

        # Ensure no negative values
        distance_matrix = np.where(distance_matrix > 0, distance_matrix, 0)

        # Calculate the shortest path between all satellite pairs in the network
        # graph = nx.from_numpy_array(distance_matrix, create_using=nx.MultiGraph)
        graph = csr_array(distance_matrix)

        # for i in range(num_satellites):
        #     for j in range(i+1, num_satellites):
        #         # Check if there is a path through the network between satellites i and j in topology
        #         graph = nx.from_numpy_array(new_topology, create_using=nx.MultiGraph)
        #         if nx.has_path(graph, i, j) is False:
        #             dist_graph = nx.from_numpy_array(distance_matrix, create_using=nx.MultiGraph)
        #             path = nx.shortest_path(dist_graph, source=i, target=j)
        #             for p in range(1, len(path)):
        #                 new_topology[path[p], path[p - 1]] = 1
        #                 new_topology[path[p - 1], path[p]] = 1
        #                 print('hi')

        # result = dict(nxp.all_pairs_shortest_path(graph))

        # start = time.time()
        # result = dict(nxp.all_pairs_dijkstra_path(graph))
        #
        # print(time.time() - start)

        # print(np.fromiter((itertools.pairwise(np.array([5, 6, 7, 8, 9, ]))), dtype=np.dtype((int, 2))))

        _, result = dijkstra(graph, directed=False, return_predecessors=True)

        # For each node, finds tree of shortest paths that connects node to every other node in graph (ignore paths that
        # cannot be found (at the moment - this problem with MDTD is discussed in the final paper)
        result = np.array([reconstruct_path(csgraph=graph, predecessors=result[t], directed=False).todense() > 0 for t
                           in range(num_satellites)], dtype=np.int32)

        print("stage 1")

        tmp = result[0]

        # Element-wise add all arrays
        for v in range(1, len(result)):
            tmp += result[v]

        # Find where edges exist in sum across topology
        new_topology = (tmp > 0).astype(np.int32)

        print(len(np.argwhere(new_topology == 0)))

        # Find union of all shortest paths in network
        for sat_i in range(num_satellites):
            for sat_j in range(sat_i+1, num_satellites):
                # Look at all edges on path between satellites i and j and add to new topology
                path_edges = np.fromiter(itertools.pairwise(result[sat_i, sat_j]), dtype=np.dtype((int, 2))).T
                new_topology[path_edges[0], path_edges[1]] = 1
                new_topology[path_edges[1], path_edges[0]] = 1

        # # Find union of all shortest paths in network
        # for i in range(num_satellites):
        #     for j in range(i+1, num_satellites):
        #
        #         # Look at all edges on path between satellites i and j and add to new topology
        #         path = result[i, j]
        #
        #         print(path)
        #
        #         for p in range(1, len(path)):
        #             new_topology[path[p], path[p-1]] = 1
        #             new_topology[path[p - 1], path[p]] = 1

        print("stage 2")

        # Calculate degree of all vertices
        degree = [np.sum(new_topology[k]) for k in range(num_satellites)]

        # Attempt to delete edges to prevent violation of degree constraint
        for v in range(num_satellites):
            if degree[v] > constraints[v]:

                # Find all links where v is endpoint
                connected_satellites = np.argwhere(new_topology[v] > 0).flatten()

                # Find associated costs
                links = np.array([[distance_matrix[v, u], v, u] for u in connected_satellites])

                # Sort according to distance (cost) in decreasing order - CHECK
                links = np.flip(links[links[:, 0].argsort()], 0)

                # In decreasing order, check if deleting component leads to disconnected graph
                for link in links:

                    new_topology[int(link[1]), int(link[2])] = 0
                    new_topology[int(link[2]), int(link[1])] = 0

                    # Convert to format
                    deleted_edge_graph = nx.from_numpy_array(new_topology)

                    # If causes disconnection, add edge back in
                    if nx.is_connected(deleted_edge_graph) is False:
                        new_topology[int(link[1]), int(link[2])] = 1
                        new_topology[int(link[2]), int(link[1])] = 1
                    else:
                        # Decrease degree of relevant satellite vertices
                        degree[int(link[1])] -= 1
                        degree[int(link[2])] -= 1

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
