# Libraries
import data_handling
import dcmst_construction_algorithms as topology_build
from metrics import propagation_delay
import networkx as nx
import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import dijkstra, reconstruct_path
import warnings


def minimum_delay_topology_design_algorithm(constellation_name: str, num_snapshots: int, num_satellites: int,
                                            constraints: np.ndarray, constant: float = 1.0):
    """
    Satellite ISL topology construction method (originally presented in paper by Lang et al. (found here:
    https://link.springer.com/chapter/10.1007/978-981-15-3442-3_8).
    :param constellation_name: name of satellite network constellation, e.g. Starlink-550
    :param num_snapshots: the number of snapshots of the network over one orbit for which a topology is constructed
    :param num_satellites: the number of satellites within the network
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param constant: factor by which subsequent topology's propagation delay must be quicker for subsequent topology to
     replace current topology (see paper for full explanation)
    """
    # N.B. constant set to 1 as default like original paper - 'A Novel Topology Design Method for Multi-layered Optical
    # Satellite Networks' (see report for full reference)

    # Stores previous snapshot's topology and average propagation delay
    # N.B. Cannot parallelise this function, as relies on results of previous topology calculations
    former_topology = 0

    # Set to value that is guaranteed to be larger than initial propagation delay
    previous_propagation_delay = 1000000

    # Generates topology for each snapshot
    for k in range(num_snapshots):

        # Define output filename
        output_filename = './Results/mdtd/' + constellation_name.lower()

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

        # Convert to networkx Graph
        G = nx.from_scipy_sparse_array(csr_array(new_topology))

        # Attempt to delete edges to prevent violation of degree constraint
        for v in range(num_satellites):

            if G.degree[v] > constraints[v]:

                # Find all links where v is endpoint
                connected_satellites = G.edges(v)

                # Find associated costs AND only keeps edges were u greater than v (otherwise those edges will have
                # already been considered)

                links = np.array([[distance_matrix[v, u[1]], v, u[1]] for u in connected_satellites])

                # Sort according to distance (cost) in decreasing order
                links = np.flip(links[links[:, 0].argsort()], 0).astype(int)

                # In decreasing order, check if deleting component leads to disconnected graph

                for link in links:

                    a = link[1]
                    b = link[2]

                    G.remove_edge(a, b)

                    if nx.has_path(G, a, b) is False:
                        G.add_edge(a, b, weight=link[0])

        new_topology = nx.to_numpy_array(G)
        degree = np.array(G.degree).T[1]

        # Increase connectivity of topology (adding in edges such that degree constraints not violated)
        new_topology = topology_build.increase_connectivity(new_topology, constraints, degree, distance_matrix,
                                                            num_satellites)

        # If first snapshot, this is the new topology
        if k == 0:
            former_topology = new_topology
            _, current_prop_delay = propagation_delay(new_topology, distance_matrix, num_satellites)
            previous_propagation_delay = current_prop_delay

            # Write topology to file
            data_handling.write_topology_to_file(output_filename, new_topology, k)

        else:
            _, current_prop_delay = propagation_delay(new_topology, distance_matrix, num_satellites)
            # Assuming that if links are switched, changes are done concurrently, so there is no need to take into
            # account link switch time (it will be constant no matter the number of link switches). If no link switches
            # required, topology is same and will remain the same
            if current_prop_delay < constant * previous_propagation_delay:
                former_topology = new_topology
                previous_propagation_delay = current_prop_delay
                data_handling.write_topology_to_file(output_filename, new_topology, k)
            else:
                data_handling.write_topology_to_file(output_filename, former_topology, k)
