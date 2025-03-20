# Import Relevant Libraries
import data_handling
import dcmst_construction_algorithms as topology_build
import numpy as np
import random
import satellite_network_attribute_functions as satnet


# Returns heuristic approximation of degree constrained minimum spanning tree of network, using primal-cut branch
# algorithm (see paper for references)
def dcmst(cost_matrix, constraints, total_satellites):

    # Construct initial DCMST (Degree-Constrained Spanning Tree) using modified version of Prim's algorithm (modified so
    # number of edges incident to any given vertex cannot be greater than constraint (maximum degree) of given vertex)
    tree, degree = topology_build.modified_prims_algorithm(cost_matrix, constraints, total_satellites,
                                                           random.randint(0, total_satellites))

    # Run second stage where edges are swapped if better connection found
    tree, degree = topology_build.edge_exchange(cost_matrix, constraints, total_satellites, tree, degree)

    return tree, degree


# Builds ISL topology for a single snapshot
def heuristic_topology_design_algorithm_isls(arguments):

    (input_file_name, constellation_name, total_satellites, orbit_period, num_snapshot, max_comm_dist,
     degree_constraints, snapshot_id, params, output_filename_isls, method) = arguments

    # ### TIME VISIBILITY MATRIX ###

    # Calculate time visibility matrix for current snapshot - need to rearrange order of visibility matrices to get time
    # visibility matrices of other snapshots
    time_visibility_matrix = satnet.time_visibility_function(num_snapshot, total_satellites, snapshot_id,
                                                             constellation_name)

    # COST MATRIX #

    # Calculates cost matrix for current snapshot

    # Set hyperparameters
    alpha, beta, gamma = params

    # Calculate cost matrix. Calculating for current snapshot, so visibility matrix chosen is at pos
    # 0 in array (same for distance matrix
    cost_matrix = satnet.cost_function(np.load("./"+constellation_name+"/visibility_matrices/visibility_matrix_" +
                                               str(snapshot_id)+".npy"), time_visibility_matrix,
                                       np.load("./" + constellation_name+"/distance_matrices/dist_matrix_" +
                                               str(snapshot_id)+".npy"),
                                       np.load("./" + constellation_name+"/sunlight_matrices/sunlight_matrix_" +
                                               str(snapshot_id)+".npy"), alpha, beta, gamma, total_satellites)

    # BUILD DEGREE CONSTRAINED MINIMUM SPANNING TREE (HEURISTICALLY) #

    # Using a heuristic algorithm, find the (approx.) degree constrained minimum spanning tree (where weights on edges
    # are costs calculated by cost function).
    # See original paper on DCMST and Prim's Algorithm (https://en.wikipedia.org/wiki/Prim%27s_algorithm). Tree holds a
    # graphical representation of the network topology (1 where an ISL exists between satellites i and j, 0 otherwise).
    # Current ISL number holds the degree of each node in the graph - i.e. the number of active ISLs each satellite
    # possesses

    tree, current_isl_number = dcmst(cost_matrix, degree_constraints, total_satellites)

    # INCREASE CONNECTIVITY #

    # Add edges in increasing order of cost (experiment with decreasing cost) until no longer possible to increase
    # connectivity (and, therefore, reliability/fault tolerance) at expense of energy efficiency
    isls = topology_build.increase_connectivity(tree, degree_constraints, current_isl_number, cost_matrix,
                                                total_satellites)

    # SAVE TOPOLOGY #

    # Convert final topology for given snapshot to correct format and save algorithm results to file.
    data_handling.write_topology_to_file(output_filename_isls, isls, method)

    return
