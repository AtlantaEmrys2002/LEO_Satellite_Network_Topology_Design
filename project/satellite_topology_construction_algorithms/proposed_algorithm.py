# Libraries
import data_handling
import dcmst_construction_algorithms as topology_build
import numpy as np
import random
import satellite_network_attribute_functions as satnet
import time


# Returns heuristic approximation of degree constrained minimum spanning tree of network, using primal-cut branch
# algorithm (see paper for references)
def dcmst(cost_matrix, constraints, total_satellites, method):

    if method == "primal":

        # Construct DCMST according to Primal Algorithm (see function file for more details)
        tree, degree = topology_build.primal_algorithm(cost_matrix, constraints, total_satellites,
                                                       random.randint(0, total_satellites))

    elif method == "aco":

        # Construct DCMST using ant-based algorithm
        tree, degree = topology_build.ant_colony(cost_matrix, constraints, total_satellites)

    elif method == "ga":
        # Construct DCMST using genetic algorithm
        tree, degree = topology_build.genetic_algorithm(cost_matrix, constraints, total_satellites)

    else:
        raise ValueError("That DCMST algorithm does not exist.")

    return tree, degree


# Builds ISL topology for a single snapshot
def heuristic_topology_design_algorithm_isls(arguments):

    (constellation_name, total_satellites, num_snapshot, degree_constraints, snapshot_id, params,
     output_filename_location, method) = arguments

    # Location of network attributes values (i.e. distance matrices, visibility matrices, etc.)
    network_attributes_location = "./" + constellation_name

    # TIME VISIBILITY MATRIX #

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
    cost_matrix = satnet.cost_function(np.load(network_attributes_location + "/visibility_matrices/visibility_matrix_" +
                                               str(snapshot_id) + ".npy"), time_visibility_matrix,
                                       np.load(network_attributes_location + "/distance_matrices/dist_matrix_" +
                                               str(snapshot_id) + ".npy"),
                                       np.load(network_attributes_location + "/sunlight_matrices/sunlight_matrix_" +
                                               str(snapshot_id) + ".npy"), alpha, beta, gamma, total_satellites)

    # BUILD DEGREE CONSTRAINED MINIMUM SPANNING TREE (HEURISTICALLY) #

    # Using a heuristic algorithm, find the (approx.) degree constrained minimum spanning tree (where weights on edges
    # are costs calculated by cost function).

    tree, current_isl_number = dcmst(cost_matrix, degree_constraints, total_satellites, method)

    # INCREASE CONNECTIVITY #

    # Add edges in increasing order of cost (experiment with decreasing cost) until no longer possible to increase
    # connectivity (and, therefore, reliability/fault tolerance) at expense of energy efficiency
    isls = topology_build.increase_connectivity(tree, degree_constraints, current_isl_number, cost_matrix,
                                                total_satellites)

    # SAVE TOPOLOGY #

    # Convert final topology for given snapshot to correct format and save algorithm results to file.
    data_handling.write_topology_to_file(output_filename_location, isls, snapshot_id)

    return
