# Libraries
import data_handling
import dcmst_construction_algorithms as topology_build
import numpy as np
import random
import satellite_network_attribute_functions as satnet


def dcmst(cost_matrix: np.ndarray, constraints: np.ndarray, total_satellites: int, method: str) -> (
        tuple)[np.ndarray, np.ndarray]:
    """
    Returns heuristically constructed degree constrained minimum spanning tree of network, using primal-cut branch
    algorith, ant colony optimisation algorithm, or genetic algorithm.
    :param cost_matrix: an adjacency matrix, such that element cost_matrix[i][j] represents the cost of the graph edge
     ij
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param total_satellites: the number of satellites within the network
    :param method: the method with which to construct the initial degree-constrained minimum spanning tree (either
     'primal', 'aco', or 'ga')
    :return: a DCMST constructed by the appropriate algorithm and the degree of each node in the graph (i.e. the number
     of active ISLs of each satellite in the constellation)
    """
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


def heuristic_topology_design_algorithm_isls(arguments: list):
    """
    Builds ISL topology for a single snapshot utilising a variation of the novel algorithm proposed by this research
    project.
    :param arguments: arguments to be passed to novel topology construction algorithm, which include the name of the
     mega-constellation, number of satellites within the constellation, the number of snapshots over an orbital period
     for which a topology is constructed, the number of ISL terminals for each satellite, the ID number of the snapshot
     for which a topology is constructed, the cost function coefficients (alpha, beta, and gamma), the name of the files
     in which the resulting topologies are saved, and the DCMST construction method (primal, ACO, or GA) utilised during
     the algorithm.
    :return:
    """
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
