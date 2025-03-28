# Libraries
from analysis import measure
from satellite_topology_construction_algorithms import heuristic_topology_design_algorithm_isls
from data_handling import write_optimisation_results_to_csv
from multiprocessing import Pool
import numpy as np
import os

# Seed Random so results can be reproduced
np.random.seed(42)


def random_search(constellation_name: str, num_snapshots: int, num_param_sets: int, num_sat: int,
                  degree_constraints: list[int], dcmst_method: str, output_directory: str):
    """
    Runs a random search optimisation function to find near-optimal values for alpha, beta, and gamma weights (can
    easily be adapted to include more weights), generates the topologies for randomly generated weight sets and saves
    the metrics, along with the best topologies.
    :param constellation_name: name of satellite constellation for which a topology is being built
    :param num_snapshots: the number of snapshots over one orbital period for which to generate topology
    :param num_param_sets: the number of randomly generated sets of weights for which to build topologies
    :param num_sat: the number of satellites within the satellite network
    :param degree_constraints: the degree constraints of each satellite within the network
    :param dcmst_method: the method with which a degree-constrained minimum spanning tree is being built
    :param output_directory:
    """
    # Randomly sample sets of parameters (where alpha, beta, and gamma can be random variables in range [0, 1])
    parameter_sets = np.random.rand(num_param_sets, 3)

    # Temporary to store results before they are written to files
    results = []

    # Make directory for parameter sets
    if os.path.isdir(output_directory) is False:

        # Create directory in which to store random search results
        try:
            os.makedirs(output_directory)
        except OSError:
            print("Directory to store results of random search optimisation could not be created.")

    # For the number of parameter sets that have been randomly generated, build topology utilising these weights and
    # save recorded metrics
    for t in range(num_param_sets):

        # Generate arguments for functions - results file path will already exist

        snapshot_arguments = [[constellation_name, num_sat, num_snapshots, degree_constraints, s, parameter_sets[t],
                               output_directory, dcmst_method] for s in range(num_snapshots)]

        # Run experiments with given parameters - assume that all "network attribute matrices", e.g. distance matrices
        # have already been created
        pool = Pool(processes=os.cpu_count())

        pool.map(heuristic_topology_design_algorithm_isls, snapshot_arguments)

        pool.terminate()

        # with Pool(os.cpu_count()) as pool:
        #     pool.map(heuristic_topology_design_algorithm_isls, snapshot_arguments)

        # Generate results files (metrics)
        results.append(list(parameter_sets[t]) + [measure.measure_dynamic(constellation_name, output_directory,
                                                                          num_sat, num_snapshots)][1:])

        print(results)

    # Write Results to CSV Format - this code was adapted from documentation
    write_optimisation_results_to_csv(output_directory, "novel", results)


# References
# CSV Documentation - https://docs.python.org/3/library/csv.html#csv.DictWriter
