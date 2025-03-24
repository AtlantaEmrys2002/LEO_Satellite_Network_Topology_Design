# Libraries
from analysis import measure
from build import heuristic_topology_design_algorithm_isls
import csv
from multiprocessing import Pool
import numpy as np
import os

# Seed Random so results can be reproduced
np.random.seed(42)


def random_search(file_name: str, constellation_name: str, num_snapshots: int, num_param_sets: int, num_sat: int,
                  degree_constraints: list[int], dcmst_method: str):
    """
    Runs a random search optimisation function to find near-optimal values for alpha, beta, and gamma weights (can
    easily be adapted to include more weights), generates the topologies for randomly generated weight sets and saves
    the metrics, along with the best topologies.
    :param file_name: name of TLE file storing satellite information
    :param constellation_name: name of satellite constellation for which a topology is being built
    :param num_snapshots: the number of snapshots over one orbital period for which to generate topology
    :param num_param_sets: the number of randomly generated sets of weights for which to build topologies
    :param num_sat: the number of satellites within the satellite network
    :param degree_constraints: the degree constraints of each satellite within the network
    :param dcmst_method: the method with which a degree-constrained minimum spanning tree is being built
    """
    # Randomly sample sets of parameters (where alpha, beta, and gamma can be random variables in range [0, 1])
    parameter_sets = np.random.rand(num_param_sets, 3)

    # Temporary to store results before they are written to files
    results = []

    # File location - where results stored
    results_location = "./Results/novel/" + constellation_name.lower() + "/" + dcmst_method + ("/optimisation/"
                                                                                               "parameter_set_")

    # For the number of parameter sets that have been randomly generated, build topology utilising these weights and
    # save recorded metrics
    for t in range(num_param_sets):

        # Make directory for parameter sets
        if os.path.isdir(results_location + str(t)) is False:

            # Create directory in which to store random search results
            try:
                os.makedirs(results_location + str(t))
            except OSError:
                print("Directory to store results of random search optimisation could not be created.")

        # Generate arguments for functions - results file path will already exist
        snapshot_arguments = [[file_name, constellation_name, num_sat, num_snapshots, degree_constraints, s,
                               parameter_sets[t], results_location + str(t) + "/isls_" + str(s) + ".txt", dcmst_method]
                              for s in range(num_snapshots)]

        # Run experiments with given parameters - assume that all "network attribute matrices", e.g. distance matrices
        # have already been created
        pool = Pool(processes=os.cpu_count())

        pool.map(heuristic_topology_design_algorithm_isls, snapshot_arguments)

        pool.terminate()

        # Generate results files (metrics)
        results.append(list(parameter_sets[t]) + [measure.measure_dynamic(constellation_name, results_location + str(t),
                                                                          num_sat, num_snapshots)])

    # Write Results to CSV Format - this code was adapted from documentation
    with (open("./Results/novel/" + constellation_name.lower() + "/" + dcmst_method + "/optimisation", 'w', newline='')
          as csvfile):
        fieldnames = ['alpha', 'beta', 'gamma', 'max_latency', 'mean_latency', 'average_hop_count', 'link_churn']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

# References
# CSV Documentation - https://docs.python.org/3/library/csv.html#csv.DictWriter
