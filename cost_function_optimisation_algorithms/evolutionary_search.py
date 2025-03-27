# Libraries
# from build import heuristic_topology_design_algorithm_isls
from satellite_topology_construction_algorithms import heuristic_topology_design_algorithm_isls
from data_handling import write_optimisation_results_to_csv
import os
from analysis.measure import measure_dynamic
from multiprocessing import Pool
import numpy as np


def in_0_1(parameter_set):
    """
    Determines if alpha, beta, and gamma all within [0, 1] interval
    :param parameter_set:
    :return:
    """
    for i in range(0, 3):
        if parameter_set[i] < 0 or parameter_set > 1:
            return False
    return True


def evolutionary_search(constellation_name: str, num_snapshots: int, num_sat: int, degree_constraints: list[int],
                        dcmst_method: str, output_directory: str, num_iterations: int = 1000, mu: int = 2,
                        pop_size: int = 10, step_size: float = 0.05):
    """
    Runs an evolutionary search optimisation function (based on evolutionary strategy) to find near-optimal values for
    alpha, beta, and gamma weights (can easily be adapted to include more weights), generates the topologies for
    a given network using an evolutionary strategy algorithm and saves the metrics, along with the best topologies.
    :param constellation_name:
    :param num_snapshots:
    :param num_sat:
    :param degree_constraints:
    :param dcmst_method:
    :param output_directory:
    :param num_iterations:
    :param mu:
    :param pop_size:
    :param step_size:
    """

    # N.B. lambda / mu should have no remainder

    # Temporary to store results before they are written to files
    results = []

    if os.path.isdir(output_directory) is False:

        # Create directory in which to store evolutionary optimisation search results
        try:
            os.makedirs(output_directory)
        except OSError:
            print("Directory to store results of evolutionary search optimisation could not be created.")

    # EVOLUTIONARY PARAMETER OPTIMISATION SEARCH #

    current_iteration = 0

    # Initialise population with random values for alpha, beta, and gamma
    candidates = np.random.rand(pop_size, 3)

    # Initialise array in which fitness (according to each metric) is stored
    fitness = np.zeros((pop_size, 3))

    # Evaluate fitness of individuals within initial population
    for c in range(pop_size):

        # Generate arguments for topology build
        snapshot_arguments = [(constellation_name, num_sat, num_snapshots,
                               degree_constraints, snapshot_id, [candidates[c][0], candidates[c][1], candidates[c][2]],
                               output_directory, dcmst_method) for snapshot_id in range(num_snapshots)]

        # Build topologies with given candidate values for alpha, beta, and gamma

        pool = Pool(processes=os.cpu_count())

        pool.map(heuristic_topology_design_algorithm_isls, snapshot_arguments)

        pool.terminate()

        # Calculate the fitness (metrics) of initial population
        _, mean_delay, hop_count, link_churn = measure_dynamic(constellation_name, output_directory, num_sat,
                                                               num_snapshots)

        # Assign calculated fitness to population individuals
        fitness[c, 0] = mean_delay
        fitness[c, 1] = hop_count
        fitness[c, 2] = link_churn

        # Add values to results array
        results.append([candidates[c][0], candidates[c][1], candidates[c][2], mean_delay, hop_count, link_churn])

    # For a set number of iterations, evolve solutions
    while current_iteration < num_iterations:

        # SELECT PARENTS #

        # Select parents (truncation selection - select subset of the best solutions as parents) - divide the number of
        # parents to select between three metrics

        # Best solutions in terms of propagation delay
        parents_1 = np.argsort(fitness.T[0])[:mu//3]

        # Best solutions in terms of hop count
        parents_2 = np.argsort(fitness.T[1])[:mu//3]

        # Best solutions in terms of link churn
        parents_3 = np.argsort(fitness.T[2])[:mu//3]

        # Define parents
        parents = np.concatenate((np.concatenate((parents_1, parents_2), axis=0), parents_3), axis=0)

        # CREATE CHILDREN #

        children = []

        # Number of children per parent
        num_children = pop_size // mu

        # Generate viable parameter sets
        for k in range(len(parents)):
            for _ in range(num_children):
                while True:
                    # Generate possible child - selected randomly from Gaussian distribution with mu equal to current
                    # values of alpha, beta, and gamma. Sigma is derived from step_size (user-assigned). This is the
                    # mutation section of the Evolutionary Strategy
                    child = np.random.normal(parents[k], step_size, 3)
                    if in_0_1(child):
                        children.append(child)
                        break

        # EVALUATE CHILDREN #

        # Calculate fitness from topology of each child
        child_fitness = []

        for child in children:

            # BUILD TOPOLOGY WITH CHILD VALUES #

            # Generate arguments for topology build
            snapshot_arguments = [
                (constellation_name, num_sat, num_snapshots, degree_constraints, snapshot_id, [child[0], child[1],
                 child[2]], output_directory, dcmst_method) for snapshot_id in range(num_snapshots)]

            # Build topologies with given candidate values for alpha, beta, and gamma

            pool = Pool(processes=os.cpu_count())

            pool.map(heuristic_topology_design_algorithm_isls, snapshot_arguments)

            pool.terminate()

            # Calculate fitness values for topology returned by child
            _, mean_delay, hop_count, link_churn = measure_dynamic(constellation_name, output_directory, num_sat,
                                                                   num_snapshots)

            child_fitness.append([mean_delay, hop_count, link_churn])

        # Evaluate fitness of individuals within initial population

        for c in range(len(child_fitness)):
            # Add values to results array
            results.append([candidates[c][0], candidates[c][1], candidates[c][2], child_fitness[c][0],
                            child_fitness[c][1], child_fitness[c][2]])

        # Prepare for next iteration
        candidates = children
        fitness = child_fitness

        # Move to next iteration
        current_iteration += 1

    # WRITE RESULTS TO CSV #

    write_optimisation_results_to_csv(output_directory, "novel", results)

# References:
# Evolutionary Search - https://en.wikipedia.org/wiki/Evolutionary_algorithm#Monte-Carlo_methods
# Evolutionary Strategy - https://machinelearningmastery.com/evolution-strategies-from-scratch-in-python/
# Evolutionary Strategy - https://en.wikipedia.org/wiki/Evolution_strategy
