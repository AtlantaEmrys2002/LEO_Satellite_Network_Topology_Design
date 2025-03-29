# Libraries
import copy
from itertools import chain
import numpy as np
import warnings


def prufer_encode(tree: np.ndarray) -> np.ndarray:
    """
    Encodes tree as its Prufer Number - uses Cayley's representation (see DCMST comparison paper).
    :param tree: adjacency matrix representing a spanning tree of the graph
    :return: a Prufer number (encoding) of the given tree
    """
    # Will store Prufer encoding
    P = []

    # Find all edges in tree
    edges = np.argwhere(tree > 0)

    # Sort (smaller node first) and remove duplicates (undirected graph)
    edges = np.unique(np.sort(edges), axis=0)

    while len(np.unique(edges)) > 2:

        # Find leaf nodes (will only appear once)
        leaf = np.argwhere(np.bincount(edges.flatten()) == 1)[0, 0]

        # Find edge that connects leaf to another node
        leaf_location = np.argwhere(np.logical_or(edges[:, 1] == leaf, edges[:, 0] == leaf))[0, 0]

        # Add leaf's neighbour to encoding
        if leaf == edges[leaf_location, 0]:
            P.append(edges[leaf_location, 1])
        else:
            P.append(edges[leaf_location, 0])

        # Delete edge connecting leaf to neighbour from temporary edges
        edges = np.delete(edges, leaf_location, axis=0)

    return np.asarray(P, dtype=np.int32)


def prufer_decode(prufer_number: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Decodes Prufer number and returns tree encoded in Prufer number as adjacency matrix, as well as the degree of each
    node in the resulting tree.
    :param prufer_number: a Prufer number (encoding) of a tree
    :return: an adjacency matrix representing a spanning tree of the graph, as well as the degree of each vertex within
     the tree
    """
    # Tree will have 2 more nodes than length of list
    num_sat = len(prufer_number) + 2

    # Initialise tree (using adjacency matrix representation) as array of all 0s
    tree = np.zeros((num_sat, num_sat))

    # Stores degree of each vertex
    degree = np.ones(num_sat, dtype=np.int32)

    np.add.at(degree, prufer_number, 1)

    for i in prufer_number:

        j = np.nonzero(degree == 1)[0][0]

        tree[i, j] = 1
        tree[j, i] = 1

        degree[i] -= 1
        degree[j] -= 1

    # Final two nodes in tree (not already included)
    degrees_of_1 = np.argwhere(degree == 1)
    u, v = degrees_of_1[0], degrees_of_1[1]

    tree[u, v] = 1
    tree[v, u] = 1

    # Return tree in adjacency matrix format and degree of each node
    return tree.astype(np.int32), np.sum(tree, axis=1).astype(np.int32)


def check_degree(tree_encoding: np.ndarray, constraints: np.ndarray, num_sat: int) -> bool:
    """
    Checks that all nodes within tree (encoded as a Prufer number) meet degree constraints. The degree of each node in
    the tree is 1 more than the number of times it appears in the Prufer number encoding.
    :param tree_encoding: a Prufer number (encoding) of a tree
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param num_sat: the number of satellites within the network
    :return: determines if all degree constraints on vertices are met
    """
    degrees = np.bincount(tree_encoding, minlength=num_sat) + 1

    comparison = degrees <= constraints

    return np.array_equal(np.full_like(constraints, True), comparison)


def fitness(chromosome: np.ndarray, cost_matrix: np.ndarray) -> float:
    """
    Calculates the total sum cost of all edges in a given degree-constrained minimum spanning tree.
    :param chromosome: a list of possible solutions or "genes" to the DCMST problem
    :param cost_matrix: an adjacency matrix, such that element cost_matrix[i][j] represents the cost of the graph edge
     ij
    :return: the fitness/suitability of each solution or "gene" within the chromosome, as determined by the total sum
     cost of all edges within the tree (aim to minimise)
    """
    # Decode tree from Prufer number to adjacency matrix
    tree, _ = prufer_decode(chromosome)

    # Select all edges in the tree
    edges = np.argwhere(tree > 0)

    # Sort (smaller node first) and remove duplicates (undirected graph)
    edges = np.unique(np.sort(edges), axis=0)

    # Find costs associated with those edges and sum together
    total_cost = np.sum(np.asarray([cost_matrix[edge[0], edge[1]] for edge in edges]))

    return total_cost


def random_trees(num_sat: int, constraints: np.ndarray, pop_size: int) -> list:
    """
    Randomly generates a set number (indicated by population size) of valid degree-constrained spanning trees of
    satellite network.
    :param num_sat: the number of satellites within the network
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param pop_size: number of solutions generated each iteration of the algorithm
    :return: list of randomly generated degree constrained minimum spanning trees (encoded as prufer numbers)
    """
    random_chromosomes = []

    while len(np.unique(np.array(random_chromosomes), axis=0)) < pop_size:
        # Randomly select chromosome genes from uniform distribution - discard any trees that violate degree constraints
        # candidate = np.random.randint(0, num_sat, num_sat - 2)

        candidate = np.random.permutation(np.array(list(chain.from_iterable((i, i, i, i) for i in
                                                                            range(num_sat)))))[:num_sat-2]

        if check_degree(candidate, constraints, num_sat) is True:
            random_chromosomes.append(candidate)

    return random_chromosomes


def genetic_algorithm(cost_matrix: np.ndarray, constraints: np.ndarray, num_sat: int, population_size=20):
    """
    Returns a degree-constrained minimum spanning tree of the network built using a genetic algorithm presented in
    'Comparison of Algorithms for the Degree Constrained Minimum Spanning Tree'. https://doi.org/10.1023/A:1011977126230
    :param cost_matrix: an adjacency matrix, such that element cost_matrix[i][j] represents the cost of the graph edge
     ij
    :param constraints: list that describes the maximum number of ISLs each satellite can establish at a given
     point in time
    :param num_sat: the number of satellites within the network
    :param population_size: number of solutions generated each iteration of the algorithm
    :return: a DCMST and the degree of each vertex within the tree
    """
    # CALCULATE TERMINATION CONDITION #
    # Justification in original paper - larger degree constraints makes problem significantly easier. However, larger
    # problems in general are harder to solve
    # termination_condition = (50 * num_sat) / (np.sum(constraints) / num_sat)
    termination_condition = (20 * num_sat) / (np.sum(constraints) / num_sat)

    # Initialise counter for number of iterations
    iteration_count = 0

    # GENERATE INITIAL POPULATION OF CHROMOSOMES #

    # Initialise population with degree constrained minimum spanning trees created by a modified version of Prim's
    # algorithm (also used in the primal method for DCMST construction). This is the only change from the genetic
    # algorithm presented in the original paper (which random initialises DCMSTs). Randomly select initial nodes.
    # Encode DCMSTs in Prufer number representation.

    # Randomly generate valid and unique degree-constrained trees as initial population
    chromosomes = np.array(random_trees(num_sat, constraints, population_size))

    # Calculate fitness of each chromosome
    fitness_values = np.array([fitness(chromosome, cost_matrix) for chromosome in chromosomes])

    if len(chromosomes) == 1:
        raise ValueError("Only 1 unique DCMST could be constructed.")

    # Successful evolution count - provides user with understanding of "how much" evolution has taken place
    unsuccessful_evol_count = 0

    previous_chromosomes = copy.deepcopy(chromosomes)

    while iteration_count < termination_condition:

        # CROSSOVER METHODS #

        # Use a standard single point method to create a new population of children

        # Randomly select sets of parents
        parents = np.random.choice(len(chromosomes), (len(chromosomes)//2, 1, 2), replace=False)
        parents = np.array([parents[k][0] for k in range(len(parents))])

        for parent_pair in parents:

            # Select a point at which to crossover
            crossover_point = np.random.randint(0, num_sat - 2)

            # Create 2 children using simple single-point crossover
            child_1 = np.concatenate((chromosomes[parent_pair[0]][crossover_point:],
                                      chromosomes[parent_pair[1]][:crossover_point]))
            child_2 = np.concatenate((chromosomes[parent_pair[1]][crossover_point:],
                                      chromosomes[parent_pair[0]][:crossover_point]))

            # MUTATION #
            # With probability 0.05, randomly change gene in child chromosome to another randomly chosen value
            mutation_values = np.random.rand(2, 1)
            if mutation_values[0] <= 0.05:
                # Randomly select gene and gene's new value
                gene = np.random.randint(0, num_sat - 2, 1)
                gene_new_val = np.random.randint(0, num_sat, 1)
                child_1[gene] = gene_new_val
            if mutation_values[1] <= 0.05:
                # Randomly select gene and gene's new value
                gene = np.random.randint(0, num_sat - 2, 1)
                gene_new_val = np.random.randint(0, num_sat, 1)
                child_2[gene] = gene_new_val

            # Check children are viable solutions (i.e. they are DCMSTs and do not violate degree constraints) and
            # calculate their associated fitness if they are. Find the "fitter" of the two new solutions.

            child_1_check = check_degree(child_1, constraints, num_sat)
            child_2_check = check_degree(child_2, constraints, num_sat)

            if child_1_check is True and child_2_check is True:
                fitness_child_1 = fitness(child_1, cost_matrix)
                fitness_child_2 = fitness(child_2, cost_matrix)

                if fitness_child_1 > fitness_child_2:
                    best_new_fitness, best_new_solution = fitness_child_2, child_2
                else:
                    best_new_fitness, best_new_solution = fitness_child_1, child_1

            elif child_1_check is True:
                best_new_fitness, best_new_solution = fitness(child_1, cost_matrix), child_1
            else:
                best_new_fitness, best_new_solution = fitness(child_2, cost_matrix), child_2

            # BREEDING #

            if child_1_check is True or child_2_check is True:
                # Find weaker solution of the two parent solutions
                if fitness_values[parent_pair[0]] > fitness_values[parent_pair[1]]:
                    weaker_parent = parent_pair[0]
                else:
                    weaker_parent = parent_pair[1]

                # Check if new solutions are better than current parent solution, if so, then replace
                if fitness_values[weaker_parent] > best_new_fitness:
                    chromosomes[weaker_parent] = best_new_solution
                    fitness_values[weaker_parent] = best_new_fitness

        # Every 10 iterations, the best solution is merged with all other solutions in the population
        if iteration_count % 10 == 0:

            # Select random points at which to cross over each chromosome with best solution
            crossover_points = np.random.randint(0, num_sat - 2, len(chromosomes))

            # Randomly select order of crossover (with the best solution) for each chromosome
            coin_flips = np.random.rand(len(chromosomes))

            # Select best current solution
            pos = np.argmin(np.asarray(fitness_values))

            for chromosome in range(len(chromosomes)):

                # Don't cross over best solution with best solution
                if chromosome == pos:
                    continue
                else:
                    crossover_point = crossover_points[chromosome]
                    # Cross over
                    if coin_flips[pos] > 0.5:

                        potential_new_chromosome = np.concatenate((chromosomes[chromosome][crossover_point:],
                                                                  chromosomes[pos][:crossover_point]))

                        if check_degree(potential_new_chromosome, constraints, num_sat) is True:
                            chromosomes[chromosome] = potential_new_chromosome

                            # Calculate new fitness
                            fitness_values[chromosome] = fitness(chromosomes[chromosome], cost_matrix)
                    else:

                        potential_new_chromosome = np.concatenate((chromosomes[pos][crossover_point:],
                                                                  chromosomes[chromosome][:crossover_point]))

                        if check_degree(potential_new_chromosome, constraints, num_sat) is True:
                            chromosomes[chromosome] = potential_new_chromosome

                            # Calculate new fitness
                            fitness_values[chromosome] = fitness(chromosomes[chromosome], cost_matrix)

        iteration_count += 1

        if np.array_equal(previous_chromosomes, chromosomes):
            unsuccessful_evol_count += 1
            previous_chromosomes = copy.deepcopy(chromosomes)

    percentage_unsuccessful_iterations = 100 * (unsuccessful_evol_count / termination_condition)

    message = ("Percentage of unsuccessful evolutions of total number of evolutions: " +
               str(percentage_unsuccessful_iterations) + "%.")

    if percentage_unsuccessful_iterations > 0:
        warnings.warn(message)

    # Find the best solution and return corresponding DCMST with degree of each node
    pos = np.argmin(np.asarray(fitness_values))

    return prufer_decode(chromosomes[pos])


# References
# Adding 1 to Specific Indices - https://stackoverflow.com/questions/66315038/add-1-to-numpy-array-from-a-list-of-
# indices
# Comparison of Algorithms for the Degree Constrained Minimum Spanning Tree Paper (see report for full reference)
# Counting Occurrences in numpy - https://stackoverflow.com/questions/28663856/how-do-i-count-the-occurrence-of-a-
# certain-item-in-an-ndarray
# Crossover in Evolutionary Algorithms - https://en.wikipedia.org/wiki/Crossover_(evolutionary_algorithm)
# Depth Error in numpy - https://stackoverflow.com/questions/15923081/what-does-valueerror-object-too-deep-for-desired-
# array-mean-and-how-to-fix-it
# First occurrences - https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-
# something-in-an-array
# Generating Random Pairs - https://stackoverflow.com/questions/30890434/how-to-generate-random-pairs-of-numbers-in-
# python-including-pairs-with-one-entr
# Genetic Algorithms - https://en.wikipedia.org/wiki/Genetic_algorithm
# Prufer Sequence - https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
# Raising Warnings - https://stackoverflow.com/questions/3891804/raise-warning-in-python-without-interrupting-program
# Selecting Rows - https://stackoverflow.com/questions/23359886/selecting-rows-in-numpy-ndarray-based-on-the-value-of-
# two-columns
# Sets and Numpy - https://stackoverflow.com/questions/33196102/how-to-turn-numpy-array-to-set-efficiently
