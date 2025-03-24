# Libraries
import numpy as np
from primal_method import modified_prims_algorithm
import warnings


def prufer_encode(tree):
    """
    Encodes tree as its Prufer Number - uses Cayley's representation (see DCMST comparison paper).
    :param tree:
    :return:
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


def prufer_decode(prufer_number):
    """
    Decodes Prufer number and returns tree encoded in Prufer number as adjacency matrix, as well as the degree of each
    node in the resulting tree.
    :param prufer_number:
    :return:
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


def check_degree(tree_encoding, constraints, num_sat: int):
    """
    Checks that all nodes within tree (encoded as a Prufer number) meet degree constraints. The degree of each node in
    the tree is 1 more than the number of times it appears in the Prufer number encoding.
    :param tree_encoding:
    :param constraints:
    :param num_sat:
    :return:
    """
    degrees = np.bincount(tree_encoding, minlength=num_sat) + 1

    comparison = degrees <= constraints

    return np.array_equal(np.full_like(constraints, True), comparison)


def fitness(chromosome, cost_matrix):
    """
    Calculates the total sum cost of all edges in a given degree-constrained minimum spanning tree.
    :param chromosome:
    :param cost_matrix:
    :return:
    """
    # Decode tree from Prufer number to adjacency matrix
    tree, _ = prufer_decode(chromosome)

    # Select all edges in the tree
    edges = np.argwhere(tree > 0)

    # Sort (smaller node first) and remove duplicates (undirected graph)
    edges = np.unique(np.sort(edges), axis=0)

    # Find costs associated with those edges and sum together
    total_cost = sum([cost_matrix[edge[0], edge[1]] for edge in edges])

    return total_cost


def genetic_algorithm(cost_matrix, constraints, num_sat: int, population_size=30):
    """
    Returns a degree-constrained minimum spanning tree of the network built using a genetic algorithm presented in
    'Comparison of Algorithms for the Degree Constrained Minimum Spanning Tree'. https://doi.org/10.1023/A:1011977126230
    :param cost_matrix:
    :param constraints:
    :param num_sat:
    :param population_size:
    """
    # CALCULATE TERMINATION CONDITION #
    # Justification in original paper - larger degree constraints makes problem significantly easier. However, larger
    # problems in general are harder to solve
    termination_condition = (50 * num_sat) / (np.sum(constraints) / num_sat)

    # Initialise counter for number of iterations
    iteration_count = 0

    # GENERATE INITIAL POPULATION OF CHROMOSOMES #

    # Initialise population with degree constrained minimum spanning trees created by a modified version of Prim's
    # algorithm (also used in the primal method for DCMST construction). This is the only change from the genetic
    # algorithm presented in the original paper (which random initialises DCMSTs). Randomly select initial nodes.
    # Encode DCMSTs in Prufer number representation.

    initial_nodes = np.random.choice(num_sat, size=population_size, replace=False).tolist()

    chromosomes = np.array([prufer_encode(modified_prims_algorithm(cost_matrix, constraints, num_sat, initial_nodes[k]))
                            for k in range(population_size)])

    # Select all unique chromosomes
    chromosomes = np.unique(chromosomes, axis=0)

    # Calculate fitness of each chromosome
    fitness_values = [fitness(chromosome, cost_matrix) for chromosome in chromosomes]

    if len(chromosomes) == 1:
        raise ValueError("Only 1 unique DCMST could be constructed.")

    # Successful evolution count - provides user with understanding of "how much" evolution has taken place
    unsuccessful_evol_count = 0

    previous_chromosomes = chromosomes

    while iteration_count < termination_condition:

        # CROSSOVER METHODS #

        # Use a standard single point method to create a new population of children

        # Randomly select sets of parents
        parents = np.random.choice(len(chromosomes) // 2, (2, 1), replace=False)

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
            crossover_points = np.random.randint(0, num_sat - 2, len(chromosomes) - 1)

            # Randomly select order of crossover (with the best solution) for each chromosome
            coin_flips = np.random.rand(len(chromosomes) - 1)

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
                        chromosomes[chromosome] = np.concatenate((chromosomes[chromosome][crossover_point:],
                                                                  chromosomes[pos][:crossover_point]))
                        # Calculate new fitness
                        fitness_values[chromosome] = fitness(chromosomes[chromosome], cost_matrix)
                    else:
                        chromosomes[chromosome] = np.concatenate((chromosomes[pos][crossover_point:],
                                                                  chromosomes[chromosome][:crossover_point]))
                        # Calculate new fitness
                        fitness_values[chromosome] = fitness(chromosomes[chromosome], cost_matrix)

        iteration_count += 1

        if np.equal(previous_chromosomes, chromosomes):
            unsuccessful_evol_count += 1
            warnings.warn("Unsuccessful evolution (no new solutions found). Percentage of unsuccessful evolutions of "
                          "total number of evolutions: " + str(100 * (unsuccessful_evol_count/termination_condition)) +
                          str("%."))

    # Find the best solution and return corresponding DCMST with degree of each node
    pos = np.argmin(np.asarray(fitness_values))

    return prufer_decode(chromosomes[pos])


# References
# Adding 1 to Specific Indices - https://stackoverflow.com/questions/66315038/add-1-to-numpy-array-from-a-list-of-
# indices
# Comparison of Algorithms for the Degree Constrained Minimum Spanning Tree Paper (see report for full reference)
# Counting Occurrences in numpy - https://stackoverflow.com/questions/28663856/how-do-i-count-the-occurrence-of-a-
# certain-item-in-an-ndarray
# First occurrences - https://stackoverflow.com/questions/432112/is-there-a-numpy-function-to-return-the-first-index-of-
# something-in-an-array
# Prufer Sequence - https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
