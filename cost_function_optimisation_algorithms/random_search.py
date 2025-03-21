# Libraries
import numpy as np

def random_search(num_param_sets, num_snapshots, num_snapshots_to_sample, num_satellites, method, constellation_name):

    # Randomly sample sets of parameters (where alpha, beta, and gamma can be random variables in range [0, 1])
    parameter_sets = np.random.rand(num_param_sets, 3)

    # Randomly sample snapshots to analyse for each set of parameters (where satellite IDs are in range (num_snapshots)
    snapshot_ids = np.random.randint(0, num_snapshots, (num_param_sets, num_snapshots_to_sample))

    # Run experiments with given parameters
    for k in range(num_param_sets):
