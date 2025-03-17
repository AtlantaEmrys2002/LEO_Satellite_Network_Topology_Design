#TODO
# When running tests, check all data_generation params are correct (including defaults, e.g. does perigee arg change for
# Kuiper?)
# Find maximum transmission dist for Starlink, Kuiper and telesat - 27000 paper sets at 5014 km - also mentioned in code
# Improve visibility by looking at NSGA-III paper and the antenna direction - will need a way to get antenna direction
# CHECK CORRECTNESS OF ALGORITHM
# INCLUDE OTHER HYPERPARAMETERS - BANDWIDTH
# EXPERIMENT WITH INCREASE CONNECTIVITY FUNC - IN FREE OPTICAL SPACE NETWORKS PAPER (BROADBAND NOT SATELLITE), THEY
# CONNECT LARGEST COST EDGES - REDUCES GRAPH DIAMETER AT EXPENSE OF ENERGY EFFICIENCY
# ADDED 1 to TIME VISIBILITY MATRIX IN COST FUNCTION - IS THIS JUST TEMPORARY FIX OR IS IT CORRECT - NEED TO CHECK!!!!
# Randomly select edge from all edges with same cost rather than just selecting first one

# Import Relevant Libraries

from collections import deque
import data_handling
from multiprocessing import Pool
import numpy as np
import os
import random
import satellite_network_attribute_functions as satnet
from scipy.spatial.distance import cdist  # imported due to https://vaghefi.medium.com/fast-distance-calculation-in-
# python-bb2bc9810ea5
from skyfield.api import EarthSatellite, load  # recommended by astropy for calculating information about
# satellites described by TLE (convert to TEME), timescale is only defined once
# (https://rhodesmill.org/skyfield/api-time.html#skyfield.timelib.Time)
import sys
import time
from build import heuristic_topology_design_algorithm_isls

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()

# Sun's Ephemeris - used to calculate whether satellite is in sunlight or not
# if os.path.isfile('./de421.bsp') is False:
eph = load('./de421.bsp')

# Seed Random so results can be reproduced
random.seed(42)

# Returns the maximum transmission distance for satellite in network (values found through research). TEMPORARILY SET TO LARGE VALUES!!!!!!!
def maximum_transmission_distance(name):
    if 'Starlink' in name:
        return 10000
    elif 'Telesat' in name:
        return 10000
    else:  # Kuiper
        return 10000

# Given time in orbit (based on snapshot number) in seconds, convert the time to TDB time format. Calculate time in
# orbit and reformat for satellite. Assume 60 minutes in an hour and 60 seconds in a minute (no leap seconds), etc.
def snapshot_time_stamp(time_stamp):
    hours = 0
    minutes = 0

    if time_stamp // 3600 >= 1:
        hours = time_stamp // 3600
        time_stamp -= (hours * 3600)
    if time_stamp // 60 >= 1:
        minutes = time_stamp // 60
        time_stamp -= (minutes * 60)
    seconds = time_stamp

    return ts.tdb(2000, 1, 1, hours, minutes, seconds)

# NEED TO FINISH IMPLEMENTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# This performs the edge exchange portion of the primal cut algorithm
def edge_exchange(cost_matrix, constraints, total_satellites, tree, degree):

    # List edges in the tree
    tree_edges = np.argwhere(tree > 0)

    # Sort edges and remove duplicates (undirected edges)

    tree_edges = np.unique(np.sort(tree_edges), axis=0)

    # Evaluate each edge
    for m in range(total_satellites - 1):

        # Edge in tree
        edge = tree_edges[m]

        cost_of_edge = cost_matrix[edge[0], edge[1]]

        ### CREATE TWO SUBTREES CREATED BY EDGE REMOVAL ###

        # Tree without edge (i.e. two subtrees with edge removed)
        temp_tree_edges = np.delete(tree_edges, m, axis=0)

        # Identify 2 subtrees created by deleting edge
        subtree_i = set([edge[0]])
        subtree_j = set([edge[1]])

        current_i = deque([edge[0]])
        current_j = deque([edge[1]])

        # While 2 subtrees do not contain all vertices
        while len(current_i) != 0 or len(current_j) != 0:

            # If all remaining vertices in subtree j
            if len(current_i) == 0:
                subtree_j.update(set(range(total_satellites)) - subtree_i - subtree_j)
                break
            # If all remaining vertices in subtree i
            elif len(current_j) == 0:
                subtree_j.update(set(range(total_satellites)) - subtree_i - subtree_j)
                break
            else:

                # Fetch current_i[0] and current_j[0] from their respective queues
                current_i_first_val = current_i.popleft()
                current_j_first_val = current_j.popleft()

                # Select all edges where vertex endpoint of edge is connected to current vertex
                next_i = np.append(temp_tree_edges[temp_tree_edges[:,0] == current_i_first_val],
                                   temp_tree_edges[temp_tree_edges[:, 1] == current_i_first_val], axis=0)

                next_j = np.append(temp_tree_edges[temp_tree_edges[:, 0] == current_j_first_val],
                                   temp_tree_edges[temp_tree_edges[:, 1] == current_j_first_val], axis=0)

                # Select all points not in subtree i or j
                next_i_tmp = set(next_i.flatten()) - subtree_i
                next_i = np.fromiter(next_i_tmp, int, len(next_i_tmp))

                next_j_tmp = set(next_j.flatten()) - subtree_j
                next_j = np.fromiter(next_j_tmp, int, len(next_j_tmp))

                # Add unexplored vertices to queues and subtrees
                current_i.extend(next_i)
                current_j.extend(next_j)

                subtree_i.update(next_i)
                subtree_j.update(next_j)

        # Convert sets to numpy arrays
        subtree_i = np.fromiter(subtree_i, int, len(subtree_i))
        subtree_j = np.fromiter(subtree_j, int, len(subtree_j))

        ### ANALYSE EDGE COSTS ###

        # Look at all edges connecting subtree i to subtree j
        # potential_better_edges = np.fromiter(product(subtree_i, subtree_j), dtype=np.dtype((int, 2)))
        potential_better_edges = np.array(np.meshgrid(subtree_i, subtree_j)).T.reshape(-1, 2)

        # Create list of edges with their associated costs - only take costs and edges where costs >= 0

        # potential_better_edges_costs = np.asarray([[cost_matrix[e[0], e[1]], e[0], e[1]] for e in potential_better_edges if cost_matrix[e[0], e[1]] >= 0])

        # test = np.asarray(
        #     [[cost_matrix[e[0], e[1]], e[0], e[1]] for e in potential_better_edges])

        # Select all costs for relevant edges - stack them with corresponding edges
        potential_better_edge_costs = cost_matrix[potential_better_edges.T[0], potential_better_edges.T[1]]

        # CHANGED SO NO TRANSPOSE
        tmp = np.vstack((potential_better_edge_costs, potential_better_edges.T[0], potential_better_edges.T[1])).T
        # tmp = np.vstack((potential_better_edge_costs, potential_better_edges.T[0], potential_better_edges.T[1]))

        ### EXPERIMENT STARTS
        # TRUE
        # print(np.array_equal(np.argsort(tmp.T[0]), tmp[:, 0].argsort()))

        # FALSE
        # print(np.array_equal(np.argsort(test, axis=0), np.argsort(tmp.T[0])))
        #
        # print(np.argsort(test, axis=0).T[0])
        # print(tmp[:, 0].argsort())

        # Proof of sorting
        # test = np.array([[0.1, 2, 3], [0, 5, 4], [-1, 7, 3]])
        # print(test[test[:, 0].argsort()])
        # print(test[np.argsort(test, axis=0).T[0]])
        ### EXPERIMENT ENDS

        # CHANGED SO NO TRANSPOSE
        # Sort according to cost in increasing order
        tmp = tmp[tmp[:, 0].argsort()]
        # tmp = tmp[tmp[0].argsort()]


        # Remove all costs less than 0
        # CHANGED SO NO TRANSPOSE
        costs_less_than_zero = np.searchsorted(tmp.T[0], 0)
        # costs_less_than_zero = np.searchsorted(tmp[0], 0)

        # CHANGED SO USES TRANSPOSE
        sorted_costs = tmp[costs_less_than_zero:]
        # sorted_costs = (tmp.T[costs_less_than_zero:]).T

        # Sort potential edges costs (in increasing order)
        # sorted_costs_2 = potential_better_edges_costs[np.argsort(potential_better_edges_costs, axis=0).T[0]]

        ### EXPERIMENT STARTS
        # print(sorted_costs_2.shape)
        # print(sorted_costs.shape)
        #
        # for k in range(len(sorted_costs)):
        #     if False in (sorted_costs[k] == sorted_costs_2[k]):
        #         print("FALSE")
        #         print(sorted_costs[k])
        #         print(sorted_costs_2[k])
        ### EXPERIMENT ENDS

        # Ignore all costs < 0 and take all edges with cost less than OR equal to current cost of edge connecting subtree i to
        # subtree j

        # costs_less_than_zero, costs_less_than_current_cost = np.searchsorted(sorted_costs.T[0], [-1, cost_of_edge], side='right')

        # WORKING VERSION BELOW
        costs_less_than_current_cost = np.searchsorted(sorted_costs.T[0], cost_of_edge, side='right')
        # CHANGED SO NO TRANSPOSE
        # costs_less_than_current_cost = np.searchsorted(sorted_costs[0], cost_of_edge, side='right')

        # sorted_costs = sorted_costs[costs_less_than_zero:costs_less_than_current_cost]
        # WORKING VERSION BELOW
        sorted_costs = sorted_costs[:costs_less_than_current_cost]
        # CHANGED SO NO TRANSPOSE
        # sorted_costs = (sorted_costs.T[:costs_less_than_current_cost]).T

        # Stores potential new edge
        new_edge = np.array([])

        # If edge exists with cost smaller than or equal to current edge's cost
        if sorted_costs.size > 1:
            # If there exists edge with smaller cost than current edge

            # CHANGED SO NO TRANSPOSE
            pos_of_smaller_cost = np.searchsorted(sorted_costs.T[0], cost_of_edge)
            # pos_of_smaller_cost = np.searchsorted(sorted_costs[0], cost_of_edge)

            if pos_of_smaller_cost != 0:
                # Iterate over all edges with smaller cost than current edge

                # CHANGED SO NO TRANSPOSE
                for x in sorted_costs[:pos_of_smaller_cost].T[1:].T.astype(int):
                # for x in sorted_costs[:pos_of_smaller_cost][1:].astype(int):
                    if degree[x[0]] != constraints[x[0]] and degree[x[1]] != constraints[x[1]]:
                        new_edge = x
                        break

            else:
                # If degrees of either vertex are at maximum, see if edge with equal cost that does not have max degree
                # for one or both vertices
                if degree[edge[0]] == constraints[edge[0]] or degree[edge[1]] == constraints[edge[1]]:
                    for x in sorted_costs.T[1:].T.astype(int):
                    # for x in sorted_costs[1:].astype(int):
                        if degree[x[0]] != constraints[x[0]] and degree[x[1]] != degree[x[1]]:
                            new_edge = x
                            break

            # If new (better) edge has been found, update tree and degree values
            if new_edge.size > 0:
                # Update tree
                tree[edge[0], edge[1]] = 0
                tree[edge[1], edge[0]] = 0

                tree[new_edge[0], new_edge[1]] = 1
                tree[new_edge[1], new_edge[0]] = 1

                # Update degree
                degree[edge[0]] -= 1
                degree[edge[1]] -= 1

                degree[new_edge[0]] += 1
                degree[new_edge[1]] += 1

    return tree, degree

# Used for testing
# Data from: https://github.com/AtlantaEmrys2002/hypatia/tree/master/paper/satellite_networks_state
# main("starlink-constellation_tles.txt.tmp", "Starlink-550", 72, 22, 53, 15.19, [3, 7])

# main("starlink-constellation_tles.txt.tmp", "Starlink-550", 72, 22, 53, 15.19, [2, 6, 10, 14], [1, 1, 0.2])

start = time.time()
if __name__ == "__main__":

    # Get inputs
    if len(sys.argv) == 1:
        file_name, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree, mean_motion_rev_per_day, snapshot_ids, params, multi_shell = "starlink-constellation_tles.txt.tmp", "Starlink-550", 72, 22, 53, 15.19, [2, 6, 10, 14], [1, 1, 0.2], False
    elif len(sys.argv) == 9:
        file_name, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree, mean_motion_rev_per_day, snapshot_ids, params = sys.argv
        multi_shell = False
    else:
        file_name, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree, mean_motion_rev_per_day, snapshot_ids, params, multi_shell = sys.argv

    # Calculate the number of satellites in the network
    total_sat = num_sats_per_orbit * num_orbits

    # Generate test data using network description from https://github.com/snkas/hypatia/blob/master/satgenpy/tests/test
    # _tles.py
    data_handling.data_generation(file_name, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree,
                                  mean_motion_rev_per_day)

    # Read test data into appropriate data structure (dictionary)
    # data = format_tle_data(file_name)
    data = data_handling.format_tle_data(file_name)

    # Extract description of satellite positions and unique orbits from data
    satellite_data = data["satellites"]

    # Calculate orbital period of network (or maximum orbital period if satellites orbit at different altitudes)
    if multi_shell is False:
        orbital_period = satnet.orbital_period_calculation(satellite_data[0], total_sat)
    else:
        orbital_period = satnet.orbital_period_calculation(satellite_data, total_sat)

    # Find the maximum communication distance between two satellites (may vary as satellite altitudes vary)
    # max_communication_dist_1 = maximum_communication_distance(file_name, total_sat)
    max_communication_dist = satnet.maximum_communication_distance(file_name, total_sat)

    # This is the maximum distance a satellite can establish signal (transmission power) - need to research for Kuiper
    # and StarLink satellites
    max_transmission_dist = maximum_transmission_distance(constellation_name)

    # Find smaller of these two numbers - the maximum distance before Earth is in the way or the max transmission
    # distance (due to satellite power constraints)
    max_communication_dist = min(max_communication_dist, max_transmission_dist)

    # Initialise degree constraint for each satellite - can be changed based on technical specifications of satellites
    satellite_degree_constraints = [3 for _ in range(len(satellite_data))]

    # Check satellite network has a sufficient number of satellites and orbits
    if total_sat < 3:
        raise ValueError("Number of satellites must be greater than 3.")

    # In original Hypatia paper, snapshots of 100ms were utilised - this is continued here (all times are in seconds)
    # snapshot_interval = 0.1

    # TEMPORARY CHANGE
    snapshot_interval = 60

    # The number of snapshots over an orbital period
    num_snapshot = int(orbital_period / snapshot_interval)

    # TEMPORARY CHANGE
    # num_snapshot = 10

    # Generate arguments for functions
    snapshot_arguments = [[file_name, constellation_name, total_sat, orbital_period, num_snapshot, max_communication_dist,
                           satellite_degree_constraints, t, params,
                           "./" + constellation_name + "_isls" + str(t) + ".txt"] for t in snapshot_ids]

    # Get TLEs-formatted data
    tles_data = data_handling.read_file(file_name)

    # Convert TLES data to Skyfield EarthSatellite objects - used to convert satellite positions to geocentric
    # coordinates - all measurements in km
    earth_satellite_objects = [EarthSatellite(satellite_i[1], satellite_i[2], satellite_i[0], ts) for satellite_i in
                               tles_data]

    ### DISTANCE MATRICES ###

    # Calculate distance matrices - set time t (using TDB time format) to indicate point in orbital period when snapshot
    # taken

    # Calculate the time (in TDB format) at which each snapshot is taken
    snapshot_times = [snapshot_time_stamp(snapshot_interval * k) for k in range(num_snapshot)]

    # Calculate distance matrices for each snapshot

    # Create directory to store the distance matrix for each snapshot in individual file within directory (can't process
    # all at once otherwise)
    # if os.path.isdir("./distance_matrices") is False:
    if os.path.isdir("./" + constellation_name + "/distance_matrices") is False:

        # Create directory in which to store distance matrices
        try:
            # os.mkdir("./distance_matrices")
            os.makedirs("./" + constellation_name + "/distance_matrices")
        except OSError:
            print("Directory to store distance matrices could not be created.")

        # Keep track of file ID
        file_id = 0

        # Calculate the distance matrix (symmetric) for each snapshot of the network. Distance between satellite and
        # itself = 0km.
        for k in snapshot_times:
            # Calculates position of all satellites in the network at snapshot time k
            satellites_at_k = [i.at(k).position.km for i in earth_satellite_objects]

            # Calculate distance (Euclidean) between all satellite pairs i and j in the network at snapshot time k
            dist_matrix = cdist(satellites_at_k, satellites_at_k, metric='euclidean')

            # Save distance matrix to .npy file
            # np.save("./distance_matrices/dist_matrix_"+str(file_id) + ".npy", dist_matrix)
            np.save("./" + constellation_name + "/distance_matrices/dist_matrix_" + str(file_id) + ".npy", dist_matrix)

            # Increment ID counter
            file_id += 1

    ### VISIBILITY AND TIME VISIBILITY MATRICES ###

    # Calculate visibility and time visibility matrices for all snapshots

    # Create directory to store the visibility matrix for each snapshot in individual file within directory (can't
    # process all at once otherwise) - within a visibility matrix, element [i][j] is set to 0 if satellites are not
    # visible to one another
    if os.path.isdir("./" + constellation_name + "/visibility_matrices") is False:

        # Create directory in which to store distance matrices
        try:
            os.makedirs("./" + constellation_name + "/visibility_matrices")
        except OSError:
            print("Directory to store visibility matrices could not be created.")

        # Calculate all visibility matrices
        for k in range(num_snapshot):
            # Calculate visibility matrix for snapshot and save to .npy file - load distance from corresponding distance
            # matrix file
            visibility_matrix = satnet.visibility_function(
                np.load("./" + constellation_name + "/distance_matrices/dist_matrix_" + str(k) + ".npy"),
                max_communication_dist)

            np.save("./" + constellation_name + "/visibility_matrices/visibility_matrix_" + str(k) + ".npy",
                    visibility_matrix)

    ### SUNLIGHT MATRICES ###

    # Calculate whether satellites are in sunlight (i.e. vulnerable to solar flares) or on the opposite side of the
    # Earth

    if os.path.isdir("./" + constellation_name + "/sunlight_matrices") is False:

        # Create directory in which to store distance matrices
        try:
            os.makedirs("./" + constellation_name + "/sunlight_matrices")
        except OSError:
            print("Directory to store sunlight matrices could not be created.")

        file_id = 0

        # Calculate all distance matrices
        for k in snapshot_times:
            # Calculate whether each satellite is in sunlight or not
            satellites_in_sun = [i.at(k).is_sunlit(eph) for i in earth_satellite_objects]

            # Update matrix such that element sunlight_matrix[i][j] is set to 1 if i or j is in sunlight and save to
            # file
            np.save("./" + constellation_name + "/sunlight_matrices/sunlight_matrix_" + str(file_id) + ".npy",
                    satnet.sunlight_function(satellites_in_sun,
                                      total_sat))

            file_id += 1

    # Run topology generation algorithm for each specified snapshot

    pool = Pool(processes=os.cpu_count())

    pool.map(heuristic_topology_design_algorithm_isls, snapshot_arguments)

    pool.terminate()

print(time.time() - start)


# References
# Argsort Error - https://stackoverflow.com/questions/53923914/weird-wrong-outpout-of-np-argsort
# Asserting Numpy Equality - https://stackoverflow.com/questions/3302949/best-way-to-assert-for-numpy-array-equality
# Building Project Documentation - https://medium.com/@pratikdomadiya123/build-project-documentation-quickly-with-the-sphinx-python-2a9732b66594
# Combinations of Arrays - https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
# Combinations with Numpy - https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
# Common Elements of Numpy Array - https://stackoverflow.com/questions/44265572/find-common-elements-in-2d-numpy-arrays
# Conditional List Comprehension - https://stackoverflow.com/questions/4260280/if-else-in-a-list-comprehension
# Conversion with Numpy - https://stackoverflow.com/questions/17973507/why-is-converting-a-long-2d-list-to-numpy-array-so-slow
# Conversion with Numpy 2 - https://stackoverflow.com/questions/8466014/how-to-convert-a-python-set-to-a-numpy-array
# Critical Sections - https://stackoverflow.com/questions/29798635/python-locking-critical-section
# Current Directory - https://stackoverflow.com/questions/5137497/find-the-current-directory-and-files-directory
# DCMST Primal Algorithm - https://www-sciencedirect-com.ezphost.dur.ac.uk/science/article/pii/0305054880900222?via%3Dihub
# Deleting Rows - https://stackoverflow.com/questions/3877491/deleting-rows-in-numpy-array
# Docstrings - https://realpython.com/documenting-python-code/
# Docstrings Tutorial - https://www.datacamp.com/tutorial/docstrings-python
# Earth Radius - https://en.wikipedia.org/wiki/Earth_radius
# Empty Numpy Array - https://stackoverflow.com/questions/11295609/how-can-i-check-whether-a-numpy-array-is-empty-or-not
# Empty Numpy Array 2 - https://stackoverflow.com/questions/11295609/how-can-i-check-whether-a-numpy-array-is-empty-or-not
# Fast Calculation of Euclidean Distance - https://vaghefi.medium.com/fast-distance-calculation-in-python-bb2bc9810ea5
# Feature Scaling - https://en.wikipedia.org/wiki/Feature_scaling
# File Existence - https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists-without-exceptions
# File Handling - https://www.w3schools.com/python/python_file_handling.asp
# First Occurrence of Value - https://stackoverflow.com/questions/16243955/numpy-first-occurrence-of-value-greater-than-existing-value
# Float vs Integer Comparisons - https://stackoverflow.com/questions/30100725/why-are-some-float-integer-comparisons-four-times-slower-than-others
# Hypatia - https://github.com/snkas/hypatia/tree/master
# Hyperparameter Tuning Introduction - https://aws.amazon.com/what-is/hyperparameter-tuning/
# Importing - https://stackoverflow.com/questions/4383571/importing-files-from-different-folder
# Import Error - https://stackoverflow.com/questions/338768/python-error-importerror-no-module-named
# __init__ Files - https://stackoverflow.com/questions/448271/what-is-init-py-for
# Internal Imports - https://www.reddit.com/r/learnpython/comments/yusnr0/python_modulespackages_internal_imports/
# Intersecting Arrays - https://stackoverflow.com/questions/25220975/find-the-non-intersecting-values-of-two-arrays
# Iterating over Numpy Columns - https://stackoverflow.com/questions/10148818/numpy-how-to-iterate-over-columns-of-array
# Itertools and Numpy Product - https://stackoverflow.com/questions/28684492/numpy-equivalent-of-itertools-product
# Loop Speed-Up - https://medium.com/@nirmalya.ghosh/13-ways-to-speedup-python-loops-e3ee56cd6b73
# Main Function Call - https://stackoverflow.com/questions/42012140/if-name-main-function-call
# __main__ Use - https://stackoverflow.com/questions/4042905/what-is-main-py
# Multiple Indices - https://stackoverflow.com/questions/14162026/how-to-get-the-values-from-a-numpy-array-using-multiple-indices
# Multiprocessing Errors - https://www.reddit.com/r/learnpython/comments/yo8pab/why_does_multiprocessing_not_work_properly_in/
# Multiprocessing Imap - https://stackoverflow.com/questions/19063238/in-what-situation-do-we-need-to-use-multiprocessing-pool-imap-unordered
# Multiplying Numpy Array with Scalar - https://stackoverflow.com/questions/53485221/numpy-multiply-array-with-scalar
# Multiprocessing Pool Errors - https://stackoverflow.com/questions/41169146/python-multiprocessing-pool-map-seems-to-not-call-function-at-all
# Multiprocessing - https://stackoverflow.com/questions/35675114/multiprocessing-in-python-going-from-a-for-loop-to-multiprocessing-for-loop?rq=3
# Nested For Loop List Comprehension - https://stackoverflow.com/questions/3633140/nested-for-loops-using-list-comprehension
# Numpy Documentation - https://numpy.org/doc/2.2/reference/index.html
# Numpy Equality - https://stackoverflow.com/questions/10580676/comparing-two-numpy-arrays-for-equality-element-wise
# Numpy Linalg Overhead - https://stackoverflow.com/questions/49866638/why-is-numpy-linalg-norm-slow-when-called-many-times-for-small-size-data
# Numpy savetxt Format - https://stackoverflow.com/questions/71499463/how-to-export-numpy-array-without-brackets
# Orbital Distance - https://space.stackexchange.com/questions/27872/how-to-calculate-the-orbital-distance-between-2-satellites-given-the-tles
# Or Order - https://stackoverflow.com/questions/55503078/why-does-the-order-of-statements-in-an-or-condition-matter
# OS Library Commands - https://stackoverflow.com/questions/8933237/how-do-i-check-if-a-directory-exists-in-python
# Prim's Algorithm - https://en.wikipedia.org/wiki/Prim%27s_algorithm
# Pyephem Code - https://github.com/brandon-rhodes/pyephem
# Pyephem Documentation - https://rhodesmill.org/pyephem/quick
# Python Documentation - https://docs.python.org/3/
# Reading Files - https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
# Reading Files Line By Line - https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
# Relative Imports - https://stackoverflow.com/questions/16981921/relative-imports-in-python-3
# Saving Numpy Data to Files - https://stackoverflow.com/questions/28439701/how-to-save-and-load-numpy-array-data-properly
# Scaling vs Normalisation - https://www.kaggle.com/code/alexisbcook/scaling-and-normalization
# SciPy Documentation - https://docs.scipy.org/doc/scipy/reference/index.html
# Selecting Columns - https://stackoverflow.com/questions/22927181/selecting-specific-rows-and-columns-from-numpy-array
# Selecting Numpy Rows - https://stackoverflow.com/questions/58079075/numpy-select-rows-based-on-condition
# Selecting Rows - https://stackoverflow.com/questions/73733370/select-rows-of-numpy-array-based-on-column-values
# Selecting Rows 2 - https://stackoverflow.com/questions/56419519/selecting-rows-of-numpy-array-that-contains-all-the-given-values
# Selecting Rows 3 - https://stackoverflow.com/questions/23359886/selecting-rows-in-numpy-ndarray-based-on-the-value-of-two-columns
# Selecting Specific Numpy Rows - https://stackoverflow.com/questions/22927181/selecting-specific-rows-and-columns-from-numpy-array
# Set vs List Search - https://stackoverflow.com/questions/5993621/fastest-way-to-search-a-list-in-python
# SkyField Documentation - https://rhodesmill.org/skyfield/ & https://rhodesmill.org/skyfield/toc.html
# Skyfield Timings - https://stackoverflow.com/questions/49494082/skyfield-achieve-sgp4-results-with-1-second-periodicity-for-given-time-interval
# Speedup For Loops - https://stackoverflow.com/questions/8770791/how-to-speed-up-python-loop
# Spg4 Import Error - https://stackoverflow.com/questions/72639365/no-module-named-sgp4-but-requirement-already-satisfied-sgp4-in-c
# Sphinx - https://www.sphinx-doc.org/en/master/tutorial/automatic-doc-generation.html
# Sorting by Column - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
# Subset of Rows - https://stackoverflow.com/questions/70341932/select-subset-of-rows-of-numpy-array-based-on-a-selection-of-rows-in-another-arr
# Sys Arguments - https://stackoverflow.com/questions/11853508/how-to-make-sys-argv-arguments-optional
# Take (Numpy) - https://stackoverflow.com/questions/60192041/python-numpy-np-take-with-2-dimensional-array
# Threshold Values with Numpy - https://stackoverflow.com/questions/37973135/numpy-argmin-for-elements-greater-than-a-threshold
# Threshold Zero - https://stackoverflow.com/questions/28430904/set-numpy-array-elements-to-zero-if-they-are-above-a-specific-threshold
# TLE Definitions - https://platform-cdn.leolabs.space/static/files/tle_definition.pdf?7ba94f05897b4ae630a3c5b65be7396c642d9c72
# Type Hinting - https://stackoverflow.com/questions/35673895/type-hinting-annotation-pep-484-for-numpy-ndarray
# Type Hinting with Numpy - https://stackoverflow.com/questions/58877370/pycharm-typehint-for-numpy-array-of-objects
# Type Hinting with Numpy - https://stackoverflow.com/questions/70714087/how-to-typehint-numpy-array-of-floats-and-or-integers
# Typing - https://stackoverflow.com/questions/14113187/how-do-you-set-a-conditional-in-python-based-on-datatypes
# Unique Values from List - https://stackoverflow.com/questions/12897374/get-unique-values-from-a-list-in-python
# World Geodetic System - https://en.wikipedia.org/wiki/World_Geodetic_System#Definition
# Writing Packages - https://stackoverflow.com/questions/15746675/how-to-write-a-python-module-package
# Zero-Size Reduction - https://www.reddit.com/r/learnpython/comments/lhljgh/valueerror_zerosize_array_to_reduction_operation/

# Extra Code Used During Analysis, Optimisation, and Debugging

# Used to check that visibility matrices change over time (confirming theory)
# for k in range(1, num_snapshot):
#     previous = np.load("./visibility_matrices/visibility_matrix_" + str(k) + ".npy")
#     now = np.load("./visibility_matrices/visibility_matrix_" + str(k - 1) + ".npy")
#     if np.array_equal(previous, now) is False:
#         print("yes")
