# Import Relevant Libraries

import data_handling
import dcmst_construction_algorithms as topology_build
import numpy as np
import os
import random
import satellite_network_attribute_functions as satnet
from scipy.spatial.distance import cdist  # imported due to https://vaghefi.medium.com/fast-distance-calculation-in-
# python-bb2bc9810ea5
from skyfield.api import EarthSatellite, load  # recommended by astropy for calculating information about
# satellites described by TLE (convert to TEME), timescale is only defined once
# (https://rhodesmill.org/skyfield/api-time.html#skyfield.timelib.Time)

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()

# Sun's Ephemeris - used to calculate whether satellite is in sunlight or not
if os.path.isfile('./de421.bsp') is False:
    eph = load('./de421.bsp')

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

# Calculate whether each satellite is within sunlight or not (i.e. vulnerable to damage from solar flares) then add 1 to
# all elements where edge connected to a satellite vulnerable to solar flares (i.e. in sunlight) - this an approximation
# used to determine if satellites are vulnerable to solar flares
def sunlight_function(satellites_in_sun, total_satellites):

    sunlight_matrix = np.zeros((total_satellites, total_satellites))

    # in_sun = np.argwhere(np.asarray(satellites_in_sun)).T[0]
    in_sun = np.argwhere(satellites_in_sun).T[0]

    sunlight_matrix[np.ix_(in_sun, in_sun)] = 1

    return sunlight_matrix

# Returns heuristic approximation of degree constrained minimum spanning tree of network, using primal-cut branch
# algorithm (see paper for references)
def dcmst(cost_matrix, constraints, total_satellites):

    # Construct initial DCST (Degree-Constrained Spanning Tree) using modified version of Prim's algorithm (modified so
    # number of edges incident to any given vertex cannot be greater than constraint (maximum degree) of given vertex)
    # tree, degree = prims_algorithm(cost_matrix, constraints, total_satellites)
    tree, degree = topology_build.modified_prims_algorithm(cost_matrix, constraints, total_satellites, random.randint(0, total_satellites))

    print("TWO")

    # Run second stage where edges are swapped if better connection found

    # tree, degree = edge_exchange(cost_matrix, constraints, total_satellites, tree, degree)

    return tree, degree

# Builds ISL topology for a single snapshot
# def heuristic_topology_design_algorithm_isls(input_file_name, constellation_name, total_satellites, orbit_period, max_comm_dist, degree_constraints, snapshot_id, params, output_filename_isls):
def heuristic_topology_design_algorithm_isls(arguments):

    input_file_name, constellation_name, total_satellites, orbit_period, num_snapshot, max_comm_dist, degree_constraints, snapshot_id, params, output_filename_isls = arguments

    # # Check satellite network has a sufficient number of satellites and orbits
    # if total_satellites < 3:
    #     raise ValueError("Number of satellites must be greater than 3.")
    #
    # # In original Hypatia paper, snapshots of 100ms were utilised - this is continued here (all times are in seconds)
    # snapshot_interval = 0.1
    #
    # # TEMPORARY CHANGE
    # snapshot_interval = 60
    #
    # # The number of snapshots over an orbital period
    # num_snapshot = int(orbit_period/snapshot_interval)
    #
    # print(num_snapshot)
    #
    # # TEMPORARY CHANGE
    # # num_snapshot = 10
    #
    # # Get TLEs-formatted data
    # tles_data = data_handling.read_file(input_file_name)
    #
    # # Convert TLES data to Skyfield EarthSatellite objects - used to convert satellite positions to geocentric
    # # coordinates - all measurements in km
    # earth_satellite_objects = [EarthSatellite(satellite_i[1], satellite_i[2], satellite_i[0], ts) for satellite_i in tles_data]
    #
    # ### DISTANCE MATRICES ###
    #
    # # Calculate distance matrices - set time t (using TDB time format) to indicate point in orbital period when snapshot
    # # taken
    #
    # # Calculate the time (in TDB format) at which each snapshot is taken
    # snapshot_times = [snapshot_time_stamp(snapshot_interval * k) for k in range(num_snapshot)]
    #
    # # Calculate distance matrices for each snapshot
    #
    # # Create directory to store the distance matrix for each snapshot in individual file within directory (can't process
    # # all at once otherwise)
    # # if os.path.isdir("./distance_matrices") is False:
    # if os.path.isdir("./"+constellation_name+"/distance_matrices") is False:
    #
    #     # Create directory in which to store distance matrices
    #     try:
    #         # os.mkdir("./distance_matrices")
    #         os.makedirs("./"+constellation_name+"/distance_matrices")
    #     except OSError:
    #         print("Directory to store distance matrices could not be created.")
    #
    #     # Keep track of file ID
    #     file_id = 0
    #
    #     # Calculate the distance matrix (symmetric) for each snapshot of the network. Distance between satellite and
    #     # itself = 0km.
    #     for k in snapshot_times:
    #
    #         # Calculates position of all satellites in the network at snapshot time k
    #         satellites_at_k = [i.at(k).position.km for i in earth_satellite_objects]
    #
    #         # Calculate distance (Euclidean) between all satellite pairs i and j in the network at snapshot time k
    #         dist_matrix = cdist(satellites_at_k, satellites_at_k, metric='euclidean')
    #
    #         # Save distance matrix to .npy file
    #         # np.save("./distance_matrices/dist_matrix_"+str(file_id) + ".npy", dist_matrix)
    #         np.save("./"+constellation_name+"/distance_matrices/dist_matrix_" + str(file_id) + ".npy", dist_matrix)
    #
    #         # Increment ID counter
    #         file_id += 1
    #
    # ### VISIBILITY AND TIME VISIBILITY MATRICES ###
    #
    # # Calculate visibility and time visibility matrices for all snapshots
    #
    # # Create directory to store the visibility matrix for each snapshot in individual file within directory (can't
    # # process all at once otherwise) - within a visibility matrix, element [i][j] is set to 0 if satellites are not
    # # visible to one another
    # # if os.path.isdir("./visibility_matrices") is False:
    # if os.path.isdir("./"+constellation_name+"/visibility_matrices") is False:
    #
    #     # Create directory in which to store distance matrices
    #     try:
    #         # os.mkdir("./visibility_matrices")
    #         os.makedirs("./"+constellation_name+"/visibility_matrices")
    #     except OSError:
    #         print("Directory to store visibility matrices could not be created.")
    #
    #     # Calculate all distance matrices
    #     for k in range(num_snapshot):
    #
    #         # Calculate visibility matrix for snapshot and save to .npy file - load distance from corresponding distance
    #         # matrix file
    #         # visibility_matrix = visibility_function(np.load("./"+constellation_name+"/distance_matrices/dist_matrix_" + str(k) + ".npy"),
    #         #                                         max_comm_dist, total_satellites)
    #         visibility_matrix = satnet.visibility_function(np.load("./"+constellation_name+"/distance_matrices/dist_matrix_" + str(k) + ".npy"),
    #                                                 max_comm_dist)
    #
    #         # np.save("./visibility_matrices/visibility_matrix_" + str(k) + ".npy", visibility_matrix)
    #         np.save("./"+constellation_name+"/visibility_matrices/visibility_matrix_" + str(k) + ".npy", visibility_matrix)
    #
    #
    # ### SUNLIGHT MATRICES ###
    #
    # # Calculate whether satellites are in sunlight (i.e. vulnerable to solar flares) or on the opposite side of the
    # # Earth
    #
    # # if os.path.isdir("./sunlight_matrices") is False:
    # if os.path.isdir("./"+constellation_name+"/sunlight_matrices") is False:
    #
    #     # Create directory in which to store distance matrices
    #     try:
    #         # os.mkdir("./sunlight_matrices")
    #         os.makedirs("./"+constellation_name+"/sunlight_matrices")
    #     except OSError:
    #         print("Directory to store sunlight matrices could not be created.")
    #
    #     file_id = 0
    #
    #     # Calculate all distance matrices
    #     for k in snapshot_times:
    #
    #         # Calculate whether each satellite is in sunlight or not
    #         satellites_in_sun = [i.at(k).is_sunlit(eph) for i in earth_satellite_objects]
    #
    #         # Update matrix such that element sunlight_matrix[i][j] is set to 1 if i or j is in sunlight and save to
    #         # file
    #         # np.save("./sunlight_matrices/sunlight_matrix_" + str(file_id) + ".npy", sunlight_function(satellites_in_sun,
    #         #                                                                                           total_satellites))
    #         np.save("./"+constellation_name+"/sunlight_matrices/sunlight_matrix_" + str(file_id) + ".npy", sunlight_function(satellites_in_sun,
    #                                                                                                   total_satellites))
    #
    #         file_id += 1

    # Calculate time visibility matrix for current snapshot - need to rearrange order of visibility matrices to get time
    # visibility matrices of other snapshots
    # time_visibility_matrix = time_visibility_function(num_snapshot, total_satellites)
    # time_visibility_matrix = time_visibility_function(num_snapshot, total_satellites, snapshot_id, constellation_name)
    print("HI")

    time_visibility_matrix = satnet.time_visibility_function(num_snapshot, total_satellites, snapshot_id, constellation_name)

    ### COST MATRIX ###

    # Calculates cost matrix for current snapshot

    # Set hyperparameters
    # alpha, beta, gamma = 1, 1, 0.2
    alpha, beta, gamma = params

    # Calculate cost matrix. Calculating for current snapshot, so visibility matrix chosen is at pos
    # 0 in array (same for distance matrix

    cost_matrix = satnet.cost_function(np.load("./"+constellation_name+"/visibility_matrices/visibility_matrix_"+str(snapshot_id)+".npy"), time_visibility_matrix,
                                np.load("./"+constellation_name+"/distance_matrices/dist_matrix_"+str(snapshot_id)+".npy"),
                                np.load("./"+constellation_name+"/sunlight_matrices/sunlight_matrix_"+str(snapshot_id)+".npy"), alpha, beta, gamma,
                                total_satellites)




    ### BUILD DEGREE CONSTRAINED MINIMUM SPANNING TREE (HEURISTICALLY) ###

    # Using a heuristic algorithm, find the (approx.) degree constrained minimum spanning tree (where weights on edges
    # are costs calculated by cost function).
    # See original paper on DCMST and Prim's Algorithm (https://en.wikipedia.org/wiki/Prim%27s_algorithm). Tree holds a
    # graphical representation of the network topology (1 where an ISL exists between satellites i and j, 0 otherwise).
    # Current ISL number holds the degree of each node in the graph - i.e. the number of active ISLs each satellite
    # possesses

    print("ONE")

    tree, current_isl_number = dcmst(cost_matrix, degree_constraints, total_satellites)

    print("THREE")

    ### INCREASE CONNECTIVITY ###

    # Add edges in increasing order of cost (experiment with decreasing cost) until no longer possible to increase
    # connectivity (and, therefore, reliability/fault tolerance) at expense of energy efficiency
    isls = topology_build.increase_connectivity(tree, degree_constraints, current_isl_number, cost_matrix, total_satellites)

    print("FOUR")

    ### SAVE TOPOLOGY ###

    # Convert final topology for given snapshot to correct format and save algorithm results to file.
    data_handling.write_topology_to_file(output_filename_isls, isls)

    return