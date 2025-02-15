#TODO: test data
# unit tests for each function - need to ensure correctness
# check through notes
# this function will be placed in hypatia/satgenpy/satgen/isls
# properly document your code with doxygen or similar
# get orbital period from satellite orbiting at maximum height?
# this function may be called once a snapshot, so check through code for Hypatia - might just need to pass topology for given time
# Check it returns satellite IDs and not matrix indices - are they the same?
# what is ISL shift?
# Think this is the ns3 file that calls on ISL topology file built by this program - hypatia/ns3-sat-sim/simulator/contrib/satellite-network/model
# /topology-satellite-network.cc - this is function called by something else - need to edit so reads ISL file with relevant counter (mybe pass argument to function)
# above only contains functions, so need to find where the ns3 namespace is called (looking for function TopologySatelliteNetwork::ReadISLs()) - reading from
# file in directory m_satellite_network_dir/isls.txt (are there multiple directories?). Files are created with ptr val basicSimulation
# (passed to TopologySatelliteNetwork::TopologySatelliteNetwork constructor)
# hypatia/ns3-sat-sim/simulator/scratch/main_satnet
# /main_satnet.cc calls the <TopologySatelliteNetwork> object
# Could be key for altitude and orbit decay - https://github.com/snkas/hypatia/blob/master/paper/satellite_networks_state/main_starlink_550.py
# Need to use paper code - they have all the networks you want to test/work with (Starlink, Kuiper, Telesat, etc) - READ THROUGH
# LIST OF INBUILT RESTRAINTS - e.g altitude, cone radius, etc.
# Some satellites will always be visible to one another
# ADAPT SO PER SNAPSHOT


# Libraries
import generate_tles_from_scratch as hypatia_data
import read_tles as hypatia_read_data
import os
from skyfield.api import wgs84  # https://rhodesmill.org/skyfield/toc.html - recommended by astropy for calculating information about
# satellites described by TLE format (will need to confirm everything is accurate and maybe write my own function
# versions, but this will allow me to create a working version of the code)
# Used to calculate TEME Coordinates from TLE
from skyfield.api import EarthSatellite, load
import math
from astropy.time import Time
from astropy import units as u
import random

# Function calls on Hypatia software function generate_tles_from_scratch_with_sgp
def data_generation():
    # Generates Sample File
    # Test data taken from https://github.com/snkas/hypatia/blob/master/satgenpy/tests/test_tles.py
    # This generates the physical layout of the Starlink-550 satellite network
    if os.path.isfile("./Starlink_Test_TLEs.txt.tmp") is False:
        hypatia_data.generate_tles_from_scratch_with_sgp("Starlink_Test_TLEs.txt.tmp", "Starlink-550", 72, 22, True, 53,
                                                         0.0000001, 0.0, 13.66)
    return


# https://docs.astropy.org/en/latest/coordinates/satellites.html
def maximum_communication_distance(data_file):

    # Get sample satellite to calculate maximum communication distance
    with open(data_file, 'r') as f:
        _ = f.readline()
        name = f.readline()
        s = f.readline().strip()
        t = f.readline().strip()

    # Convert TLES to Geocentric Coordinates
    ts = load.timescale()
    sample_satellite = EarthSatellite(s, t, name, ts)

    # Assuming satellite is travelling in a circular orbit, assume orbital height is same at all times
    t = ts.now()

    # Convert to geocentric coordinate system
    geocentric_coords = sample_satellite.at(t)

    height_above_earth = float(wgs84.height_of(geocentric_coords).km)

    # From NSGA-III Paper, use formula to calculate maximum communication distance between two satellites at that height
    # Value for Earth Radius = 6378137.0 m - https://en.wikipedia.org/wiki/Earth_radius
    earth_radius = 6378.1370
    max_distance = 2 * math.sqrt(pow(earth_radius + height_above_earth, 2) - pow(earth_radius, 2))

    return max_distance


# Function calls on Hypatia Software function read_tles
def format_tle_data(filename):
    return hypatia_read_data.read_tles(filename)


# Hope to improve this by looking at NSGA-III paper and the antenna direction - will need a way to get antenna direction
def satellite_visibility_function(distance_matrix, max_dist):
    total_satellites = len(distance_matrix)
    visibility_matrix = [[0 for _ in range(total_satellites)] for _ in range(total_satellites)]
    for i in range(total_satellites):
        for j in range(i+1, total_satellites):
            # Accept max distance - although considering rejecting if >= max_dist
            if distance_matrix[i][j] <= max_dist:
                visibility_matrix[i][j] = 1
                visibility_matrix[j][i] = 1
    return visibility_matrix

def time_visibility_function(visibility_matrices, snapshot_num):
    total_satellites = len(visibility_matrices[0][0])
    tv_matrix = [[[0 for _ in range(total_satellites)] for _ in range(total_satellites)] for _ in range(snapshot_num)]
    for t in range(snapshot_num):
        for i in range(total_satellites):
            for j in range(i+1, total_satellites):
                future_snapshot = t + 1
                for k in range(snapshot_num):
                    if future_snapshot == snapshot_num:
                        future_snapshot = 0
                    if visibility_matrices[future_snapshot][i][j] == 1:
                        tv_matrix[t][i][j] += 1
                        tv_matrix[t][j][i] += 1
                    else:
                        break
    return tv_matrix

def cost_function(visibility, time_visibility, distance, alpha, beta):
    total_satellites = len(visibility)
    cost_matrix = [[0 for _ in range(total_satellites)] for _ in range(total_satellites)]
    for i in range(total_satellites):
        for j in range(total_satellites):
            if visibility[i][j] == 0:
                cost_matrix[i][j] = -1
            else:
                cost_matrix[i][j] = (alpha * time_visibility) + (beta * distance)

    return cost_matrix


# This file is used to store raw TLE data rather than convert it to pyephem format
# Please note that this function (the following 35 lines of code) is adapted from Hypatia code (read_tles.py in the
# satgenpy module)
def read_file(file_name):

    tles_data = []
    with open(file_name, 'r') as f:
        n_orbits, n_sats_per_orbit = [int(n) for n in f.readline().split()]
        universal_epoch = None
        i = 0
        for tles_line_1 in f:
            tles_line_2 = f.readline()
            tles_line_3 = f.readline()

            # Retrieve name and identifier
            name = tles_line_1
            sid = int(name.split()[1])
            if sid != i:
                raise ValueError("Satellite identifier is not increasing by one each line")
            i += 1

            # Fetch and check the epoch from the TLES data
            # In the TLE, the epoch is given with a Julian data of yyddd.fraction
            # ddd is actually one-based, meaning e.g. 18001 is 1st of January, or 2018-01-01 00:00.
            # As such, to convert it to Astropy Time, we add (ddd - 1) days to it
            # See also: https://www.celestrak.com/columns/v04n03/#FAQ04
            epoch_year = tles_line_2[18:20]
            epoch_day = float(tles_line_2[20:32])
            epoch = Time("20" + epoch_year + "-01-01 00:00:00", scale="tdb") + (epoch_day - 1) * u.day
            if universal_epoch is None:
                universal_epoch = epoch
            if epoch != universal_epoch:
                raise ValueError("The epoch of all TLES must be the same")

            # Finally, store the satellite information
            # tles_data.append(ephem.readtle(tles_line_1, tles_line_2, tles_line_3))
            tles_data.append([tles_line_1.strip(), tles_line_2.strip(), tles_line_3.strip()])

    return tles_data

# Calculates the distance between satellites in the network
# Need to check distance calculation is correct (research wsg84)
def dist(satellite_i, satellite_j, time_stamp):

    # NEED TO CHECK POINT IN ORBIT BASED ON TIME IS CORRECT
    # Assuming 60 mins in hour, 60 seconds in minute, etc
    hours = 0
    minutes = 0

    # Convert time to
    while time_stamp - 3600 > 0:
        hours += 1
        time_stamp -= 3600
    while time_stamp - 60 > 0:
        minutes += 1
        time_stamp -= 60
    seconds = time_stamp

    # Looked at dates included with Starlink data and based on data (need to check time format is correct and is using jpl)
    # print(satellite_i.epoch.utc)

    # Set time
    ts = load.timescale()
    t = ts.utc(2000, 1, 1, hours, minutes, seconds)

    # Convert satellites i and j to geocentric coordinates - get measurement in km
    satellite_i = EarthSatellite(satellite_i[1], satellite_i[2], satellite_i[0], ts)
    satellite_j = EarthSatellite(satellite_j[1], satellite_j[2], satellite_j[0], ts)

    difference = (satellite_i - satellite_j).at(t).position.km

    distance = math.sqrt(sum([pow(difference[i], 2) for i in range(3)]))

    return distance


# Calculates the distance matrix for a snapshot of the network
def distance_function(satellites, total_satellites, snapshot_time):

    dist_matrix = [[0.0 for _ in range(total_satellites)] for _ in range(total_satellites)]

    # Alter satellite description to describe its position at snapshot_time t

    for i in range(total_satellites):
        for j in range(i + 1, total_satellites):
            distance = dist(satellites[i], satellites[j], snapshot_time)
            dist_matrix[i][j] = distance
            dist_matrix[j][i] = distance

    return dist_matrix

# Returns degree constrained minimum spanning tree of network - uses primal-cut branch algorithm (see paper for references)
# Greedy Algorithm
def degree_constrained_mst(cost_matrix, constraints):
    total_satellites = len(cost_matrix)
    tree = [[0 for _ in range(total_satellites)] for _ in range(total_satellites)]

    tree_vertices = []

    degrees = [0 for _ in range(total_satellites)]

    # Select arbitrary vertex to put into tree vertex
    tree_vertices.append(random.randint(0, total_satellites))

    minimum = -1
    # Stores closest vertex not in tree
    closest_vertex = 0
    connected_to = 0

    while len(tree_vertices) < total_satellites:

        for k in range(total_satellites):
            if (k not in tree_vertices) or (degrees[k] == constraints[k]):
                continue
            else:
                min_row = min(cost_matrix[k])
                closest_vertex_temp = cost_matrix[k].index(min_row)
                if ((minimum < 0) or (min_row < minimum)) and (closest_vertex_temp not in tree_vertices) and (degrees[closest_vertex_temp] != constraints[closest_vertex_temp]):
                    minimum = min_row
                    closest_vertex = closest_vertex_temp
                    connected_to = k


        tree_vertices.append(closest_vertex)

        tree[connected_to][closest_vertex] = 1
        tree[closest_vertex][connected_to] = 1

        degrees[closest_vertex] += 1

    return tree

degree_constrained_mst([[0, 56, 91, 350], [12, 0, 45, 61], [78, 45, 0, 31], [3, 2, 1, 0]], [3, 3, 3, 3])









# Add ISLs to increase network connectivity
def increase_connectivity(network, constraints, current_connection_number, costs):

    available_isls = []

    # Ignore all ISLs already established
    for i in range(len(network)):
        for j in range(i+1, len(network)):
            if network[i][j] == 1:
                costs[i][j] = -1
            else:
                if costs[i][j] != -1:
                    available_isls.append([costs[i][j], [i, j]])

    # Sort list of available ISLs in order of increasing cost
    available_isls.sort()

    for k in available_isls:
        satellites = k[1]
        if (current_connection_number[satellites[0]] == constraints[satellites[0]]) or (current_connection_number[satellites[1]] == constraints[satellites[1]]):
            continue
        else:
            network[satellites[0]][satellites[1]] = 1
            network[satellites[1]][satellites[0]] = 1

    return network

# Writes results to file (in correct format) for hypatia to process
def write_results_to_file(output_file, results, snapshot_id):
    # Write results to file for snapshot k
    with open(output_file + "_snapshot_" + str(snapshot_id) + ".txt", 'w+') as f:
        total_satellites = len(results)
        for i in range(total_satellites):
            for j in range(i+1, total_satellites):
                if results[i][j] == 1:
                    f.write(str(i) + ' ' + str(j) + '\n')

# def heuristic_topology_design_algorithm_isls(output_filename_isls, n_orbits, n_sats_per_orbit, orbital_period, num_snapshot, isl_shift, idx_offset=0):
def heuristic_topology_design_algorithm_isls(satellites, n_orbits, n_sats_per_orbit, orbit_period, num_snapshot, max_comm_dist, degree_constraints, output_filename_isls):

    # Check satellite network has a sufficient number of satellites and orbits
    if n_orbits < 3 or n_sats_per_orbit < 3:
        raise ValueError("Number of x and y must each be at least 3")

    # Calculate the time between each snapshot
    snapshot_interval = orbit_period/num_snapshot

    # Calculate the total number of satellites in the network
    total_satellites = n_orbits * n_sats_per_orbit

    # Get TLEs-formatted data
    tles_data = read_file("Starlink_Test_TLEs.txt.tmp")

    # Calculate distance matrices
    # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE TIME VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
    distance_matrices = [distance_function(tles_data, total_satellites, snapshot_interval*k) for k in range(num_snapshot)]

    # Calculate visibility and time visibility matrices for all snapshots

    # Calculate visibility matrix
    # [i][j] is set to 0 if satellites are not visible to one another
    visibility_matrices = [satellite_visibility_function(distance_matrices[k], max_comm_dist) for k in range(num_snapshot)]

    # Calculate time visibility matrix for each snapshot
    time_visibility_matrices = time_visibility_function(visibility_matrices, num_snapshot)

    for snapshot in range(0, num_snapshot):    # NOTE TO SELF - DON'T FORGET TIME!!!

        # Initialise list of the current number of active ISLs each satellite has
        current_isl_number = [0 for _ in range(len(satellites))]

        # Calculate cost matrix
        # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE TIME VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
        # Include bandwidth matrix here - look for more metrics

        # Setting hyperparameters to 1 initially
        alpha, beta = 1, 1

        cost_matrix = cost_function(visibility_matrices[snapshot], time_visibility_matrices[snapshot],
                                    distance_matrices[snapshot], alpha, beta)

        # Calculate degree constrained minimum spanning tree
        # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE TIME VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
        # Original Paper on DCMST and the Primal Algorithm:
        # https://www-sciencedirect-com.ezphost.dur.ac.uk/science/article/pii/0305054880900222?via%3Dihub (NEED TO CITE
        # IN PAPER). Also, consulted https://en.wikipedia.org/wiki/Prim%27s_algorithm, as the paper references Prim's
        # algorithm. NOTE TO SELF - WALK THROUGH PAPER - CHECK CORRECTNESS OF ADAPTATION

        tree = degree_constrained_mst(cost_matrix, degree_constraints)

        # Adding Edges
        # Cost matrix has -1 where link is not possible
        isls = increase_connectivity(tree, degree_constraints, current_isl_number, cost_matrix)

        # Convert list_isls to correct format and save results of algorithm in file
        # Note to self - may want to throw error if list_isls is None (and catch)
        write_results_to_file(output_filename_isls, isls, snapshot)


# Testing Code
if __name__ == '__main__':

    # Create test data
    data_generation()

    # Read test data into a variable
    data = format_tle_data("Starlink_Test_TLEs.txt.tmp")

    # Read description of satellites and their unique orbits
    satellite_data = data["satellites"]

    # Take the first satellite and calculate (using mean motion) the orbital period of a satellite
    # Used physics formula - https://en.wikipedia.org/wiki/Mean_motion. Mean motion given in revolutions per day.
    # Convert to seconds following astropy's assumption that 1 day ==86400 seconds - https://docs.astropy.org/en/stable/time/index.html
    orbital_period = 86400 / (satellite_data[0]).n

    # Find the maximum communication distance between two satellites (may vary as satellite altitudes vary)
    max_communication_dist = maximum_communication_distance("Starlink_Test_TLEs.txt.tmp")

    # This is the maximum distance a satellite can establish signal (transmission power) - need to research for Kuiper
    # and StarLink satellites
    max_transmission_dist = 1000000

    # Find smaller of these two numbers - the maximum distance before Earth is in the way or the max transmission
    # distance (due to satellite power constraints)
    max_communication_dist = min(max_communication_dist, max_transmission_dist)

    # Initialise degree constraint for each satellite
    satellite_degree_constraints = [3 for _ in range(len(satellite_data))]

    # Initially, let us say that there are 1000 snapshots over the orbital period - will eventually change this so it
    # the time between snapshots is the time when network has no visibility changes

    # Run topology generation algorithm
    heuristic_topology_design_algorithm_isls(satellite_data, data["n_orbits"], data["n_sats_per_orbit"], orbital_period, 1000, max_communication_dist, satellite_degree_constraints, "isls")



# NOTE TO SELF - astropy may calculate orbital period for you
# https://docs.astropy.org/en/latest/coordinates/satellites.html

# References:
# Hypatia - https://github.com/snkas/hypatia/tree/master
# Pyephem - https://github.com/brandon-rhodes/pyephem
# Pyephem Documentation - https://rhodesmill.org/pyephem/quick
# SkyField Documentation - https://rhodesmill.org/skyfield/
# Prim's Algorithm - https://en.wikipedia.org/wiki/Prim%27s_algorithm
# DCMST Primal Algorithm - https://www-sciencedirect-com.ezphost.dur.ac.uk/science/article/pii/0305054880900222?via%3Dihub
# TLE Definitions - https://platform-cdn.leolabs.space/static/files/tle_definition.pdf?7ba94f05897b4ae630a3c5b65be7396c642d9c72
# Orbital Distance - https://space.stackexchange.com/questions/27872/how-to-calculate-the-orbital-distance-between-2-satellites-given-the-tles
# Earth Radius - https://en.wikipedia.org/wiki/Earth_radius
# World Geodetic System - https://en.wikipedia.org/wiki/World_Geodetic_System#Definition

# Notes
# A description of each satellite is stored in an ephem EarthSatellite object with the following
# parameters (see pyephem documentation):
# From Pyephem documentation:
# epoch — Reference epoch
# n — Mean motion, in revolutions per day
# inc — Inclination (°)
# raan — Right Ascension of ascending node (°)
# e — Eccentricity
# ap — Argument of perigee at epoch (°)
# M — Mean anomaly from perigee at epoch (°)
# decay — Orbit decay rate in revolutions per day, per day
# drag — Object drag coefficient in per earth radii
# orbit — Integer orbit number of epoch
