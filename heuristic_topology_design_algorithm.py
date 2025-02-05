#TODO
# When running tests, check all data_generation params are correct (including defaults, e.g. does perigee arg change for
# Kuiper?)
# import ephem or pyephem
# Find maximum transmission dist for Starlink, Kuiper and telesat - 27000 paper sets at 5014 km - also mentioned in code
# Check correctness of distance measure between satellites (changing so calculated position of each and then found dist
# between them - need to check!!)
# Improve visibility by looking at NSGA-III paper and the antenna direction - will need a way to get antenna direction
# Visibility - if less than max dist, are they always visible - check papers!!!!!!!!!!!
# SAVE DISTANCE MATRICES AND VISIBILITY MATRICES TO FILES SO DON'T HAVE TO CALCULATE AGAIN

# Import Relevant Libraries
from astropy.time import Time
from astropy import units as u
import ephem
import generate_tles_from_scratch as hypatia_data
import math
import numpy as np
import os
import random
import read_tles as hypatia_read_data
from skyfield.api import EarthSatellite, wgs84, load  # recommended by astropy for calculating information about
# satellites described by TLE (convert to TEME), timescale is only defined once
# (https://rhodesmill.org/skyfield/api-time.html#skyfield.timelib.Time)

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()


# Function calls on Hypatia software function generate_tles_from_scratch_with_sgp to build physical description of
# satellite network using TLE coordinate system and stores description in temporary file. Orbit eccentricity must be > 0
# (according to Hypatia) - set close to 0 to ensure orbit approx. circular
def data_generation(file_name, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree,
                    mean_motion_rev_per_day, phase_diff=True, eccentricity=0.0000001, arg_of_perigee_degree=0.0):
    # Check file does not already exist
    if os.path.isfile("./" + file_name) is False:
        hypatia_data.generate_tles_from_scratch_with_sgp(file_name, constellation_name, num_orbits, num_sats_per_orbit,
                                                         phase_diff, inclination_degree, eccentricity,
                                                         arg_of_perigee_degree, mean_motion_rev_per_day)
    return

# Function calls on Hypatia software function read_tles to build dictionary containing the total orbits in the
# constellation, as well as the number of satellites per orbit, the epoch of the network, and a description of each
# satellite's position
def format_tle_data(file_name):
    return hypatia_read_data.read_tles(file_name)


# Calculates (using mean motion) the orbital period of a satellite, using physics formula -
# https://en.wikipedia.org/wiki/Mean_motion. Mean motion given in revolutions per day and converted to seconds following
# astropy's assumption that 1 day == 86400 seconds - https://docs.astropy.org/en/stable/time/index.html
def orbital_period_calculation(satellite_description, num_sat):
    if isinstance(satellite_description, ephem.EarthSatellite):
        return 86400 / satellite_description.n
    else:
        # Return longest orbital period (i.e. for the satellites that orbit at the greatest altitude
        return max([86400/satellite_description[x].n for x in range(num_sat)])


# Calculates a satellite's altitude above Earth given TLES coordinates of satellite
def satellite_height_above_earth(name, s, t):

    # Convert TLES to Geocentric Coordinates
    sample_satellite = EarthSatellite(s, t, name, ts)

    # Assuming satellite is travelling in a circular orbit (eccentricity is close to 0), and assume orbital height
    # is approximately the same at all times. Set time to be fixed, e.g. 1st January 2000 Midnight (ensures
    # deterministic behaviour when calculating position, as opposed to using ts.now()). Find satellite position at given
    # time (as circular orbit, should not matter at what time the satellite position is recorded and retrieve height
    # above earth value in km

    height_above_earth = float(wgs84.height_of(sample_satellite.at(ts.tdb(2000, 1, 1, 0, 0))).km)

    return height_above_earth

# Based on recommendations in https://docs.astropy.org/en/latest/coordinates/satellites.html
# From NSGA-III Paper, use formula to calculate maximum communication distance between two satellites at given altitude
def maximum_communication_distance(data_file, num_sat):

    # Constant - source of value: https://en.wikipedia.org/wiki/Earth_radius and https://github.com/AtlantaEmrys2002/hypatia/blob/master/paper/satellite_networks_state/main_starlink_550.py
    earth_radius = 6378.135 # 6378.1370

    # Get sample satellite TLES coordinates to calculate maximum communication distance
    with open(data_file, 'r') as f:
        lines = [line.strip() for line in f]

    # Get data needed for calculating satellite position and calculate each satellite's altitude. Find the lowest
    # orbital altitude - this will provide the "smallest" maximum communication distance, as closer to Earth.
    lowest_satellite_altitude = min([satellite_height_above_earth(lines[line_index], lines[line_index + 1],
                                                                  lines[line_index + 2]) for line_index in
                                                                  range(1, num_sat, 3)])

    # Return the maximum communication distance between two satellites at the lowest satellite altitude in the network
    return 2 * math.sqrt(pow(earth_radius + lowest_satellite_altitude, 2) - pow(earth_radius, 2))

# Returns the maximum transmission distance for satellite in network (values found through research). TEMPORARILY SET TO LARGE VALUES!!!!!!!
def maximum_transmission_distance(name):
    if 'Starlink' in name:
        return 10000
    elif 'Telesat' in name:
        return 10000
    else:  # Kuiper
        return 10000


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
            tles_data.append([tles_line_1.strip(), tles_line_2.strip(), tles_line_3.strip()])

    return tles_data

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


# Calculates the distance matrix for a snapshot of the network
# def distance_function(satellites, total_satellites, snapshot_time):
def distance_function(satellites, total_satellites):

    # Calculate distance (Euclidean) between all satellite pairs i and j in the network and alter satellite description
    # to describe its position at snapshot_time t
    dist_matrix = [[np.linalg.norm(satellites[i] - satellites[j]) for j in range(i + 1, total_satellites)]
                   for i in range(total_satellites)]

    # Distance between satellite and itself is 0km
    for i in dist_matrix:
        i.insert(0, 0)

    # Add in previously calculated values (to create symmetric matrix)
    for i in range(total_satellites):
        dist_matrix[i] = [dist_matrix[j][i] for j in range(i)] + dist_matrix[i]

    return dist_matrix


# Calculate whether each satellite pair i, j can establish a connection (i.e. if they are visible to one another)
def visibility_function(distance_matrix, max_dist, total_satellites):

    # Calculates whether satellites are within visible distance to each other
    visibility_matrix = np.asarray(distance_matrix) <= max_dist
    visibility_matrix = np.where(visibility_matrix == True, visibility_matrix, 0)

    # Satellites cannot be visible to themselves
    for i in range(total_satellites):
        visibility_matrix[i][i] = 0

    return visibility_matrix



def heuristic_topology_design_algorithm_isls(input_file_name, satellites, total_satellites, orbit_period, max_comm_dist, degree_constraints, output_filename_isls):

    # Check satellite network has a sufficient number of satellites and orbits
    if total_satellites < 3:
        raise ValueError("Number of satellites must be greater than 3.")

    # In original Hypatia paper, snapshots of 100ms were utilised - this is continued here (all times are in seconds)
    snapshot_interval = 0.001

    # The number of snapshots over an orbital period
    num_snapshot = int(orbit_period/snapshot_interval)

    # TEMPORARY CHANGE - TAKES A LONG TIME TO TAKE SNAPSHOTS
    num_snapshot = 5

    # Get TLEs-formatted data
    tles_data = read_file(input_file_name)

    # Convert TLES data to Skyfield EarthSatellite objects - used to convert satellite positions to geocentric
    # coordinates - all measurements in km
    earth_satellite_objects = [EarthSatellite(satellite_i[1], satellite_i[2], satellite_i[0], ts) for satellite_i in tles_data]

    ### DISTANCE MATRIX ###

    # Calculate distance matrices - set time t (using TDB time format) to indicate point in orbital period when snapshot
    # taken

    # Calculate the time (in TDB format) at which each snapshot is taken
    snapshot_times = [snapshot_time_stamp(snapshot_interval * k) for k in range(num_snapshot)]

    # Calculate distance matrices for each snapshot
    distance_matrices = [distance_function([i.at(k).position.km for i in earth_satellite_objects], total_satellites) for k in
                         snapshot_times]

    ### VISIBILITY MATRICES ###

    # Calculate visibility and time visibility matrices for all snapshots

    # Calculate visibility matrix - [i][j] is set to 0 if satellites are not visible to one another
    visibility_matrices = [visibility_function(distance_matrices[k], max_comm_dist, total_satellites) for k in range(num_snapshot)]

    # GOT TO HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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



# Main Function used to test code - constellation name specified name of network to build topology for
# Multishell indicates if satellites orbit at different heights
def main(file_name, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree, mean_motion_rev_per_day, multi_shell=False):

    # Calculate the number of satellites in the network
    total_sat = num_sats_per_orbit * num_orbits

    # Generate test data using network description from https://github.com/snkas/hypatia/blob/master/satgenpy/tests/test
    # _tles.py
    data_generation(file_name, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree,
                    mean_motion_rev_per_day)

    # Read test data into appropriate data structure (dictionary)
    data = format_tle_data(file_name)

    # Extract description of satellite positions and unique orbits from data
    satellite_data = data["satellites"]

    # Calculate orbital period of network (or maximum orbital period if satellites orbit at different altitudes)
    if multi_shell is False:
        orbital_period = orbital_period_calculation(satellite_data[0], total_sat)
    else:
        orbital_period = orbital_period_calculation(satellite_data, total_sat)

    # Find the maximum communication distance between two satellites (may vary as satellite altitudes vary)
    max_communication_dist = maximum_communication_distance(file_name, total_sat)

    # This is the maximum distance a satellite can establish signal (transmission power) - need to research for Kuiper
    # and StarLink satellites
    max_transmission_dist = maximum_transmission_distance(constellation_name)

    # Find smaller of these two numbers - the maximum distance before Earth is in the way or the max transmission
    # distance (due to satellite power constraints)
    max_communication_dist = min(max_communication_dist, max_transmission_dist)

    # Initialise degree constraint for each satellite - can be changed based on technical specifications of satellites
    satellite_degree_constraints = [3 for _ in range(len(satellite_data))]

    # Initially, let us say that there are 1000 snapshots over the orbital period - will eventually change this so it
    # the time between snapshots is the time when network has no visibility changes.

    # GOT TO HERE - GOT TO CHANGE SO THAT IT ONLY RETURNS TOPOLOGY FOR CURRENT SNAPSHOT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Run topology generation algorithm
    heuristic_topology_design_algorithm_isls(file_name, satellite_data, total_sat, orbital_period, max_communication_dist, satellite_degree_constraints, "isls")


# Used for testing
# Data from: https://github.com/AtlantaEmrys2002/hypatia/tree/master/paper/satellite_networks_state
main("constellation_tles.txt.tmp", "Starlink-550", 72, 22, 53, 15.19)


# References
# DCMST Primal Algorithm - https://www-sciencedirect-com.ezphost.dur.ac.uk/science/article/pii/0305054880900222?via%3Dihub
# Earth Radius - https://en.wikipedia.org/wiki/Earth_radius
# Hypatia - https://github.com/snkas/hypatia/tree/master
# Orbital Distance - https://space.stackexchange.com/questions/27872/how-to-calculate-the-orbital-distance-between-2-satellites-given-the-tles
# Prim's Algorithm - https://en.wikipedia.org/wiki/Prim%27s_algorithm
# Pyephem Code - https://github.com/brandon-rhodes/pyephem
# Pyephem Documentation - https://rhodesmill.org/pyephem/quick
# SkyField Documentation - https://rhodesmill.org/skyfield/ & https://rhodesmill.org/skyfield/toc.html
# TLE Definitions - https://platform-cdn.leolabs.space/static/files/tle_definition.pdf?7ba94f05897b4ae630a3c5b65be7396c642d9c72
# World Geodetic System - https://en.wikipedia.org/wiki/World_Geodetic_System#Definition