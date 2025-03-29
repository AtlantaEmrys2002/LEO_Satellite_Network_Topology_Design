# Libraries
import numpy as np
import metrics


def read_isl_file(isl_file_name: str, num_satellites: int) -> np.ndarray:
    """
    Reads description of topology from file for analysis according to multiple metrics. Often used by cost function
    optimisation code.

    :param isl_file_name: name of .npy file that contains ISl pairs describing a topology
    :param num_satellites: the number of satellites within the network
    :return:
    """
    # Read in topology built for given snapshot
    topology_isls = np.loadtxt(isl_file_name).astype(int).T

    # Create topology matrix for ISLs
    topology_matrix = np.zeros((num_satellites, num_satellites))
    topology_matrix[topology_isls[0], topology_isls[1]] = 1

    return topology_matrix


def measure_static(constellation_name: str, topology_file_name: str, num_satellites: int) -> (
        tuple)[float, float, float, int]:
    """
    Calculates the maximum propagation delay, mean propagation delay, and average hop count (as well as returning link
    churn for completeness) for static satellite network topology.
    :param constellation_name: name of satellite network constellation, e.g. Starlink-550
    :param topology_file_name: name of .npy file that contains ISl pairs describing a topology
    :param num_satellites: the number of satellites within the network
    :return:
    """
    # Read in topology
    topology = read_isl_file(topology_file_name, num_satellites)

    # Calculates values using distance matrix for snapshot 0
    distance = np.load("./" + constellation_name + "/distance_matrices/dist_matrix_0" + ".npy")

    # Calculate propagation delay (max and mean) for static topology
    max_pd, mean_pd = metrics.propagation_delay(topology, distance, topology[0].size)

    # Calculate hop count
    av_hop_count = metrics.hop_count(topology, distance, topology[0].size)

    # No need to calculate link churn
    link_churn = 0

    return max_pd, mean_pd, av_hop_count, link_churn


def measure_dynamic(constellation_name: str, topology_file_location: str, num_satellites: int, num_snapshots: int) -> (
        tuple)[float, float, float, int]:
    """
    Calculates the maximum propagation delay, mean propagation delay, average hop count, and link churn for dynamic
    satellite network topology.
    :param constellation_name: name of satellite network constellation, e.g. Starlink-550
    :param topology_file_location: location of .npy file in which file(s) containing topology description(s) stored
    :param num_satellites: the number of satellites within the network
    :param num_snapshots: the number of snapshots of the network over one orbit for which a topology is constructed
    """
    # Initialise values
    av_max_pd, av_mean_pd, av_av_hop_count, link_churn = 0, 0, 0, 0

    # For each topology over the course of one orbital period
    for k in range(num_snapshots):

        # Read in current topology
        topology = read_isl_file(topology_file_location + "/isls_" + str(k) + ".txt", num_satellites)

        # Read in distance matrix for current snapshot
        distance = np.load("./" + constellation_name + "/distance_matrices/dist_matrix_" + str(k) + ".npy")

        # Calculate propagation delay (max and mean) for static topology
        max_pd, mean_pd = metrics.propagation_delay(topology, distance, topology[0].size)

        av_max_pd += max_pd
        av_mean_pd += mean_pd

        # Calculate hop count
        av_hop_count = metrics.hop_count(topology, distance, topology[0].size)

        av_av_hop_count += av_hop_count

    # Calculate link churn
    link_churn = metrics.link_churn(topology_file_location, num_snapshots, num_satellites)

    # Find average of delays and hop count across all snapshots
    av_max_pd /= num_snapshots
    av_mean_pd /= num_snapshots
    av_av_hop_count /= num_snapshots

    if av_max_pd == 0 or av_mean_pd == 0 or av_av_hop_count == 0:
        raise ValueError("Check metric calculations - it is not possible for satellites in two different positions to "
                         "have a propagation delay value equal to 0. Similarly, hop count is always >= 1.")

    return av_max_pd, av_mean_pd, av_av_hop_count, link_churn
