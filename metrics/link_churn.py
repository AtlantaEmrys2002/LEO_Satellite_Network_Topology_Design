# Libraries
import numpy as np


def link_churn(location, snapshot_num, num_satellites):
    """
    Calculates link churn for satellite network over one orbital period (as determined by number of snapshots of network
     taken).
    :param location:
    :param snapshot_num:
    :param num_satellites:
    :return:
    """
    # LOAD INITIAL TOPOLOGY MATRICES #

    # Read in topology built for given snapshot
    previous_isls = np.loadtxt(location + "/isls_0.txt").astype(int).T

    # Create topology matrix for previous
    previous = np.zeros((num_satellites, num_satellites))
    previous[previous_isls[0], previous_isls[1]] = 1

    # Read in topology built for given snapshot
    current_isls = np.loadtxt(location + "/isls_1.txt").astype(int).T

    # Create topology matrix for current
    current = np.zeros((num_satellites, num_satellites))
    current[current_isls[0], current_isls[1]] = 1

    # Keeps track of current topology comparing to previous
    current_id = 1

    # Keeps track of total number of link changes over the course of one orbital period
    total = 0

    # CALCULATE NUMBER OF LINK CHANGES BETWEEN EACH PAIRWISE SNAPSHOTS #

    # For each snapshot, find the changed links (i.e. ISLs that have connected between previous snapshot and current
    # snapshot
    for t in range(1, snapshot_num):

        # Find changed links
        total += np.sum(np.logical_xor(previous, current))

        # Calculate next relevant topology
        if current_id == snapshot_num - 1:
            current_id = 0
        else:
            current_id += 1

        # Load next topology matrices

        # Create topology matrix for previous
        previous = current

        # Read in topology built for given snapshot
        current_isls = np.loadtxt(location + "/isls_" + str(current_id) + ".txt").astype(int).T

        # Create topology matrix for current
        current = np.zeros((num_satellites, num_satellites))
        current[current_isls[0], current_isls[1]] = 1

    return total
