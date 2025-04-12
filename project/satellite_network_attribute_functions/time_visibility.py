# Libraries
import numpy as np


def time_visibility_function(snapshot_num: int, total_satellites: int, initial_id_num: int, constellation_name: str):
    """
    For a given snapshot, calculates how long each pair of satellites within the network will remain visible to one
    another (in terms of number of future snapshots).

    :param snapshot_num: the number of satellites within the network
    :param total_satellites:
    :param initial_id_num:
    :param constellation_name:
    :return:
    """

    # Declare common file name of all visibility matrices
    file_name = "./" + constellation_name + "/visibility_matrices/visibility_matrix_"

    # If changed[i][j] == 1, indicates that satellite visibility has not changed
    changed = np.ones((total_satellites, total_satellites))

    # Initialise time visibility matrix for current snapshot
    tv = np.zeros((total_satellites, total_satellites))

    # Looks at next snapshot after current one being analysed (i.e. the snapshot with id == id_num)
    current_id = initial_id_num + 1

    # If the next snapshot is equal to the number of snapshots in the simulation. out of range and go back to start of
    # orbital period (i.e. the snapshot with id == 0
    if current_id == snapshot_num:
        current_id = 0

    # Load relevant visibility matrices
    previous = np.load(file_name + str(initial_id_num) + ".npy")
    current = np.load(file_name + str(current_id) + ".npy")

    # For each snapshot, calculate what has changed visibility-wise between last snapshot and current snapshot (i.e.
    # what has become invisible)
    for t in range(1, snapshot_num):

        # Find all satellite pairs that are visible to one another in both snapshots
        changed_step = np.logical_and(previous, current)

        # Make sure hasn't changed in the past and now visible again - only important that visible satellites do not
        # change
        changed = np.logical_and(changed_step, changed)

        # Add 1 to every satellite time visibility where visibility of satellite has not changed
        tv[changed] += 1

        # Calculate next relevant visibility matrix
        if current_id == snapshot_num - 1:
            current_id = 0
        else:
            current_id += 1

        # Load next visibility matrices
        previous = current
        current = np.load(file_name + str(current_id) + ".npy")

    return tv.astype(int)
