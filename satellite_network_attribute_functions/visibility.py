# Libraries
import numpy as np


def visibility_function(distance_matrix, max_dist) -> np.ndarray:
    """
    For each pair of satellites within the network, calculates whether two satellites are visible to one another based
    on current positions.

    :param distance_matrix:
    :param max_dist:
    :return: visibility_matrix

    """
    # Calculates whether satellites are within visible distance of each other by seeing how close they are - max
    # distance depends on the height at which the satellites orbit
    visibility_matrix = np.asarray(distance_matrix) <= max_dist

    # Satellites cannot be visible to themselves - therefore, set diagonal of matrix to zero
    np.fill_diagonal(visibility_matrix, False)

    return visibility_matrix
