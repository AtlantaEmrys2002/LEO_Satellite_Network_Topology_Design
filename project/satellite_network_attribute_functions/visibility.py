# Libraries
import numpy as np


def visibility_function(distance_matrix: np.ndarray, max_dist: float) -> np.ndarray:
    """
    For each pair of satellites within the network, calculates whether two satellites are visible to one another based
    on current positions.

    :param distance_matrix: an adjacency matrix where each element represents the distance between a satellite pair.
    :param max_dist: the maximum distance between satellites over which an ISL may be established.
    :return: adjacency matrix such that if satellites i and j are visible to one another, visibility[i][j] == 1, else 0

    """
    # Calculates whether satellites are within visible distance of each other by seeing how close they are - max
    # distance depends on the height at which the satellites orbit
    visibility_matrix = np.asarray(distance_matrix) <= max_dist

    # Satellites cannot be visible to themselves - therefore, set diagonal of matrix to zero
    np.fill_diagonal(visibility_matrix, False)

    return visibility_matrix
