# Libraries
import numpy as np


def sunlight_function(satellites_in_sun: list, total_satellites: int) -> np.ndarray:
    """
    Calculates whether each satellite is within sunlight or not at given time (i.e. naive approximation of how
    vulnerable to solar flares). All edges that incident to at least one node (satellite) that is in sunlight is
    assigned 1 to indicate satellite is vulnerable to solar flares. This approximation may be improved in future work.

    :param satellites_in_sun: matrix indicating whether satellite in path of sunlight
    :param total_satellites: the number of satellites within the network
    :return: matrix where 1 indicates that satellite or satellite it could connect to is vulnerable to solar flares.
    """
    sunlight_matrix = np.zeros((total_satellites, total_satellites))

    in_sun = np.argwhere(satellites_in_sun).T[0]

    sunlight_matrix[np.ix_(in_sun, in_sun)] = 1

    return sunlight_matrix
