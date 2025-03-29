# Libraries
from .read_tles import read_tles


def format_tle_data(file_name: str) -> dict:
    """
    Function calls on Hypatia read_tles to build dictionary containing the total number of orbits in the
    constellation, as well as the number of satellites per orbit, the epoch of the network, and a description of each
    satellite's position.

    :param file_name: the name of the file that holds the TLE description of all satellites within the network
    :return: dictionary containing formatted descriptions of satellite network. The 'satellites' key provides
     descriptions used in several functions.
    """

    return read_tles(file_name)
