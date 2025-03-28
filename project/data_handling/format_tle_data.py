# Libraries
# import read_tles as hypatia_read_data
from .read_tles import read_tles


def format_tle_data(file_name: str):
    """
    # Function calls on Hypatia read_tles to build dictionary containing the total number of orbits in the
    constellation, as well as the number of satellites per orbit, the epoch of the network, and a description of each
    satellite's position.

    :param file_name:
    :return:
    """
    # return hypatia_read_data.read_tles(file_name)
    return read_tles(file_name)
