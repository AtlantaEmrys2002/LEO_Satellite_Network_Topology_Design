# Libraries
import math
from skyfield.api import EarthSatellite, load, wgs84

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()

earth_radius = 6378.135


def satellite_height_above_earth(name: str, s: str, t: str) -> float:
    """
    Calculates the height (in km) of a satellite above the Earth's surface, assuming the satellite has an approximately
    circular orbit (eccentricity close to 0).

    :param name: the name of the constellation
    :param s: the first line of the TLE coordinate format
    :param t: the second line of the TLE coordinate format
    :return: the height/altitude (in km) of a satellite above Earth
    """
    # Convert TLES to Geocentric Coordinates
    sample_satellite = EarthSatellite(s, t, name, ts)

    # Assume orbital height is approximately the same at all times. Set time to be fixed, e.g. 1st January 2000 Midnight
    # (ensures deterministic behaviour when calculating position, as opposed to using ts.now()). Find satellite position
    # at given time (as circular orbit, should not matter at what time the satellite position is recorded and retrieve
    # height above earth value in km)

    height_above_earth = float(wgs84.height_of(sample_satellite.at(ts.tdb(2000, 1, 1, 0, 0))).km)

    return height_above_earth


def maximum_communication_distance(data_file: str, num_sat: int) -> float:
    """
    Calculates the maximum possible communication distance between two satellites orbiting at the lowest possible orbit
    permitted in the satellite constellation.

    :param data_file: the location of the file containing satellite descriptions in TLE format
    :param num_sat: the number of satellites in the network
    :return: the maximum possible communication distance between two satellites
    """
    # Get sample satellite TLES coordinates to calculate maximum communication distance
    with open(data_file, 'r') as f:
        lines = [line.strip() for line in f]

    # Get data needed for calculating satellite position and calculate each satellite's altitude. Find the lowest
    # orbital altitude - this will provide the "smallest" maximum communication distance, as closer to Earth.
    lowest_satellite_altitude = min([satellite_height_above_earth(lines[line_index], lines[line_index + 1],
                                                                  lines[line_index + 2]) for line_index in
                                     range(1, num_sat, 3)])

    # Return the maximum communication distance between two satellites orbiting at the lowest satellite altitude in the
    # network
    return 2 * math.sqrt(pow(earth_radius + lowest_satellite_altitude, 2) - pow(earth_radius, 2))

# References:
# Astropy Documentation - https://docs.astropy.org/en/latest/coordinates/satellites.html
# NSGA-III Paper (see report references)
# Sources of Earth Radius Constants:
# - https://en.wikipedia.org/wiki/Earth_radius
# - https://github.com/AtlantaEmrys2002/hypatia/blob/master/paper/satellite_networks_state/main_starlink_550.py
