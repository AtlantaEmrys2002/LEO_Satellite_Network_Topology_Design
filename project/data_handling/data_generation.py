# Libraries
from .generate_tles_from_scratch import generate_tles_from_scratch_with_sgp
import os


def data_generation(file_name: str, constellation_name: str, num_orbits: int, num_sats_per_orbit: int,
                    inclination_degree: float, mean_motion_rev_per_day, phase_diff=True, eccentricity=0.0000001,
                    arg_of_perigee_degree=0.0):
    """
    Function calls on Hypatia function to build physical description of satellite network using TLE coordinate system
    and stores description in temporary file. Orbit eccentricity must be > 0 (according to Hypatia) - set close to 0 to
    ensure orbit is approx. circular.

    :param file_name: the name of the file that holds the TLE description of all satellites within the network
    :param constellation_name: name of satellite network constellation, e.g. Starlink-550
    :param num_orbits: the number of orbits (i.e. orbital paths) of satellites within the network
    :param num_sats_per_orbit: the number of satellites per orbit
    :param inclination_degree: the inclination degree ("angle of orbit") of the satellite network
    :param mean_motion_rev_per_day: the number of orbits of Earth each satellite completes per day
    :param phase_diff: the difference in phase between satellites in adjacent orbits
    :param eccentricity: the eccentricity of the orbit - assume orbits are circular as default
    :param arg_of_perigee_degree: point in orbit at which satellite gets closest to Earth
    """
    # Check file does not already exist
    if os.path.isfile("./" + file_name) is False:
        # Generate TLE coordinates of satellites
        generate_tles_from_scratch_with_sgp(file_name, constellation_name, num_orbits, num_sats_per_orbit, phase_diff,
                                            inclination_degree, eccentricity, arg_of_perigee_degree,
                                            mean_motion_rev_per_day)
