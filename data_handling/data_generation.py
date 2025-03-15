# Libraries
import generate_tles_from_scratch as hypatia_data
import os

# Function calls on Hypatia software function generate_tles_from_scratch_with_sgp to build physical description of
# satellite network using TLE coordinate system and stores description in temporary file. Orbit eccentricity must be > 0
# (according to Hypatia) - set close to 0 to ensure orbit approx. circular
def data_generation(file_name:str, constellation_name:str, num_orbits:int, num_sats_per_orbit:int, inclination_degree:float, mean_motion_rev_per_day, phase_diff=True, eccentricity=0.0000001, arg_of_perigee_degree=0.0):
    """
    Function calls on Hypatia function to build physical description of satellite network using TLE coordinate system and stores description in temporary file. Orbit eccentricity must be > 0 (according to Hypatia) - set close to 0 to ensure orbit is approx. circular.

    :param file_name:
    :param constellation_name:
    :param num_orbits:
    :param num_sats_per_orbit:
    :param inclination_degree:
    :param mean_motion_rev_per_day:
    :param phase_diff:
    :param eccentricity:
    :param arg_of_perigee_degree:
    :return:
    """
    # Check file does not already exist
    if os.path.isfile("./" + file_name) is False:
        # Generate TLE coordinates of satellites
        hypatia_data.generate_tles_from_scratch_with_sgp(file_name, constellation_name, num_orbits, num_sats_per_orbit,
                                                         phase_diff, inclination_degree, eccentricity,
                                                         arg_of_perigee_degree, mean_motion_rev_per_day)
    return