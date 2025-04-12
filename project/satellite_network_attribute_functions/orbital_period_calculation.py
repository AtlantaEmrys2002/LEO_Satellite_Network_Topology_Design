# Libraries
import ephem


def orbital_period_calculation(satellite_description, num_sat: int) -> float:
    """
    Calculates (using mean motion) the orbital period of a satellite. Mean motion has units in revolutions per day. Uses
     astropy assumption that 1 day == 86400 seconds. Returns the maximum orbital period of any satellite in the network
     (i.e. the orbital period of the satellite orbited at the greatest altitude above the Earth).

    :param satellite_description:
    :param num_sat:
    :return: the maximum length (in seconds) of an orbital period of any satellite within the mega-constellation.
    """
    # For one satellite
    if isinstance(satellite_description, ephem.EarthSatellite):
        return 86400 / satellite_description.n
    else:
        # Return longest orbital period
        return max([86400/satellite_description[x].n for x in range(num_sat)])

# References:
# Astropy Assumptions: https://docs.astropy.org/en/stable/time/index.html
# Mean Motion Formula: https://en.wikipedia.org/wiki/Mean_motion
