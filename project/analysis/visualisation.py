# Libraries
import data_handling as data_handling
import matplotlib.pyplot as plt
import numpy as np
import os
from skyfield.api import EarthSatellite, load

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()


def generate_Earth():

    radius_of_earth = 6378.135

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = radius_of_earth * np.outer(np.cos(u), np.sin(v))
    y = radius_of_earth * np.outer(np.sin(u), np.sin(v))
    z = radius_of_earth * np.outer(np.ones(np.size(u)), np.cos(v))

    return x, y, z


def visualise_static(location, tle_file):

    # INPUTS #

    # Check if any ISL topology exists
    if os.path.isfile(location + "/isls_0.txt") is False:
        raise ValueError("At least one ISL topology must have been built in order to calculate a network's link churn.")

    # Read in topology built for given snapshot
    isls = np.loadtxt(location + "/isls_0.txt").astype(int).T

    # GET SATELLITE DESCRIPTIONS #
    # Get TLEs-formatted data
    tles_data = data_handling.read_file(tle_file)

    # Convert TLES data to Skyfield EarthSatellite objects - used to convert satellite positions to geocentric
    # coordinates - all measurements in km
    earth_satellite_objects = [EarthSatellite(satellite_i[1], satellite_i[2], satellite_i[0], ts) for satellite_i in
                               tles_data]

    # CREATE FIGURE #
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # MODEL EARTH

    # Generate coordinates
    x, y, z = generate_Earth()

    # Plot surface
    ax.plot_surface(x, y, z)

    # Set aspect ratio
    ax.set_aspect('equal')

    # MODEL SATELLITES #
    k = 0
    satellites_at_k = [i.at(k).position.km for i in earth_satellite_objects]

    print(satellites_at_k)

    plt.show()


visualise_static("./Results/plus_grid/kuiper-630", "kuiper-constellation_tles.txt.tmp")
