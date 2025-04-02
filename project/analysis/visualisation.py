# Libraries

from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
# from mayavi import mlab
import numpy as np
import os
from skyfield.api import EarthSatellite, load

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()


def read_file(file_name):
    """
    Stores raw TLE data rather than convert it to pyephem format. Please note that this function is adapted from Hypatia
     (read_tles.py in the satgenpy module)

    :param file_name:
    :return:
    """
    tles_data = []
    with open(file_name, 'r') as f:
        _, _ = [int(n) for n in f.readline().split()]
        universal_epoch = None
        i = 0
        for tles_line_1 in f:
            tles_line_2 = f.readline()
            tles_line_3 = f.readline()

            # Retrieve name and identifier
            name = tles_line_1
            sid = int(name.split()[1])
            if sid != i:
                raise ValueError("Satellite identifier is not increasing by one each line")
            i += 1

            # Fetch and check the epoch from the TLES data
            # In the TLE, the epoch is given with a Julian data of yyddd.fraction
            # ddd is actually one-based, meaning e.g. 18001 is 1st of January, or 2018-01-01 00:00.
            # As such, to convert it to Astropy Time, we add (ddd - 1) days to it
            # See also: https://www.celestrak.com/columns/v04n03/#FAQ04
            epoch_year = tles_line_2[18:20]
            epoch_day = float(tles_line_2[20:32])
            epoch = Time("20" + epoch_year + "-01-01 00:00:00", scale="tdb") + (epoch_day - 1) * u.day
            if universal_epoch is None:
                universal_epoch = epoch
            if epoch != universal_epoch:
                raise ValueError("The epoch of all TLES must be the same")

            # Finally, store the satellite information
            tles_data.append([tles_line_1.strip(), tles_line_2.strip(), tles_line_3.strip()])

    return tles_data


def snapshot_time_stamp(time_stamp):
    """
    Given time in orbit (based on snapshot number) in seconds, this function converts the time to TDB time format.
    Calculates a satellite's position in orbit (as time since it passed a given point) and reformats for satellite.
    Assumes 60 minutes in an hour and 60 seconds in a minute (i.e. no leap seconds), etc.
    :param time_stamp:
    :return:
    """
    hours = 0
    minutes = 0

    if time_stamp // 3600 >= 1:
        hours = time_stamp // 3600
        time_stamp -= (hours * 3600)
    if time_stamp // 60 >= 1:
        minutes = time_stamp // 60
        time_stamp -= (minutes * 60)
    seconds = time_stamp

    return ts.tdb(2000, 1, 1, hours, minutes, seconds)


def generate_Earth():

    radius_of_earth = 6378.135

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = radius_of_earth * np.outer(np.cos(u), np.sin(v))
    y = radius_of_earth * np.outer(np.sin(u), np.sin(v))
    z = radius_of_earth * np.outer(np.ones(np.size(u)), np.cos(v))

    return x, y, z


def visualise_static(location, tle_file, num_snapshot=94, snapshot_interval=60):

    # INPUTS #

    # Check if any ISL topology exists
    if os.path.isfile(location + "/isls_0.txt") is False:
        raise ValueError("At least one ISL topology must have been built in order to calculate a network's link churn.")

    # Read in topology built for given snapshot
    isls = np.loadtxt(location + "/isls_0.txt").astype(int).T

    # CREATE FIGURE #
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
    # mlab.clf()

    # MODEL EARTH

    # Generate coordinates
    x, y, z = generate_Earth()

    # Plot surface
    ax.plot_surface(x, y, z, antialiased=True, alpha=1)
    ax.plot_surface(x, y, z, antialiased=True, alpha=1)
    ax.plot_surface(x, y, z, antialiased=True, alpha=1)
    # mlab.mesh(x, y, z)

    # GET SATELLITE DESCRIPTIONS #
    # Get TLEs-formatted data
    tles_data = read_file(tle_file)

    # Convert TLES data to Skyfield EarthSatellite objects - used to convert satellite positions to geocentric
    # coordinates - all measurements in km
    earth_satellite_objects = [EarthSatellite(satellite_i[1], satellite_i[2], satellite_i[0], ts) for satellite_i in
                               tles_data]

    # TIMES #

    # Calculate times at which satellite positions observed
    snapshot_times = [snapshot_time_stamp(snapshot_interval * k) for k in range(num_snapshot)]

    for k in snapshot_times:

        # Positions of satellites at each snapshot time
        sats = np.asarray([i.at(k).position.km for i in earth_satellite_objects])

        satellite_xs = sats.T[0].tolist()
        satellite_ys = sats.T[1].tolist()
        satellite_zs = sats.T[2].tolist()

        ax.scatter(satellite_xs, satellite_ys, satellite_zs, marker='o', s=2, color='red')

        # mlab.points3d(satellite_xs, satellite_ys, satellite_zs, scale_factor=0.05)

        # for point in range(len(sats)):
        #     ax.plot(satellite_xs[point], satellite_ys[point], satellite_zs[point], marker='o', markersize=2,
        #             color='red')

        break

    # Set aspect ratio
    ax.set_aspect('equal')

    plt.show()
    # mlab.show()


visualise_static("./Results/plus_grid/kuiper-630", "kuiper-constellation_tles.txt.tmp")
