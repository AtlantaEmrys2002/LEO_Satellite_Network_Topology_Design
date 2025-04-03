# Libraries

from astropy.time import Time
from astropy import units as u
import copy
import matplotlib.pyplot as plt
# from mayavi import mlab
import numpy as np
import os
import plotly.graph_objects as go
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
    d = np.pi / 32

    radius_of_earth = 6378.135

    theta, phi = np.mgrid[0:np.pi + d:d, 0:2 * np.pi:d]

    # Convert to Cartesian coordinates
    x = radius_of_earth * np.sin(theta) * np.cos(phi)
    y = radius_of_earth * np.sin(theta) * np.sin(phi)
    z = radius_of_earth * np.cos(theta)
    points = np.vstack([x.ravel(), y.ravel(), z.ravel()])

    return points


def visualise_static_plotly(location, tle_file, num_snapshot=94, snapshot_interval=60, constellation_name="Kuiper-630"):

    # INPUTS #

    # Check if any ISL topology exists
    if os.path.isfile(location + "/isls_0.txt") is False:
        raise ValueError("At least one ISL topology must have been built in order to calculate a network's link churn.")

    # Read in topology built for given snapshot
    isls = np.loadtxt(location + "/isls_0.txt").astype(int).T

    # CREATE FIGURE #
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # MODEL EARTH

    x, y, z = generate_Earth()

    fig = go.Figure(data=[go.Scatter3d(x=[], y=[], z=[], marker=dict(color="red", size=2)),
                          go.Mesh3d(x=x, y=y, z=z, color='lightblue', opacity=1.0, alphahull=0, name="Earth")])

    # fig.add_traces = go.Scatter3d(x=[], y=[], z=[])

    # Format figure
    fig.update_layout(paper_bgcolor="black", title_font_color="white")

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

    frames = []

    for k in snapshot_times:

        # Positions of satellites at each snapshot time
        sats = np.asarray([i.at(k).position.km for i in earth_satellite_objects]).T.tolist()

        # satellites = go.Scatter3d(x=sats[0], y=sats[1], z=sats[2], mode='markers',
        #                           name=constellation_name)
        #
        # satellites.update(marker_size=2)

        sats = np.asarray(sats).T

        # Calculate the isl positions
        isl_a = [sats[int(isls[0, link])].tolist() for link in range(len(isls[0]))]
        isl_b = [sats[int(isls[1, link])].tolist() for link in range(len(isls[0]))]

        links = []

        for x in range(len(isls[0])):
            links.append(isl_a[x])
            links.append(isl_b[x])
            links.append([None, None, None])

        links = links[:-1]

        links = np.asarray(links).T

        satellites = go.Scatter3d(x=links[0], y=links[1], z=links[2], line=dict(color='red'), name=constellation_name)

        # satellites.update(marker_size=2)

        # frames.append(go.Frame(data=[satellites], traces=[0], name=f'frame{k}'))
        frames.append(go.Frame(data=[satellites], traces=[0], name=f'frame{k}'))

    fig.update(frames=frames)

    fig.update_layout(title=constellation_name, updatemenus=[dict(type="buttons", buttons=[dict(label="Play",
                                                                                                method="animate",
                                                                                                args=[None])])],
                      scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)))

    fig.show()


visualise_static_plotly("./Results/plus_grid/kuiper-630", "kuiper-constellation_tles.txt.tmp", 94, 60, "Kuiper-630")
