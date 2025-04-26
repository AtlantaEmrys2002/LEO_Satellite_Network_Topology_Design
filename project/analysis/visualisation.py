# Libraries

import argparse
from astropy.time import Time
from astropy import units as u
import geopandas as gpd
import numpy as np
import os
import plotly.graph_objects as go
from skyfield.api import EarthSatellite, load
import zipfile

# Please note that these tutorials/sources were heavily utilised in the construction of this visualisation feature
# Dataset from Kaggle - https://www.kaggle.com/datasets/poznyakovskiy/natural-earth-1110m-countries (used for texturing)
# Texture Application - https://community.plotly.com/t/applying-full-color-image-texture-to-create-an-interactive-earth-
# globe/60166
# Texture Application - https://community.plotly.com/t/create-earth-sphere-with-all-countries-in-plotly/79284

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()

# Global variable
radius_of_earth = 6378.135


def read_file(file_name: str):
    """
    Stores raw TLE data rather than convert it to pyephem format. Please note that this function is adapted from Hypatia
    (read_tles.py in the satgenpy module). This is a copy of a function (the original is in the data_handling
    directory).

    :param file_name: name of TLEs file in which description of each satellite in mega-constellation situated.
    :return: TLEs data
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


def snapshot_time_stamp(time_stamp: float) -> Time:
    """
    Given time in orbit (based on snapshot number) in seconds, this function converts the time to TDB time format.
    Calculates a satellite's position in orbit (as time since it passed a given point) and reformats for satellite.
    Assumes 60 minutes in an hour and 60 seconds in a minute (i.e. no leap seconds), etc.
    :param time_stamp:
    :return: formatted timestamp (TDB)
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


def plot_back(fig):
    """Sourced from here: https://community.plotly.com/t/create-earth-sphere-with-all-countries-in-plotly/79284"""
    clor = f'rgb(220, 220, 220)'
    R = np.sqrt(6368.134)
    u_angle = np.linspace(0, np.pi, 25)
    v_angle = np.linspace(0, np.pi, 25)
    x_dir = np.outer(R * np.cos(u_angle), R * np.sin(v_angle))
    y_dir = np.outer(R * np.sin(u_angle), R * np.sin(v_angle))
    z_dir = np.outer(R * np.ones(u_angle.shape[0]), R * np.cos(v_angle))
    fig.add_surface(z=z_dir, x=x_dir, y=y_dir, colorscale=[[0, clor], [1, clor]], opacity=1.0, showlegend=False,
                    lighting=dict(diffuse=0.1), name="Earth")


def plot_front(fig):
    """Sourced from here: https://community.plotly.com/t/create-earth-sphere-with-all-countries-in-plotly/79284"""
    clor = f'rgb(220, 220, 220)'
    R = np.sqrt(radius_of_earth)
    u_angle = np.linspace(-np.pi, 0, 25)
    v_angle = np.linspace(0, np.pi, 25)
    x_dir = np.outer(R * np.cos(u_angle), R * np.sin(v_angle))
    y_dir = np.outer(R * np.sin(u_angle), R * np.sin(v_angle))
    z_dir = np.outer(R * np.ones(u_angle.shape[0]), R * np.cos(v_angle))
    fig.add_surface(z=z_dir, x=x_dir, y=y_dir, colorscale=[[0, clor], [1, clor]], opacity=1.0, showlegend=False,
                    lighting=dict(
                        diffuse=0.1), name="Earth")


def plot_polygon(poly):
    """
    Sourced from here: https://community.plotly.com/t/create-earth-sphere-with-all-countries-in-plotly/79284
    :param poly:
    :return:
    """
    xy_coords = poly.exterior.coords.xy
    lon = np.array(xy_coords[0])
    lat = np.array(xy_coords[1])

    lon = lon * np.pi / 180
    lat = lat * np.pi / 180

    R = radius_of_earth
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)

    return x, y, z


def visualise(location: str, tle_file: str, num_snapshot: int = 94, snapshot_interval: float = 60,
              constellation_name: str = "Kuiper-630", topology_type: str = "static",
              topology_method: str = "plus_grid"):
    """
    Used to visualise constellations (based on topologies constructed by algorithms).
    :param location: location of ISL topology files within file system.
    :param tle_file: the name of the file that stores the TLE descriptions of each satellite within the
     mega-constellation.
    :param num_snapshot: the number of snapshots of the network over one orbit for which a topology is constructed
    :param snapshot_interval: the time (in seconds) between snapshots
    :param constellation_name: name of satellite network constellation, e.g. Starlink-550
    :param topology_type: 'static' or 'dynamic' - indicates whether satellite topology evolves over time (i.e. ISLs
     connect and disconnect over the course of one orbit)
    :param topology_method: the algorithm with which the topology was constructed (plus_grid, x_grid, mdtd, or novel)
    """
    # INPUTS #

    # Check if any ISL topology exists
    if os.path.isfile(location + "/isls_0.txt") is False:
        raise ValueError("At least one ISL topology must have been built in order to calculate a network's "
                         "link churn.")

    # Read in topology built for given snapshot
    # isls = np.loadtxt(location + "/isls_0.txt").astype(int).T

    # MODEL EARTH

    fig = go.Figure(data=[go.Scatter3d(x=[], y=[], z=[], marker=dict(color="red", size=2), line=dict(color="red"),
                                       name=constellation_name)])

    # Used to download dataset for Earth surface texture
    if os.path.isfile("ne_110m_admin_0_countries.shp") is False:
        # Data File
        with zipfile.ZipFile("./analysis/archive.zip", 'r') as zip_ref:
            zip_ref.extractall("./analysis")

    # Read in texture file
    gdf = gpd.read_file("./analysis/ne_110m_admin_0_countries.shp")
    plot_front(fig)
    plot_back(fig)

    for i in gdf.index:

        polys = gdf.loc[i].geometry

        if polys.geom_type == 'Polygon':
            x, y, z = plot_polygon(polys)
            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color=f'rgb(0, 0,0)'), showlegend=False,
                                       name="Earth"))

        elif polys.geom_type == 'MultiPolygon':

            for poly in polys.geoms:
                x, y, z = plot_polygon(poly)
                fig.add_trace(
                    go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color=f'rgb(0, 0,0)'), showlegend=False,
                                 name="Earth"))

    # FORMAT FIGURE #
    fig.layout.template = "plotly"
    fig.update_layout(paper_bgcolor="black", title_font_color="white", font_color="white")

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

    time_value = 0

    counter = 0

    for k in snapshot_times:

        # Positions of satellites at each snapshot time
        sats = np.asarray([i.at(k).position.km for i in earth_satellite_objects]).tolist()

        if topology_type == 'static':

            # Check if any ISL topology exists
            if os.path.isfile(location + "/isls_0.txt") is False:
                raise ValueError(
                    "At least one ISL topology must have been built in order to calculate a network's link churn.")

            # Read in topology built for given snapshot
            isls = np.loadtxt(location + "/isls_0.txt").astype(int).T

        else:

            # Check if any ISL topology exists
            if os.path.isfile(location + "/isls_" + str(counter) + ".txt") is False:
                raise ValueError(
                    "That topology does not exist.")

            # Read in topology built for given snapshot
            isls = np.loadtxt(location + "/isls_" + str(counter) + ".txt").astype(int).T

        # Calculate the isl positions
        isl_a = [sats[int(isls[0, link])] for link in range(len(isls[0]))]
        isl_b = [sats[int(isls[1, link])] for link in range(len(isls[0]))]

        links = []

        for x in range(len(isls[0])):
            links.append(isl_a[x])
            links.append(isl_b[x])
            links.append([None, None, None])

        links = np.asarray(links[:-1]).T

        # Draw links
        satellites = go.Scatter3d(x=links[0], y=links[1], z=links[2])

        frames.append(go.Frame(data=[satellites], traces=[0], name=str(time_value) + 's'))

        time_value += snapshot_interval
        counter += 1

    fig.frames = frames

    # SLIDER AND ANIMATION #

    fig.update_layout(title=constellation_name,
                      scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)),
                      updatemenus=[
                          {
                              "showactive": False,
                              "buttons": [{"args": [None, {"frame": {"duration": 300, "redraw": True}}],
                                           "label": "Play", "method": "animate"},
                                          {"args": [[None], {"frame": {"duration": 0, "redraw": False},
                                                             "mode": "immediate", "transition": {"duration": 0}, }, ],
                                           "label": "Pause", "method": "animate", }, ],
                              "type": "buttons",
                          }
                      ],
                      sliders=[{"steps": [{"args": [[f.name], {"frame": {"duration": 0, "redraw": True},
                                                               "mode": "immediate", }, ],
                                           "label": f.name, "method": "animate", }
                                          for f in frames], }])

    fig.data[0].visible = True

    if os.path.isdir("Results/visualisation") is False:
        try:
            os.makedirs("Results/visualisation")
        except OSError:
            print("Directory to store visualisations could not be created.")

    fig.write_html("Results/visualisation" + "/" + constellation_name.lower() + "_" + topology_method + ".html")

    fig.show()

    # REMOVE TEXTURE FILES #
    os.remove('./analysis/ne_110m_admin_0_countries.cpg')
    os.remove('./analysis/ne_110m_admin_0_countries.dbf')
    os.remove('./analysis/ne_110m_admin_0_countries.prj')
    os.remove('./analysis/ne_110m_admin_0_countries.README.html')
    os.remove('./analysis/ne_110m_admin_0_countries.shp')
    os.remove('./analysis/ne_110m_admin_0_countries.shx')
    os.remove('./analysis/ne_110m_admin_0_countries.VERSION.txt')


# RUN VISUALISATION CODE #

# Run visualisation function

if __name__ == "__main__":

    print("Building topology visualisation... ")

    # Parse inputs to module
    parser = argparse.ArgumentParser()

    parser.add_argument("--location", type=str, help="location of files containing topology",
                        required=True)

    parser.add_argument("--tles", type=str, help="name of tle file containing TLE description of satellite"
                                                 " positions ", required=True)

    parser.add_argument("--num_snapshots", type=int, help="number of constellation snapshots for which "
                                                         "dynamic topology was built", required=True)

    parser.add_argument("--snapshot_interval", type=float, help="time (in seconds) between snapshots",
                        required=True)

    parser.add_argument("--constellation_name", type=str, help="name of constellation to be visualised",
                        required=True)

    parser.add_argument("--topology-type", type=str, help="options: 'static' or 'dynamic' - indicates "
                                                          "whether ISLs change over time (dynamic)", required=True)

    parser.add_argument("--method", type=str, help="method with which topology was constructed - options "
                                                   "include: plus_grid, mdtd, novel", required=True)

    args = parser.parse_args()

    location = args.location
    tle_file = args.tles
    num_snapshots = args.num_snapshots
    snapshot_interval = args.snapshot_interval
    constellation_name = args.constellation_name
    topology_type = args.topology_type
    method = args.method

    visualise(location, tle_file, num_snapshots, snapshot_interval, constellation_name, topology_type, method)

# visualise("./Results/plus_grid/kuiper-630", "kuiper-constellation_tles.txt.tmp", 97,
#                  60, "Kuiper-630", topology_type="static", topology_method="plus_grid")
#
# visualise("./Results/plus_grid/starlink-550", "starlink-constellation_tles.txt.tmp", 90,
#                  60, "Starlink-550", topology_method="plus_grid")
#
# visualise("./Results/plus_grid/telesat-1015", "telesat-constellation_tles.txt.tmp", 105,
#                  60, "Telesat-1015", topology_method="plus_grid")


# print("Completed")
#
# print("Building visualisations of dynamic topologies... ", end='\r')
#
# # MDTD Benchmark
#
# visualise("./Results/mdtd/kuiper-630", "kuiper-constellation_tles.txt.tmp", 97,
#           60, "Kuiper-630", topology_type="dynamic", topology_method="mdtd")
#
# visualise("./Results/mdtd/starlink-550", "starlink-constellation_tles.txt.tmp", 90,
#           60, "Starlink-550", topology_type="dynamic", topology_method="mdtd")
#
# visualise("./Results/mdtd/telesat-1015", "telesat-constellation_tles.txt.tmp", 105,
#           60, "Telesat-1015", topology_type="dynamic", topology_method="mdtd")
#
    print("Completed")


# visualise("./Results/novel/random/primal/kuiper-630", "kuiper-constellation_tles.txt.tmp", 90,
#           60, "Kuiper-630", topology_type="dynamic", topology_method="novel")

# visualise("./novel/primal/multishell", "multishell-constellation_tles.txt.tmp", 97,
#                  60, "Multishell", topology_type="dynamic", topology_method="novel")

# References
# Animation of 3d Plots - https://stackoverflow.com/questions/68100031/animated-3d-surface-plots-with-plotly
# Animation of 3d Scatter - https://community.plotly.com/t/3d-scatter-animation/46368
# Animation of 3d Surfaces - https://stackoverflow.com/questions/76363493/plotly-animated-3d-surface-plots
# Animating Traces - https://community.plotly.com/t/how-to-simultaneously-animate-multiple-traces-of-the-same-figure/
# 64541
# Arrows - https://stackoverflow.com/questions/43164909/plotlypython-how-to-plot-arrows-in-3d
# Callback to Frames - https://community.plotly.com/t/dash-animated-3d-graph-with-callback-to-change-animation-frames/
# 46939
# Complex Animations in Plotly - https://community.plotly.com/t/how-to-use-plotly-python-for-complex-3d-animations/32983
# Cylindrical Maps - https://stackoverflow.com/questions/17365364/plotting-cylindrical-map-data-over-a-3d-sphere
# Downloading from Kaggle - https://ravi-chan.medium.com/how-to-download-any-data-set-from-kaggle-7e2adc152d7f
# Overview - https://medium.com/@alexeyyurasov/3d-modeling-with-python-c21296756db2
# Matplotlib Opacity - https://stackoverflow.com/questions/15794499/how-can-i-plot-a-graph-on-an-opaque-surface-in-
# matplotlib
# Matplotlib Option - https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html
# Plotting Orbits with Matplotlib - https://plainenglish.io/blog/plot-satellites-real-time-orbits-with-python-s-
# matplotlib
# Plotly Background Colour - https://community.plotly.com/t/how-can-i-change-the-background-color/39007
# Plotly Bug - https://github.com/plotly/plotly.py/issues/3120#issuecomment-1018308604
# Plotly Button Colour - https://stackoverflow.com/questions/75371336/how-do-i-change-the-color-of-a-button-color-in-
# plotly-for-python
# Plotly Dark Theme - https://community.plotly.com/t/plotly-default-buttons-dark-theme/72412/2
# Plotly Documentation - https://plotly.com/graphing-libraries/
# Plotly Recommendation - https://stackoverflow.com/questions/31768031/plotting-points-on-the-surface-of-a-sphere
# Plotly Removing Axes - https://community.plotly.com/t/how-to-remove-axis-planes-in-plotly-graph-objects-figure/67375
# Plotly Sliders - https://stackoverflow.com/questions/68667116/add-slider-to-plotly-heatmap-animation-python
# Plotly Switch - https://community.plotly.com/t/how-to-change-color-of-dbc-switch-button/63389/3
# Getting Started with Plotly - https://plotly.com/python/getting-started/#overview
# Plotly Textures - https://community.plotly.com/t/applying-full-color-image-texture-to-create-an-interactive-earth-
# globe/60166
# Plotly 3D Surface Bug - https://stackoverflow.com/questions/74470382/plotly-3d-surface-plot-not-appearing
# Plotting 3D Lines - https://stackoverflow.com/questions/60153352/plotting-multiple-3d-lines-in-one-figure-using-plotly
# Plotting 3D Spheres - https://stackoverflow.com/questions/70977042/how-to-plot-spheres-in-3d-with-plotly-or-another-
# library
# Plotting on Mesh - https://stackoverflow.com/questions/72650425/python-plotly-scatter-plot-on-3d-mesh
# S3DLib Documentation - https://s3dlib.org/examples/animations/anim_bases.html
# 3D Plotly - https://plotly.com/python/3d-scatter-plots/
# Plotly Documentation - https://plotly.com/python/animations/
# Plotting 3D Spheres in General - https://stackoverflow.com/questions/70977042/how-to-plot-spheres-in-3d-with-plotly-or
# -another-library
# 3D Spheres Plotly - https://stackoverflow.com/questions/70977042/how-to-plot-spheres-in-3d-with-plotly-or-another-
# library
# Satellite Visualisation - https://medium.com/spatial-data-science/8-essential-python-libraries-for-satellite-data-
# visualization-in-geospatial-analysis-ed757d4964f1
# Select Elements from Numpy Array - https://www.quora.com/How-do-I-select-elements-from-a-NumPy-array-in-Python
# Trace Labels - https://community.plotly.com/t/add-trace-s-name-in-legend-with-plotly-express-when-color-is-not-
# specified/27634
