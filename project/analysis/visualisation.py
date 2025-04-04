# Libraries

from astropy.time import Time
from astropy import units as u
import geopandas as gpd
import numpy as np
import os
import plotly.graph_objects as go
from skyfield.api import EarthSatellite, load
import zipfile

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()

# Global variable
radius_of_earth = 6378.135


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


def plot_back(fig):
    """back half of sphere"""
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
    """front half of sphere"""
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


def visualise_static(location, tle_file, num_snapshot=94, snapshot_interval=60, constellation_name="Kuiper-630",
                     topology_type="static"):
    
    # INPUTS #

    # Check if any ISL topology exists
    if os.path.isfile(location + "/isls_0.txt") is False:
        raise ValueError("At least one ISL topology must have been built in order to calculate a network's link churn.")

    # Read in topology built for given snapshot
    isls = np.loadtxt(location + "/isls_0.txt").astype(int).T

    # MODEL EARTH

    fig = go.Figure(data=[go.Scatter3d(x=[], y=[], z=[], marker=dict(color="red", size=2), line=dict(color="red"),
                                       name=constellation_name)])

    # Used to download dataset for Earth surface texture
    if os.path.isfile("ne_110m_admin_0_countries.shp") is False:
        # Data File
        with zipfile.ZipFile("./analysis/archive.zip", 'r') as zip_ref:
            zip_ref.extractall("./analysis")

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

    for k in snapshot_times:

        # Positions of satellites at each snapshot time
        sats = np.asarray([i.at(k).position.km for i in earth_satellite_objects]).tolist()

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

    fig.write_html("Results/visualisation" + "/" + constellation_name + "_plus_grid" + ".html")

    fig.show()

    # REMOVE TEXTURE FILES #
    os.remove('./analysis/ne_110m_admin_0_countries.cpg')
    os.remove('./analysis/ne_110m_admin_0_countries.dbf')
    os.remove('./analysis/ne_110m_admin_0_countries.prj')
    os.remove('./analysis/ne_110m_admin_0_countries.README.html')
    os.remove('./analysis/ne_110m_admin_0_countries.shp')
    os.remove('./analysis/ne_110m_admin_0_countries.shx')
    os.remove('./analysis/ne_110m_admin_0_countries.VERSION.txt')


# Run visualisation function
visualise_static("./Results/plus_grid/kuiper-630", "kuiper-constellation_tles.txt.tmp", 94, 60, "Kuiper-630")
