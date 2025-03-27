# Libraries
from skyfield.api import load


# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()


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
