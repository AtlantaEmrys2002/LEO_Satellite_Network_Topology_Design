def maximum_transmission_distance(constellation_name: str) -> float:
    """
    Returns the maximum transmission distance for a given constellation (based on satellite hardware specifications).
    :param constellation_name: name of satellite network constellation, e.g. Starlink-550
    :return: the maximum transmission distance (as dictated by hardware specifications) between satellites.
    """
    # Constellation name
    name = constellation_name.lower()

    if "starlink" in name:
        return 10000
    elif "telesat" in name:
        return 10000
    elif "kuiper" in name:
        return 10000
    elif "iridium" in name:
        return 10000
    else:
        raise ValueError("No constellation has that name and, therefore, no maximum transmission distance can be "
                         "returned.")
