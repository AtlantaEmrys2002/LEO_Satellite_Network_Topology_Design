#TODO: test data
# unit tests for each function
# implement all functions,
# find all input parameters for all functions (may need more than generate_plus_grid does for main function heuristic_topology_design_algorithm_isls)
# this function will be placed in hypatia/satgenpy/satgen/isls
# properly document your code with doxygen or similar

# Libraries

def satellite_visibility_function():
    pass

def time_visibility_function():
    pass

def cost_function():
    pass

def distance_function():
    pass

def degree_constrained_mst():
    pass

def conversion_from_dcmst_to_list():
    pass

def increase_connectivity():
    pass

def heuristic_topology_design_algorithm_isls(output_filename_isls, n_orbits, n_sats_per_orbit, orbital_period, num_snapshot, isl_shift, idx_offset=0):

    # Check satellite network has a sufficient number of satellites and orbits
    if n_orbits < 3 or n_sats_per_orbit < 3:
        raise ValueError("Number of x and y must each be at least 3")

    # Calculate the time between each snapshot
    snapshot_interval = orbital_period/num_snapshot

    for snapshot in range(0, orbital_period, snapshot_interval):

        # Initialise list of satellite pairs for snapshot k
        list_isls = []

        # Calculate distance matrix
        # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE TIME VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
        distance_matrix = distance_function()

        # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
        # Calculate visibility matrix
        visibility_matrix = satellite_visibility_function()

        # Calculate time visibility matrix
        # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE TIME VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
        time_visibility_matrix = time_visibility_function()

        # Calculate cost matrix
        # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE TIME VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
        cost_matrix = cost_function()

        # Calculate degree constrained minimum spanning tree
        # NOTE TO SELF - DETERMINE INPUTS TO SATELLITE TIME VISIBILITY FUNC AND ENSURE IT RETURNS A MATRIX
        tree = degree_constrained_mst()

        # Convert DCMST to list of ISLs
        list_isls = conversion_from_dcmst_to_list()

        # Adding Edges
        list_isls = increase_connectivity()


        # Write results to file for snapshot k
        with open(output_filename_isls + "_snapshot_" + str(snapshot), 'w+') as f:
            for (a, b) in list_isls:
                f.write(str(a) + " " + str(b) + "\n")

        pass

# References:
# Hypatia - https://github.com/snkas/hypatia/tree/master