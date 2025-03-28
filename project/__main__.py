#TODO
# Find maximum transmission dist for Starlink, Kuiper and telesat - 27000 paper sets at 5014 km - also mentioned in code
# EXPERIMENT WITH INCREASE CONNECTIVITY FUNC - IN FREE OPTICAL SPACE NETWORKS PAPER (BROADBAND NOT SATELLITE), THEY
# CONNECT LARGEST COST EDGES - REDUCES GRAPH DIAMETER AT EXPENSE OF ENERGY EFFICIENCY
# When make modifications to documentation: look at this first: https://medium.com/@pratikdomadiya123/build-project-docu
# mentation-quickly-with-the-sphinx-python-2a9732b66594

# Libraries
import argparse
from multiprocessing import Pool
import numpy as np
import os
import random
from scipy.spatial.distance import cdist
from skyfield.api import EarthSatellite, load  # recommended by astropy documentation
import time

# Implemented Modules
from analysis import measure
import cost_function_optimisation_algorithms
import data_handling
import satellite_network_attribute_functions as satnet
from satellite_topology_construction_algorithms import (heuristic_topology_design_algorithm_isls, plus_grid,
                                                        minimum_delay_topology_design_algorithm)

# Global Timescale (used for determining how to calculate time with Skyfield functions)
ts = load.timescale()

# Sun's Ephemeris - used to calculate whether satellite is in sunlight or not
eph = load('./de421.bsp')

# Seed Random so results can be reproduced
random.seed(42)

start = time.time()
if __name__ == "__main__":

    # PARSE INPUTS #
    print('\n')

    print("Parsing Inputs...", end='\r', flush=True)

    # Parse inputs to module
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("--tles", type=str, help="name of input file which contains TLE description of "
                                                 "satellite network", required=True)
    parser.add_argument("--constellation", type=str, help="name of satellite constellation for which a "
                                                          "topology is built (used to name output files)",
                        required=True)
    parser.add_argument("--m", type=int, help="number of orbits in constellation", required=True)
    parser.add_argument("--n", type=int, help="number of satellites per orbit", required=True)
    parser.add_argument("--i", type=float, help="inclination degree of an orbit within the constellation",
                        required=True)
    parser.add_argument("--rev", type=float, help="mean motion revolutions per day for satellite network",
                        required=True)
    parser.add_argument("--snapshot_interval", type=float, help="time intervals (s) between snapshots, e.g."
                                                                " if 60, a topology is constructed every 60s over 1 "
                                                                "orbital period", required=True)
    parser.add_argument("--multi", type=str, help="determines if the constellation contains multiple shells",
                        required=True)
    parser.add_argument("--optimise", type=str, help="determines if cost function weights should be "
                                                     "optimised and/or metrics returned", required=True)
    parser.add_argument("--topology", type=str, help="determines method with which to construct topology "
                                                     "for network (options: 'plus-grid', 'x-grid', 'mdtd', 'novel')",
                        required=True)
    parser.add_argument("--isl_terminals", type=int, nargs="+", help="specify as an int or list of ints "
                                                                     "the number of terminals each satellite in the "
                                                                     "given constellation has", required=True)

    # Optional arguments
    parser.add_argument("--snapshots", type=int, nargs="+", help="ids of the snapshots for which to build "
                                                                 "a topology (minimum of 0, maximum is number of "
                                                                 "snapshots taken over the course of one orbital period"
                                                                 ")")
    # Below are only required for novel algorithm
    parser.add_argument("--dcmst", type=str, help="if novel topology construction algorithm is used, "
                                                  "determines which dcmst construction method is used (options: 'aco', "
                                                  "'ga', 'primal')")
    parser.add_argument("--optimisation_method", type=str, help="if cost function is optimised, determines"
                                                                " method with which to find optimal weights (options: "
                                                                "'random', 'evolutionary', 'machine-learning')")
    parser.add_argument("--weights", type=float, nargs=3, help="values of weights of cost function (alpha "
                                                               "for time visibility, beta for distance, gamma for "
                                                               "probability of failure)")

    args = parser.parse_args()

    tle_file = args.tles
    constellation_name = args.constellation
    num_orbits = args.m
    num_sats_per_orbit = args.n
    inclination_degree = args.i
    mean_motion_rev_per_day = args.rev
    multi_shell = args.multi
    topology = args.topology
    snapshot_interval = args.snapshot_interval

    # Determine whether the cost function should be optimised for the novel proposed algorithm and if the topologies
    # generated by other methods should be evaluated according to given metric
    if args.optimise == "True":
        optimise = True
    elif args.optimise == "False":
        optimise = False
    else:
        raise ValueError("Argument for --optimise should be either 'True' or 'False'.")

    # Determine whether the satellite network involves multiple shells of satellites
    if args.multi == "True":
        multi_shell = True
    elif args.multi == "False":
        multi_shell = False
    else:
        raise ValueError("Argument for --multi should be either 'True' or 'False'.")

    # CALCULATE TOTAL SATELLITES

    # Calculate the number of satellites in the network
    total_sat = num_sats_per_orbit * num_orbits

    # DETERMINE IF A TOPOLOGY CAN BE BUILT FOR THE NETWORKS #

    # Check satellite network has a sufficient number of satellites and orbits - inspired by Hypatia code
    if total_sat < 3 or num_sats_per_orbit < 3:
        raise ValueError("Number of satellites must be greater than 3 and number of satellites per orbit must be "
                         "greater than 3.")

    # FIND THE MAXIMUM NUMBER OF ISLS EACH SATELLITE CAN ESTABLISH #

    # Set up degree constraints, i.e. assign the maximum number of ISLs (or maximum number of functioning/active ISL
    # terminals per satellite).

    if len(args.isl_terminals) == 1:
        try:
            satellite_degree_constraints = np.array([int(args.isl_terminals[0]) for _ in range(total_sat)],
                                                    dtype=np.int32)
        except ValueError:
            raise ValueError(
                "Satellite degree constraints must be specified as int or list of ints (length of list of ints"
                " must be equal to the number of satellites within the network.")
    else:
        if len(args.isl_terminals) != total_sat:
            raise ValueError("The number of specified degree constraints should match the number of satellites within "
                             "the network.")
        try:
            satellite_degree_constraints = [int(a) for a in args.isl_terminals]
        except ValueError:
            raise ValueError(
                "Satellite degree constraints must be specified as int or list of ints (length of list of ints"
                " must be equal to the number of satellites within the network.")

    print("Input Parsing Completed")

    print("Building Topology... ", end='\r', flush=True)

    # Generate test data using network description (data from Hypatia)
    data_handling.data_generation(tle_file, constellation_name, num_orbits, num_sats_per_orbit, inclination_degree,
                                  mean_motion_rev_per_day)

    # Read test data into appropriate data structure (dictionary) and extract description of satellite positions,
    # as well as unique orbits, from data
    satellite_data = data_handling.format_tle_data(tle_file)["satellites"]

    # Calculate orbital period of network (or maximum orbital period if satellites orbit at different altitudes -
    # i.e. in multiple shells)
    if multi_shell is False:
        orbital_period = satnet.orbital_period_calculation(satellite_data[0], total_sat)
    else:
        orbital_period = satnet.orbital_period_calculation(satellite_data, total_sat)

    # Find the maximum communication distance between two satellites (varies as satellite altitudes vary)
    max_communication_dist = satnet.maximum_communication_distance(tle_file, total_sat)

    # This is the maximum distance a satellite can establish signal (transmission power) - can differ from max
    # communication distance (which is determined by orbital positioning, rather than satellite hardware
    # specifications)
    max_transmission_dist = satnet.maximum_transmission_distance(constellation_name)

    # Find smaller of these two numbers - the maximum distance before Earth is in the way or the max transmission
    # distance (due to satellite power constraints)
    max_communication_dist = min(max_communication_dist, max_transmission_dist)

    # The number of snapshots over an orbital period for which to construct a topology
    num_snapshot = int(orbital_period / snapshot_interval)

    # Determines which snapshots to build topologies for (if specified, only construct topologies for those
    # snapshots, otherwise construct a topology for all snapshots
    if args.snapshots:
        snapshot_ids = args.snapshots
    else:
        snapshot_ids = list(range(0, num_snapshot))

    # Get TLEs-formatted data
    tles_data = data_handling.read_file(tle_file)

    # Convert TLES data to Skyfield EarthSatellite objects - used to convert satellite positions to geocentric
    # coordinates - all measurements in km
    earth_satellite_objects = [EarthSatellite(satellite_i[1], satellite_i[2], satellite_i[0], ts) for satellite_i in
                               tles_data]

    # CALCULATE NETWORK ATTRIBUTES #

    # Directory in which to store network attributes
    network_attributes_location = "./" + constellation_name

    # DISTANCE MATRICES #

    # Calculate distance matrices - set time t (using TDB time format) to indicate point in orbital period when
    # snapshot taken

    # Calculate the time (in TDB format) at which each snapshot is taken
    snapshot_times = [satnet.snapshot_time_stamp(snapshot_interval * k) for k in range(num_snapshot)]

    # Calculate distance matrices for each snapshot

    # Create directory to store the distance matrix for each snapshot in individual file within directory (can't
    # process all at once otherwise)
    if os.path.isdir(network_attributes_location + "/distance_matrices") is False:

        # Create directory in which to store distance matrices
        try:
            os.makedirs(network_attributes_location + "/distance_matrices")
        except OSError:
            print("Directory to store distance matrices could not be created.")

        # Keep track of file ID
        file_id = 0

        # Calculate the distance matrix (symmetric) for each snapshot of the network. Distance between satellite and
        # itself = 0km.
        for k in snapshot_times:
            # Calculates position of all satellites in the network at snapshot time k
            satellites_at_k = [i.at(k).position.km for i in earth_satellite_objects]

            # Calculate distance (Euclidean) between all satellite pairs i and j in the network at snapshot time k
            dist_matrix = cdist(satellites_at_k, satellites_at_k, metric='euclidean')

            # Save distance matrix to .npy file
            np.save(network_attributes_location + "/distance_matrices/dist_matrix_" + str(file_id) + ".npy",
                    dist_matrix)

            # Increment ID counter
            file_id += 1

    # STATIC ALGORITHMS #
    # Benchmark Static Topology Designs
    if topology == "plus-grid":

        # Build topology with provided parameters

        # Location to store ISL topology
        if optimise is False:
            location = "./plus_grid/" + constellation_name.lower()
        else:
            location = "./Results/plus_grid/" + constellation_name.lower()

        # Check directory for resulting topology exists
        if os.path.isdir(location) is False:
            try:
                os.makedirs(location)
            except OSError:
                print("Directory to store plus grid (+Grid) topology could not be created.")

        # Use Hypatia implementation to create +Grid topology
        plus_grid.generate_plus_grid_isls(location + "/isls_0.txt", num_orbits, num_sats_per_orbit, isl_shift=0,
                                          idx_offset=0)

        print("+Grid Topology Build Completed")

        # Return metrics if optimise is true (so topology can be evaluated)
        if optimise is True:
            print("Evaluating... ", end='\r', flush=True)

            # Calculate metrics for topology
            max_pd, mean_pd, av_hop_count, link_churn = measure.measure_static(constellation_name, location +
                                                                               "/isls_0.txt", total_sat)

            data_handling.write_optimisation_results_to_csv(location, "static", [max_pd, mean_pd,
                                                                                 av_hop_count, link_churn])

            print("+Grid Evaluation Completed")

    elif topology == "x-grid":

        print("NEEDS IMPLEMENTING")

    else:

        # VISIBILITY MATRICES ###

        # Calculate visibility matrices for all snapshots - these visibility matrices are utilised by both the MDTD and
        # novel algorithms

        # Create directory to store the visibility matrix for each snapshot in individual file within directory
        # (can't process all at once otherwise) - within a visibility matrix, element [i][j] is set to 0 if
        # satellites are not visible to one another
        if os.path.isdir(network_attributes_location + "/visibility_matrices") is False:

            # Create directory in which to store distance matrices
            try:
                os.makedirs(network_attributes_location + "/visibility_matrices")
            except OSError:
                print("Directory to store visibility matrices could not be created.")

            # Calculate all visibility matrices
            for k in range(num_snapshot):
                # Calculate visibility matrix for snapshot and save to .npy file - load distance from corresponding
                # distance matrix file
                visibility_matrix = satnet.visibility_function(
                    np.load(network_attributes_location + "/distance_matrices/dist_matrix_" + str(k) + ".npy"),
                    max_communication_dist)

                np.save(network_attributes_location + "/visibility_matrices/visibility_matrix_" + str(k) + ".npy",
                        visibility_matrix)

        # GENERATES TOPOLOGY USING BENCHMARK MDTD ALGORITHM

        # Benchmark MDTD Design
        # if topology == "mdtd":
        if topology == "mdtd":

            # Directory in which to store topologies
            if optimise is False:
                location = "./mdtd/" + constellation_name.lower()
            else:
                location = "./Results/mdtd/" + constellation_name.lower()

            # Build topologies for satellite network
            minimum_delay_topology_design_algorithm(constellation_name, num_snapshot, total_sat,
                                                    satellite_degree_constraints, topology)

            print("MDTD Topology Build Completed")

            # If optimise is true, generate metrics for topology for evaluation/comparison with novel algorithm
            if optimise is True:

                print("Evaluating...", end='\r', flush=True)

                # Evaluate generated topologies according to given metrics
                max_pd, mean_pd, av_hop_count, link_churn = measure.measure_dynamic(constellation_name, location,
                                                                                    total_sat, num_snapshot)

                # Store results accordingly
                data_handling.write_optimisation_results_to_csv(location, "dynamic", [max_pd, mean_pd,
                                                                                      av_hop_count, link_churn])

                print("MDTD Evaluation Completed")

        # NOVEL TOPOLOGY DESIGN ALGORITHM (BUILT FOR THIS PROJECT)
        elif topology == "novel":

            # PARSES OPTIONAL ARGUMENTS #
            # Parses arguments only required by novel algorithm

            # Checks DCMST construction algorithm exists and sets to dcmst
            if args.dcmst:
                dcmst = args.dcmst
            else:
                raise ValueError("DCMST construction method must be specified for novel algorithm.")

            # Directory in which to store topologies
            if optimise is False:
                location = "./novel/" + dcmst + "/" + constellation_name.lower() + "/"
            else:
                # Checks that optimisation method is specified
                if args.optimisation_method:
                    optimisation_method = args.optimisation_method
                    location = ("./Results/novel/" + optimisation_method + "/" + dcmst + "/" +
                                constellation_name.lower() + "/")
                else:
                    raise ValueError("An optimisation method must be specified.")

            # SUNLIGHT MATRICES #

            # Calculate whether satellites are in sunlight (i.e. vulnerable to solar flares) or on the opposite side of
            # the Earth

            if os.path.isdir(network_attributes_location + "/sunlight_matrices") is False:

                # Create directory in which to store distance matrices
                try:
                    os.makedirs(network_attributes_location + "/sunlight_matrices")
                except OSError:
                    print("Directory to store sunlight matrices could not be created.")

                file_id = 0

                # Calculate all distance matrices
                for k in snapshot_times:
                    # Calculate whether each satellite is in sunlight or not
                    satellites_in_sun = [i.at(k).is_sunlit(eph) for i in earth_satellite_objects]

                    # Update matrix such that element sunlight_matrix[i][j] is set to 1 if i or j is in sunlight and
                    # save to file
                    np.save(network_attributes_location + "/sunlight_matrices/sunlight_matrix_" + str(file_id) + ".npy",
                            satnet.sunlight_function(satellites_in_sun, total_sat))

                    file_id += 1

            # Run topology generation algorithm for each specified snapshot - utilise multiprocessing to make program
            # faster (dependent on the number of cores of computer run program on)

            if optimise is False:

                # Checks weights exist and set to params
                if args.weights:
                    params = args.weights
                else:
                    raise ValueError("Cost function weights must be specified for novel algorithm that does not undergo"
                                     " optimisation.")

                # Generate arguments for functions
                snapshot_arguments = [
                    [constellation_name, total_sat, num_snapshot, satellite_degree_constraints, t, params, location,
                     dcmst] for t in snapshot_ids]

                # Generate topology
                pool = Pool(processes=os.cpu_count())

                pool.map(heuristic_topology_design_algorithm_isls, snapshot_arguments)

                pool.terminate()

                print("Novel Algorithm Topology Build Completed")

            # Run cost optimisation function and calculate metrics for best topologies found
            else:

                print("Evaluating... ", end='\r', flush=True)

                # Run random search optimisation method
                if optimisation_method == "random":

                    cost_function_optimisation_algorithms.random_search(constellation_name, num_snapshot, 50, total_sat,
                                                                        satellite_degree_constraints, dcmst, location)

                # Run evolutionary strategy optimisation method
                elif optimisation_method == "evolutionary":

                    cost_function_optimisation_algorithms.evolutionary_search(constellation_name, num_snapshot,
                                                                              total_sat, satellite_degree_constraints,
                                                                              dcmst, location)

                # Run machine learning-based optimisation method
                elif optimisation_method == "machine_learning":

                    cost_function_optimisation_algorithms.machine_learning_optimisation(constellation_name, location,
                                                                                        num_snapshot, total_sat,
                                                                                        satellite_degree_constraints,
                                                                                        dcmst)

                else:
                    raise ValueError("That cost function optimisation method does not exist.")

                print("Novel Topology Algorithm Cost Function Optimisation Completed")

        else:
            raise ValueError("That topology design method does not exist.")

print("\nExecution Time: " + str(time.time() - start) + "s \n")

# References
# Argparse Terminology - https://stackoverflow.com/questions/19124304/what-does-metavar-and-action-mean-in-argparse-in-
# python
# Argsort Descending Order - https://stackoverflow.com/questions/16486252/is-it-possible-to-use-argsort-in-descending-
# order
# Argsort Error - https://stackoverflow.com/questions/53923914/weird-wrong-outpout-of-np-argsort
# Asserting Numpy Equality - https://stackoverflow.com/questions/3302949/best-way-to-assert-for-numpy-array-equality
# Building Project Documentation - https://medium.com/@pratikdomadiya123/build-project-documentation-quickly-with-the-
# sphinx-python-2a9732b66594
# cdist - https://vaghefi.medium.com/fast-distance-calculation-in-python-bb2bc9810ea5
# Check Number of Arguments - https://stackoverflow.com/questions/10698468/how-to-check-if-any-arguments-have-been-
# passed-to-argparse
# Combinations of Arrays - https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations
# -of-two-arrays
# Common Elements of Numpy Array - https://stackoverflow.com/questions/44265572/find-common-elements-in-2d-numpy-arrays
# Conditional List Comprehension - https://stackoverflow.com/questions/4260280/if-else-in-a-list-comprehension
# Conversion with Numpy - https://stackoverflow.com/questions/17973507/why-is-converting-a-long-2d-list-to-numpy-array-
# so-slow
# Conversion with Numpy 2 - https://stackoverflow.com/questions/8466014/how-to-convert-a-python-set-to-a-numpy-array
# Critical Sections - https://stackoverflow.com/questions/29798635/python-locking-critical-section
# Current Directory - https://stackoverflow.com/questions/5137497/find-the-current-directory-and-files-directory
# DCMST Primal Algorithm - https://www-sciencedirect-com.ezphost.dur.ac.uk/science/article/pii/0305054880900222?via%3Di
# hub
# Deleting Rows - https://stackoverflow.com/questions/3877491/deleting-rows-in-numpy-array
# Docstrings - https://realpython.com/documenting-python-code/
# Docstrings Tutorial - https://www.datacamp.com/tutorial/docstrings-python
# Earth Radius - https://en.wikipedia.org/wiki/Earth_radius
# Empty Numpy Array - https://stackoverflow.com/questions/11295609/how-can-i-check-whether-a-numpy-array-is-empty-or-not
# Fast Calculation of Euclidean Distance - https://vaghefi.medium.com/fast-distance-calculation-in-python-bb2bc9810ea5
# Feature Scaling - https://en.wikipedia.org/wiki/Feature_scaling
# File Existence - https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists-without-exceptions
# File Handling - https://www.w3schools.com/python/python_file_handling.asp
# First Occurrence of Value - https://stackoverflow.com/questions/16243955/numpy-first-occurrence-of-value-greater-than-
# existing-value
# Float vs Integer Comparisons - https://stackoverflow.com/questions/30100725/why-are-some-float-integer-comparisons-
# four-times-slower-than-others
# Hypatia - https://github.com/snkas/hypatia/tree/master
# Hyperparameter Tuning Introduction - https://aws.amazon.com/what-is/hyperparameter-tuning/
# Importing - https://stackoverflow.com/questions/4383571/importing-files-from-different-folder
# Import Error - https://stackoverflow.com/questions/338768/python-error-importerror-no-module-named
# __init__ Files - https://stackoverflow.com/questions/448271/what-is-init-py-for
# Internal Imports - https://www.reddit.com/r/learnpython/comments/yusnr0/python_modulespackages_internal_imports/
# Intersecting Arrays - https://stackoverflow.com/questions/25220975/find-the-non-intersecting-values-of-two-arrays
# Iterating over Numpy Columns - https://stackoverflow.com/questions/10148818/numpy-how-to-iterate-over-columns-of-array
# Itertools and Numpy Product - https://stackoverflow.com/questions/28684492/numpy-equivalent-of-itertools-product
# Loop Speed-Up - https://medium.com/@nirmalya.ghosh/13-ways-to-speedup-python-loops-e3ee56cd6b73
# Main Function Call - https://stackoverflow.com/questions/42012140/if-name-main-function-call
# __main__ Use - https://stackoverflow.com/questions/4042905/what-is-main-py
# Module Imports - https://stackoverflow.com/questions/20075884/python-import-module-from-another-directory-at-the-same-
# level-in-project-hierar
# Multiple Indices - https://stackoverflow.com/questions/14162026/how-to-get-the-values-from-a-numpy-array-using-
# multiple-indices
# Multiprocessing Errors - https://www.reddit.com/r/learnpython/comments/yo8pab/why_does_multiprocessing_not_work_
# properly_in/
# Multiprocessing Imap - https://stackoverflow.com/questions/19063238/in-what-situation-do-we-need-to-use-
# multiprocessing-pool-imap-unordered
# Multiplying Numpy Array with Scalar - https://stackoverflow.com/questions/53485221/numpy-multiply-array-with-scalar
# Multiprocessing Pool Errors - https://stackoverflow.com/questions/41169146/python-multiprocessing-pool-map-seems-to-
# not-call-function-at-all
# Multiprocessing - https://stackoverflow.com/questions/35675114/multiprocessing-in-python-going-from-a-for-loop-to-
# multiprocessing-for-loop?rq=3
# Nested For Loop List Comprehension - https://stackoverflow.com/questions/3633140/nested-for-loops-using-list-
# comprehension
# Networkx Documentation - https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.
# shortest_paths.unweighted.all_pairs_shortest_path.html
# Numpy Documentation - https://numpy.org/doc/2.2/reference/index.html
# Numpy Equality - https://stackoverflow.com/questions/10580676/comparing-two-numpy-arrays-for-equality-element-wise
# Numpy Linalg Overhead - https://stackoverflow.com/questions/49866638/why-is-numpy-linalg-norm-slow-when-called-many-
# times-for-small-size-data
# Numpy savetxt Format - https://stackoverflow.com/questions/71499463/how-to-export-numpy-array-without-brackets
# Orbital Distance - https://space.stackexchange.com/questions/27872/how-to-calculate-the-orbital-distance-between-2-
# satellites-given-the-tles
# Or Order - https://stackoverflow.com/questions/55503078/why-does-the-order-of-statements-in-an-or-condition-matter
# OS Library Commands - https://stackoverflow.com/questions/8933237/how-do-i-check-if-a-directory-exists-in-python
# Parsing Lists - https://stackoverflow.com/questions/15753701/how-can-i-pass-a-list-as-a-command-line-argument-with-
# argparse
# Parsing Integer Lists - https://stackoverflow.com/questions/15459997/passing-integer-lists-to-python
# Prim's Algorithm - https://en.wikipedia.org/wiki/Prim%27s_algorithm
# Pyephem Code - https://github.com/brandon-rhodes/pyephem
# Pyephem Documentation - https://rhodesmill.org/pyephem/quick
# Python Documentation - https://docs.python.org/3/
# Reading Files - https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
# Reading Files Line By Line - https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
# Relative Imports - https://stackoverflow.com/questions/16981921/relative-imports-in-python-3
# Satellite Test Data - https://github.com/AtlantaEmrys2002/hypatia/tree/master/paper/satellite_networks_state
# Saving Numpy Data to Files - https://stackoverflow.com/questions/28439701/how-to-save-and-load-numpy-array-data-
# properly
# Scaling vs Normalisation - https://www.kaggle.com/code/alexisbcook/scaling-and-normalization
# SciPy Documentation - https://docs.scipy.org/doc/scipy/reference/index.html
# Selecting Columns - https://stackoverflow.com/questions/22927181/selecting-specific-rows-and-columns-from-numpy-array
# Selecting Numpy Rows - https://stackoverflow.com/questions/58079075/numpy-select-rows-based-on-condition
# Selecting Rows - https://stackoverflow.com/questions/73733370/select-rows-of-numpy-array-based-on-column-values
# Selecting Rows 2 - https://stackoverflow.com/questions/56419519/selecting-rows-of-numpy-array-that-contains-all-the-
# given-values
# Selecting Rows 3 - https://stackoverflow.com/questions/23359886/selecting-rows-in-numpy-ndarray-based-on-the-value-of-
# two-columns
# Selecting Specific Numpy Rows - https://stackoverflow.com/questions/22927181/selecting-specific-rows-and-columns-from-
# numpy-array
# Set vs List Search - https://stackoverflow.com/questions/5993621/fastest-way-to-search-a-list-in-python
# SkyField Documentation - https://rhodesmill.org/skyfield/ & https://rhodesmill.org/skyfield/toc.html
# Skyfield Timings - https://stackoverflow.com/questions/49494082/skyfield-achieve-sgp4-results-with-1-second-
# periodicity-for-given-time-interval
# Speedup For Loops - https://stackoverflow.com/questions/8770791/how-to-speed-up-python-loop
# Spg4 Import Error - https://stackoverflow.com/questions/72639365/no-module-named-sgp4-but-requirement-already-
# satisfied-sgp4-in-c
# Sphinx - https://www.sphinx-doc.org/en/master/tutorial/automatic-doc-generation.html
# Sphinx Tutorial - https://medium.com/@pratikdomadiya123/build-project-documentation-quickly-with-the-sphinx-python-
# 2a9732b66594
# Sphinx Tutorial - https://stackoverflow.com/questions/74787850/generating-documentation-for-multiple-folders-with-
# sphinx
# Sorting by Column - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
# Subset of Rows - https://stackoverflow.com/questions/70341932/select-subset-of-rows-of-numpy-array-based-on-a-
# selection-of-rows-in-another-arr
# Sys Arguments - https://stackoverflow.com/questions/11853508/how-to-make-sys-argv-arguments-optional
# Take (Numpy) - https://stackoverflow.com/questions/60192041/python-numpy-np-take-with-2-dimensional-array
# Threshold Values with Numpy - https://stackoverflow.com/questions/37973135/numpy-argmin-for-elements-greater-than-a-
# threshold
# Threshold Zero - https://stackoverflow.com/questions/28430904/set-numpy-array-elements-to-zero-if-they-are-above-a-
# specific-threshold
# TLE Definitions - https://platform-cdn.leolabs.space/static/files/tle_definition.pdf?7ba94f05897b4ae630a3c5b65be739
# 6c642d9c72
# Type Hinting - https://stackoverflow.com/questions/35673895/type-hinting-annotation-pep-484-for-numpy-ndarray
# Type Hinting with Numpy - https://stackoverflow.com/questions/58877370/pycharm-typehint-for-numpy-array-of-objects
# Type Hinting with Numpy - https://stackoverflow.com/questions/70714087/how-to-typehint-numpy-array-of-floats-and-or-
# integers
# Typing - https://stackoverflow.com/questions/14113187/how-do-you-set-a-conditional-in-python-based-on-datatypes
# Typing and Int Assignment - https://stackoverflow.com/questions/21390612/typeerror-int-object-does-not-support-item-
# assignment
# Unique Values from List - https://stackoverflow.com/questions/12897374/get-unique-values-from-a-list-in-python
# World Geodetic System - https://en.wikipedia.org/wiki/World_Geodetic_System#Definition
# Writing Packages - https://stackoverflow.com/questions/15746675/how-to-write-a-python-module-package
# Zero-Size Reduction - https://www.reddit.com/r/learnpython/comments/lhljgh/valueerror_zerosize_array_to_reduction_
# operation/

# Extra Code Used During Analysis, Optimisation, and Debugging

# Used to check that visibility matrices change over time (confirming theory)
# for k in range(1, num_snapshot):
#     previous = np.load("./visibility_matrices/visibility_matrix_" + str(k) + ".npy")
#     now = np.load("./visibility_matrices/visibility_matrix_" + str(k - 1) + ".npy")
#     if np.array_equal(previous, now) is False:
#         print("yes")
