# Libraries
from data_handling import write_optimisation_results_to_csv
import os


def machine_learning_optimisation(constellation_name: str, output_directory: str, num_snapshots: int, num_sat: int,
                                  degree_constraints: list[int], dcmst_method: str):

    # CREATE LOCATIONS TO STORE RESULTS #

    # Initialise temporary array to store results
    results = []

    if os.path.isdir(output_directory) is False:

        # Create directory in which to store evolutionary optimisation search results
        try:
            os.makedirs(output_directory)
        except OSError:
            print("Directory to store results of machine learning optimisation could not be created.")

    # MACHINE LEARNING OPTIMISATION #

    # NEED TO IMPLEMENT

    # WRITE RESULTS TO CSV #
    write_optimisation_results_to_csv(output_directory, "novel", results)
