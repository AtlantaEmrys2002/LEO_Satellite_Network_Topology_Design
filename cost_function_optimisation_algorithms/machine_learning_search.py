# Libraries
import csv
import os


def machine_learning_optimisation(file_name: str, constellation_name: str, num_snapshots: int, num_sat: int,
                                  degree_constraints: list[int], dcmst_method: str):

    # CREATE LOCATIONS TO STORE RESULTS #

    # Initialise temporary array to store results
    results = []

    # Results location
    location = "./Results/novel/" + constellation_name.lower() + "/" + dcmst_method + "/machine_learning_optimisation/"

    if os.path.isdir(location) is False:

        # Create directory in which to store evolutionary optimisation search results
        try:
            os.makedirs(location)
        except OSError:
            print("Directory to store results of machine learning optimisation could not be created.")

    # MACHINE LEARNING OPTIMISATION #

    # NEED TO IMPLEMENT

    # WRITE RESULTS TO CSV #

    with (open(location + "/results.csv", 'w', newline='')
          as csvfile):
        fieldnames = ['alpha', 'beta', 'gamma', 'mean_latency', 'average_hop_count', 'link_churn']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
