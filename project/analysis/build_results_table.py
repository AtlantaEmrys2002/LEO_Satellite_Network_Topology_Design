# Libraries
import csv
import os


def build_results_overview():
    """
    Used to build an overview of the results of report.sh (topology builds of report.sh and their respective metrics).
    """
    # Stores best values for each combination
    results = []

    # for each algorithm
    for algorithm in ['plus_grid', 'mdtd', 'novel']:

        # benchmarks
        if algorithm != 'novel':

            # for each satellite mega-constellation for which a topology is constructed
            for constellation in ['iridium-780', 'kuiper-630', 'starlink-550', 'telesat-1015']:

                # location of results
                location = './Results/' + algorithm + '/' + constellation + '/results.csv'

                # check if results file exists
                if os.path.isfile(location):

                    # read results
                    with open(location) as csvfile:
                        reader = csv.DictReader(csvfile)
                        for row in reader:
                            optimisation_method = float("nan")
                            dcmst = float("nan")
                            alpha = float("nan")
                            beta = float("nan")
                            gamma = float("nan")

                            # read values
                            values = dict(algorithm=algorithm, optimisation_method=optimisation_method, dcmst=dcmst,
                                          constellation=constellation, alpha=alpha, beta=beta, gamma=gamma,
                                          mean_latency=row['mean_latency'], average_hop_count=row["average_hop_count"],
                                          link_churn=row["link_churn"])

                            results.append(values)

        else:

            # for each method used to optimise cost function
            for optimisation_method in ['evolutionary', 'random']:

                # for each method used to construct a DCMST
                for dcmst in ['aco', 'ga', 'primal']:

                    # for each satellite mega-constellation for which a topology is constructed
                    for constellation in ['iridium-780', 'kuiper-630', 'starlink-550', 'telesat-1015']:

                        # location of results
                        location = ('./Results/' + algorithm + '/' + optimisation_method + '/' + dcmst + '/' +
                                    constellation + '/results.csv')

                        # check if location exists
                        if os.path.isfile(location):

                            # temporary results storer
                            best_result = []

                            # read results
                            with open(location) as csvfile:
                                reader = csv.DictReader(csvfile)
                                for row in reader:
                                    # read values
                                    values = dict(algorithm=algorithm, optimisation_method=optimisation_method,
                                                  dcmst=dcmst, constellation=constellation, alpha=row["alpha"],
                                                  beta=row["beta"], gamma=row["gamma"],
                                                  mean_latency=row["mean_latency"],
                                                  average_hop_count=row["average_hop_count"],
                                                  link_churn=row["link_churn"])

                                    # find the best values for each metric
                                    if not best_result:
                                        # best_result = values
                                        best_result = [values, values, values]
                                    else:
                                        if best_result[0]["mean_latency"] > values["mean_latency"]:
                                            # best_result = values
                                            best_result[0] = values
                                        elif best_result[1]["average_hop_count"] > values["average_hop_count"]:
                                            # best_result = values
                                            best_result[1] = values
                                        elif best_result[2]["link_churn"] > values["link_churn"]:
                                            # best_result = values
                                            best_result[2] = values
                            for k in range(3):
                                # results.append(best_result)
                                results.append(best_result[k])

    # Write results to overview table
    with open('./results_overview.csv', 'w', newline='') as csvfile:
        fieldnames = ["algorithm", "optimisation_method", "dcmst", "constellation", "alpha", "beta", "gamma",
                      "mean_latency", "average_hop_count", "link_churn"]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)


build_results_overview()

# References:
# File Exists - https://stackoverflow.com/questions/82831/how-do-i-check-whether-a-file-exists-without-exceptions
