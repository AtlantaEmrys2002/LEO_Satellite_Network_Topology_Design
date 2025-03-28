# Libraries
import csv


def write_optimisation_results_to_csv(location: str, algorithm: str, results: list):
    """
    Write optimisation/metric results to CSV Format - this code was adapted from documentation -
    https://docs.python.org/3/library/csv.html#csv.DictWriter.
    :param location: directory in which the results of the cost function optimisation/metric evaluation are
     stored
    :param algorithm: type of algorithm with which topologies were constructed (either 'static', 'dynamic', or 'novel')
    :param results: results of topology optimisation and/or topology evaluation:
    """
    with (open(location + '/results.csv', 'w', newline='') as csvfile):

        # Benchmark algorithm results
        if algorithm == "static" or algorithm == "dynamic":
            fieldnames = ['max_latency', 'mean_latency', 'average_hop_count', 'link_churn']
            max_pd, mean_pd, av_hop_count, link_churn = results
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(dict(max_latency=max_pd, mean_latency=mean_pd, average_hop_count=av_hop_count,
                                 link_churn=link_churn))
        # Algorithm is novel (proposed in report)
        elif algorithm == "novel":
            fieldnames = ['alpha', 'beta', 'gamma', 'mean_latency', 'average_hop_count', 'link_churn']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        else:
            raise ValueError("Cannot save optimisation/measurement results for undefined algorithm.")
