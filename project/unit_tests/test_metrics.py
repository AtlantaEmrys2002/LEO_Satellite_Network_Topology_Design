# Libraries
from astropy.constants import c
from project.metrics import hop_count, link_churn, propagation_delay
import numpy as np
import os
import unittest


class TestMetricFunctions(unittest.TestCase):

    def test_av_hop_count(self):
        """
        Ensures the function that measures the average hop count of a given satellite network with a given ISL
        topology is correct.
        """

        example_topology = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0],
                                     [1, 0, 1, 0, 0, 0, 0, 0, 0],
                                     [0, 1, 0, 0, 0, 1, 0, 0, 1],
                                     [0, 0, 0, 0, 1, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 0, 1, 0, 0, 0],
                                     [0, 0, 1, 0, 1, 0, 1, 0, 0],
                                     [0, 0, 0, 0, 0, 1, 0, 1, 0],
                                     [0, 0, 0, 0, 0, 0, 1, 0, 0],
                                     [0, 0, 1, 0, 0, 0, 0, 0, 0]])

        example_dist_matrix = np.array([[0, 4.0, 0, 0, 0, 0, 0, 8.0, 0],
                                        [4.0, 0, 8.0, 0, 0, 0, 0, 11.0, 0],
                                        [0, 8.0, 0, 7.0, 0, 4.0, 0, 0, 2.0],
                                        [0, 0, 7.0, 0, 9.0, 14.0, 0, 0, 0],
                                        [0, 0, 0, 9.0, 0, 10.0, 0, 0, 0],
                                        [0, 0, 4.0, 14.0, 10.0, 0, 2.0, 0, 0],
                                        [0, 0, 0, 0, 0, 2.0, 0, 1.0, 6.0],
                                        [8.0, 11.0, 0, 0, 0, 0, 1.0, 0, 7.0],
                                        [0, 0, 2.0, 0, 0, 0, 6.0, 7.0, 0],])

        answer = 94 / 36

        self.assertTrue(hop_count(example_topology, example_dist_matrix, 9) == answer)

    def test_link_churn(self):
        """
        Ensures the function that measures the link churn of a satellite network with a given ISL topology is correct.
        """

        # Create test topologies if they have not already been created
        if os.path.isdir("./test_data/isl_topologies/test_constellation") is False:

            # Create directory in which to store test topologies
            try:
                os.makedirs("./test_data/isl_topologies/test_constellation")
            except OSError:
                print("Directory to store ISL topologies (for testing) could not be created.")

            # Create test topologies
            snapshot_0 = np.array([[0, 1, 0, 1, 0],
                                   [1, 0, 1, 0, 0],
                                   [0, 1, 0, 0, 0],
                                   [1, 0, 0, 0, 1],
                                   [0, 0, 0, 1, 0],])

            snapshot_1 = np.array([[0, 1, 0, 1, 0],
                                   [1, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0],
                                   [1, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0],])

            snapshot_2 = np.array([[0, 1, 0, 1, 0],
                                   [1, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0],
                                   [1, 0, 0, 0, 1],
                                   [0, 0, 0, 1, 0],])

            snapshot_3 = np.array([[0, 1, 0, 0, 0],
                                   [1, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0],])

            # Select all ISLs within each topology, sort, and remove duplicates (undirected edges)
            snapshot_0 = np.unique(np.sort(np.argwhere(snapshot_0 > 0)), axis=0)
            snapshot_1 = np.unique(np.sort(np.argwhere(snapshot_1 > 0)), axis=0)
            snapshot_2 = np.unique(np.sort(np.argwhere(snapshot_2 > 0)), axis=0)
            snapshot_3 = np.unique(np.sort(np.argwhere(snapshot_3 > 0)), axis=0)

            # Save test topologies in correct format
            np.savetxt("./test_data/isl_topologies/test_constellation/isls_0.txt", snapshot_0, fmt='%i %i')
            np.savetxt("./test_data/isl_topologies/test_constellation/isls_1.txt", snapshot_1, fmt='%i %i')
            np.savetxt("./test_data/isl_topologies/test_constellation/isls_2.txt", snapshot_2, fmt='%i %i')
            np.savetxt("./test_data/isl_topologies/test_constellation/isls_3.txt", snapshot_3, fmt='%i %i')

        # This the actual link churn of the satellite network over 1 orbit
        actual_link_churn = 5

        # Calculate link churn for topologies stored in test data
        result = link_churn("./test_data/isl_topologies/test_constellation", 4, 5)

        self.assertEqual(result, actual_link_churn)

    def test_propagation_delay(self):
        """
        Ensures the function that measures propagation delay for a given satellite network with a given ISL topology
        is correct.
        """
        example_topology = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0],
                                     [1, 0, 1, 0, 0, 0, 0, 0, 0],
                                     [0, 1, 0, 0, 0, 1, 0, 0, 1],
                                     [0, 0, 0, 0, 1, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 0, 1, 0, 0, 0],
                                     [0, 0, 1, 0, 1, 0, 1, 0, 0],
                                     [0, 0, 0, 0, 0, 1, 0, 1, 0],
                                     [0, 0, 0, 0, 0, 0, 1, 0, 0],
                                     [0, 0, 1, 0, 0, 0, 0, 0, 0]])

        example_dist_matrix = np.array([[0, 4.0, 0, 0, 0, 0, 0, 8.0, 0],
                                        [4.0, 0, 8.0, 0, 0, 0, 0, 11.0, 0],
                                        [0, 8.0, 0, 7.0, 0, 4.0, 0, 0, 2.0],
                                        [0, 0, 7.0, 0, 9.0, 14.0, 0, 0, 0],
                                        [0, 0, 0, 9.0, 0, 10.0, 0, 0, 0],
                                        [0, 0, 4.0, 14.0, 10.0, 0, 2.0, 0, 0],
                                        [0, 0, 0, 0, 0, 2.0, 0, 1.0, 6.0],
                                        [8.0, 11.0, 0, 0, 0, 0, 1.0, 0, 7.0],
                                        [0, 0, 2.0, 0, 0, 0, 6.0, 7.0, 0], ])

        max_result, mean_result = propagation_delay(example_topology, example_dist_matrix, 9)

        mean_actual_answer = (488 / c.to('km/s').value) / 36

        max_actual_answer = 35 / c.to('km/s').value

        # Check max propagation delay returned is correct
        self.assertEqual(max_actual_answer, max_result)

        # Check mean propagation delay returned is correct
        self.assertEqual(mean_actual_answer, mean_result)


if __name__ == '__main__':
    unittest.main()
