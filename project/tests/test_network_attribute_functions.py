import numpy as np
import os
import project.satellite_network_attribute_functions as satnet
import unittest


class TestNetworkAttributeFunctions(unittest.TestCase):
    def test_visibility_function(self):
        """
        Checks visibility function (used to calculate whether each pair of satellites within the network are visible to
        one another).
        """
        self.assertTrue(np.array_equal(satnet.visibility_function(np.asarray([[0.0, 8.0, 3.0, 4.0],
                                                                              [8.0, 0.0, 5.0, 6.0],
                                                                              [3.0, 5.0, 0.0, 7.0],
                                                                              [4.0, 6.0, 7.0, 0.0]], dtype=np.float64),
                                                                  4), np.asarray(
            [[False, False, True, True], [False, False, False, False], [True, False, False, False],
             [True, False, False, False]])))

    def test_time_visibility_function(self):
        """
        Checks time visibility function (used to calculate how long each pair of satellites is visible to one another
        from current snapshot) - both basic and advanced (time extends beyond defined orbit and "wraps around")
        """
        # Create test data
        if os.path.isdir("./test_data/visibility_matrices") is False:

            # Create directory in which to store distance matrices
            try:
                # os.mkdir("./visibility_matrices")
                os.makedirs("./test_data/visibility_matrices")
            except OSError:
                print("Directory to store visibility matrices could not be created.")

            # Create test visibility matrices
            visibility_matrix_1 = np.asarray([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
            visibility_matrix_2 = np.asarray([[0, 0, 1], [0, 0, 1], [1, 1, 0]])
            visibility_matrix_3 = np.asarray([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
            visibility_matrix_4 = np.asarray([[0, 0, 1], [0, 0, 0], [1, 0, 0]])

            np.save("./test_data/visibility_matrices/visibility_matrix_0.npy", visibility_matrix_1)
            np.save("./test_data/visibility_matrices/visibility_matrix_1.npy", visibility_matrix_2)
            np.save("./test_data/visibility_matrices/visibility_matrix_2.npy", visibility_matrix_3)
            np.save("./test_data/visibility_matrices/visibility_matrix_3.npy", visibility_matrix_4)

        # Basic Test - from snapshot zero
        self.assertTrue(np.array_equal(satnet.time_visibility_function(4, 3, 0, "test_data"),
                                       np.asarray([[0, 0, 3], [0, 0, 2], [3, 2, 0]])))

        # Advanced Test - from snapshot other than zero (current_id will have to be reset at some point to zero)
        self.assertTrue(np.array_equal(satnet.time_visibility_function(4, 3, 2, "test_data"),
                                       np.asarray([[0, 0, 3], [0, 0, 0], [3, 0, 0]])))

    def test_cost_function(self):
        """
        Tests cost function returns correct results given a series of matrices representing the different network
        attributes taken into consideration.
        """
        dist = np.array([[-1.0, 2.0, 9.0], [2.0, -1.0, 10.0], [9.0, 10.0, -1.0]])
        time_visibility = np.array([[0, 0, 2], [0, 0, 4], [2, 4, 0]])
        sunlight = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
        visibility = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])

        # Use all close here as using floating point representation
        self.assertTrue(np.allclose(satnet.cost_function(visibility, time_visibility, dist, sunlight, alpha=0.2,
                                                         beta=0.3, gamma=0.5, total_satellites=3),
                                    np.asarray([[-1., 0.2, 0.32916667], [0.2, -1., 0.84], [0.32916667, 0.84, -1.]])))


if __name__ == '__main__':
    unittest.main()
