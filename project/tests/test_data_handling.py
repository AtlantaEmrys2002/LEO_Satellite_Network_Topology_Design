# Libraries
import project.data_handling

import numpy as np
import os
import unittest


class TestDataHandlingFunctions(unittest.TestCase):
    def test_write_topology_to_file_function(self):
        """
        Checks that topologies are correctly saved and loaded to text file compatible with Hypatia software.
        """
        # Example topology to be saved
        test_topology = np.array([[0, 1, 0, 1], [1, 0, 1, 1], [0, 1, 0, 0], [1, 1, 0, 0]])

        project.data_handling.write_topology_to_file('test_topology_save.txt', test_topology,
                                                     'novel')

        # Check file exists
        self.assertTrue(os.path.isfile('./novel/isl_topologies/test_topology_save.txt'))

        # Read in topology file
        topology_isls = np.loadtxt('./novel/isl_topologies/test_topology_save.txt').astype(int).T

        # Recreate topology matrix for ISLs
        topology_matrix = np.zeros((4, 4))
        topology_matrix[topology_isls[0], topology_isls[1]] = 1
        topology_matrix[topology_isls[1], topology_isls[0]] = 1

        # Check topology creation correct
        self.assertTrue(np.array_equal(topology_matrix, test_topology))


if __name__ == '__main__':
    unittest.main()
