import numpy as np
import satellite_network_attribute_functions as satnet
import unittest


class TestNetworkAttributeFunctions(unittest.TestCase):
    def test_visibility_function(self):
        """
        Checks visibility function (used to calculate whether each pair of satellites within the network are visible to one another).
        """
        self.assertTrue(np.array_equal(satnet.visibility_function(np.asarray([[0.0, 8.0, 3.0, 4.0], [8.0, 0.0, 5.0, 6.0], [3.0, 5.0, 0.0, 7.0], [4.0, 6.0, 7.0, 0.0]], dtype=np.float64), 4), np.asarray([[False, False, True, True], [False, False, False, False], [True, False, False, False], [True, False, False, False]])))


if __name__ == '__main__':
    unittest.main()
