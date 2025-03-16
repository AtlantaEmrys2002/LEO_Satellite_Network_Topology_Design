# Libraries
import dcmst_construction_algorithms
import numpy as np
import unittest

class TestDegreeConstrainedMinimumSpanningTreeConstructionAlgorithms(unittest.TestCase):

    def test_modified_prims(self):
        """
        Checks the modified version of Prim's Algorithm (that forms the basis of the Primal DCMST construction algorithm) returns a correct DCMST. This example cost matrix confirms that a correct DCMST is constructed (node 2 would have a degree > 3 if this were a MST construction algorithm, so this example ensures all constraints are satisfied).
        """
        # Test Cost Matrix
        test_costs = np.asarray([[-1, 4, -1, -1, -1, -1, -1, 8, -1],
                                 [4, -1, 8, -1, -1, -1, -1, 11, -1],
                                 [-1, 8, -1, 7, -1, 4, -1, -1, 2],
                                 [-1, -1, 7, -1, 9, 14, -1, -1, -1],
                                 [-1, -1, -1, 9, -1, 10, -1, -1, -1],
                                 [-1, -1, 4, 14, 10, -1, 2, -1, -1],
                                 [-1, -1, -1, -1, -1, 2, -1, 1, 6],
                                 [8, 11, -1, -1, -1, -1, 1, -1, 7],
                                 [-1, -1, 2, -1, -1, -1, 6, 7, -1]])

        # Algorithm Answer (one of many possible solutions - several edges have same length)
        answer = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0],
                  [1, 0, 1, 0, 0, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0, 1, 0, 0, 1],
                  [0, 0, 0, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 1, 0, 1, 0, 0, 0],
                  [0, 0, 1, 0, 1, 0, 1, 0, 0],
                  [0, 0, 0, 0, 0, 1, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 1, 0, 0],
                  [0, 0, 1, 0, 0, 0, 0, 0, 0]])

        self.assertTrue(np.array_equal(dcmst_construction_algorithms.modified_prims_algorithm(test_costs, np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), 9, 0)[0], answer))

    def test_increase_connectivity(self):
        """
        Checks the ISLs are added correctly such that lower cost edges are added first before
        """
        # Test Cost Matrix
        test_costs = np.asarray([[-1, 4, -1, -1, -1, -1, -1, 8, -1],
                                 [4, -1, 8, -1, -1, -1, -1, 11, -1],
                                 [-1, 8, -1, 7, -1, 4, -1, -1, 2],
                                 [-1, -1, 7, -1, 9, 14, -1, -1, -1],
                                 [-1, -1, -1, 9, -1, 10, -1, -1, -1],
                                 [-1, -1, 4, 14, 10, -1, 2, -1, -1],
                                 [-1, -1, -1, -1, -1, 2, -1, 1, 6],
                                 [8, 11, -1, -1, -1, -1, 1, -1, 7],
                                 [-1, -1, 2, -1, -1, -1, 6, 7, -1]])


        current_topology, current_degree = dcmst_construction_algorithms.modified_prims_algorithm(test_costs, np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), 9, 0)

        answer = np.array([[0, 1, 0, 0, 0, 0, 0, 1, 0],
                           [1, 0, 1, 0, 0, 0, 0, 1, 0],
                           [0, 1, 0, 0, 0, 1, 0, 0, 1],
                           [0, 0, 0, 0, 1, 0, 0, 0, 0],
                           [0, 0, 0, 1, 0, 1, 0, 0, 0],
                           [0, 0, 1, 0, 1, 0, 1, 0, 0],
                           [0, 0, 0, 0, 0, 1, 0, 1, 1],
                           [1, 1, 0, 0, 0, 0, 1, 0, 0],
                           [0, 0, 1, 0, 0, 0, 1, 0, 0]])

        self.assertTrue(np.array_equal(dcmst_construction_algorithms.increase_connectivity(current_topology, np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), current_degree, test_costs, 9), answer))


if __name__ == '__main__':
    unittest.main()

# References
# Example Graph - https://www.geeksforgeeks.org/prims-minimum-spanning-tree-mst-greedy-algo-5/
