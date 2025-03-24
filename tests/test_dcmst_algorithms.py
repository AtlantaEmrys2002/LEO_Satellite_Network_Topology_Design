# Libraries
import dcmst_construction_algorithms
import numpy as np
import unittest


class TestDegreeConstrainedMinimumSpanningTreeConstructionAlgorithms(unittest.TestCase):

    def test_modified_prims(self):
        """
        Checks the modified version of Prim's Algorithm (that forms the basis of the Primal DCMST construction
        algorithm) returns a correct DCMST. This example cost matrix confirms that a correct DCMST is constructed (node
         2 would have a degree > 3 if this were a MST construction algorithm, so this example ensures all constraints
         are satisfied).
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

        self.assertTrue(np.array_equal(
            dcmst_construction_algorithms.modified_prims_algorithm(test_costs, np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), 9,
                                                                   0)[0], answer))

    def test_edge_exchange(self):
        # Tests whether function that constructs two subtrees following edge deletion is correct
        # self.assertTrue(dcmst_construction_algorithms.subtree_builder(9, np.asarray(
        #     [[0, 1], [1, 2], [2, 3], [2, 8], [3, 4], [4, 5], [5, 6], [6, 7]]), 2) == ({8, 1, 2, 0}, {3, 4, 5, 6, 7}))

        self.assertTrue(dcmst_construction_algorithms.subtree_builder(9, np.asarray(
            [[0, 1], [1, 2], [2, 3], [2, 8], [3, 4], [4, 5], [5, 6], [6, 7]])) == ({8, 1, 2, 0}, {3, 4, 5, 6, 7}))

        answer = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0],
                           [1, 0, 1, 0, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 1, 0, 0, 1],
                           [0, 0, 0, 0, 1, 0, 0, 0, 0],
                           [0, 0, 0, 1, 0, 1, 0, 0, 0],
                           [0, 0, 1, 0, 1, 0, 1, 0, 0],
                           [0, 0, 0, 0, 0, 1, 0, 1, 0],
                           [0, 0, 0, 0, 0, 0, 1, 0, 0],
                           [0, 0, 1, 0, 0, 0, 0, 0, 0]])

        # Passing tree such that edge [0, 7] could be replaced with [0, 1] to create a better degree constrained minimum
        # spanning tree
        self.assertTrue(
            np.array_equal(dcmst_construction_algorithms.edge_exchange(np.asarray([[-1, 4, -1, -1, -1, -1, -1, 8, -1],
                                                                                   [4, -1, 8, -1, -1, -1, -1, 11, -1],
                                                                                   [-1, 8, -1, 7, -1, 4, -1, -1, 2],
                                                                                   [-1, -1, 7, -1, 9, 14, -1, -1, -1],
                                                                                   [-1, -1, -1, 9, -1, 10, -1, -1, -1],
                                                                                   [-1, -1, 4, 14, 10, -1, 2, -1, -1],
                                                                                   [-1, -1, -1, -1, -1, 2, -1, 1, 6],
                                                                                   [8, 11, -1, -1, -1, -1, 1, -1, 7],
                                                                                   [-1, -1, 2, -1, -1, -1, 6, 7, -1]]),
                                                                       np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), 9,
                                                                       np.array(
                                                                           [[0, 0, 0, 0, 0, 0, 0, 1, 0],
                                                                            [0, 0, 1, 0, 0, 0, 0, 0, 0],
                                                                            [0, 1, 0, 0, 0, 1, 0, 0, 1],
                                                                            [0, 0, 0, 0, 1, 0, 0, 0, 0],
                                                                            [0, 0, 0, 1, 0, 1, 0, 0, 0],
                                                                            [0, 0, 1, 0, 1, 0, 1, 0, 0],
                                                                            [0, 0, 0, 0, 0, 1, 0, 1, 0],
                                                                            [1, 0, 0, 0, 0, 0, 1, 0, 0],
                                                                            [0, 0, 1, 0, 0, 0, 0, 0, 0]]),
                                                                       np.array([1, 2, 3, 1, 2, 3, 2, 1, 1]))[0],
                           answer))

        # Passing tree such that edge [2, 3] could be replaced with [2, 5] to create a better degree constrained minimum
        # spanning tree
        self.assertTrue(
            np.array_equal(dcmst_construction_algorithms.edge_exchange(np.asarray([[-1, 4, -1, -1, -1, -1, -1, 8, -1],
                                                                                   [4, -1, 8, -1, -1, -1, -1, 11, -1],
                                                                                   [-1, 8, -1, 7, -1, 4, -1, -1, 2],
                                                                                   [-1, -1, 7, -1, 9, 14, -1, -1, -1],
                                                                                   [-1, -1, -1, 9, -1, 10, -1, -1, -1],
                                                                                   [-1, -1, 4, 14, 10, -1, 2, -1, -1],
                                                                                   [-1, -1, -1, -1, -1, 2, -1, 1, 6],
                                                                                   [8, 11, -1, -1, -1, -1, 1, -1, 7],
                                                                                   [-1, -1, 2, -1, -1, -1, 6, 7, -1]]),
                                                                       np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), 9,
                                                                       np.array(
                                                                           [[0, 1, 0, 0, 0, 0, 0, 0, 0],
                                                                            [1, 0, 1, 0, 0, 0, 0, 0, 0],
                                                                            [0, 1, 0, 1, 0, 0, 0, 0, 1],
                                                                            [0, 0, 1, 0, 1, 0, 0, 0, 0],
                                                                            [0, 0, 0, 1, 0, 1, 0, 0, 0],
                                                                            [0, 0, 0, 0, 1, 0, 1, 0, 0],
                                                                            [0, 0, 0, 0, 0, 1, 0, 1, 0],
                                                                            [0, 0, 0, 0, 0, 0, 1, 0, 0],
                                                                            [0, 0, 1, 0, 0, 0, 0, 0, 0]]),
                                                                       np.array([1, 2, 3, 1, 2, 3, 2, 1, 1]))[0],
                           answer))

    def test_increase_connectivity(self):
        """
        Checks the ISLs are added correctly such that lower cost edges are added first before higher cost edges.
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

        current_topology, current_degree = dcmst_construction_algorithms.modified_prims_algorithm(test_costs, np.array(
            [3, 3, 3, 3, 3, 3, 3, 3, 3]), 9, 0)

        answer = np.array([[0, 1, 0, 0, 0, 0, 0, 1, 0],
                           [1, 0, 1, 0, 0, 0, 0, 1, 0],
                           [0, 1, 0, 0, 0, 1, 0, 0, 1],
                           [0, 0, 0, 0, 1, 0, 0, 0, 0],
                           [0, 0, 0, 1, 0, 1, 0, 0, 0],
                           [0, 0, 1, 0, 1, 0, 1, 0, 0],
                           [0, 0, 0, 0, 0, 1, 0, 1, 1],
                           [1, 1, 0, 0, 0, 0, 1, 0, 0],
                           [0, 0, 1, 0, 0, 0, 1, 0, 0]])

        self.assertTrue(np.array_equal(
            dcmst_construction_algorithms.increase_connectivity(current_topology, np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]),
                                                                current_degree, test_costs, 9), answer))

    def test_prufer_encoding(self):
        """
        Tests whether trees are correctly encoded by Prufer encoding function used by the GA that constructs DCMSTs.
        """
        # Example tree taken from original paper on comparison of DCMST construction algorithms (see report for full
        # citation)

        tree = np.array([[0, 0, 0, 1, 0, 0, 1, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [1, 0, 0, 0, 1, 0, 0, 0, 1],
                         [1, 0, 1, 0, 0, 1, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 1, 0, 0]])

        prufer_encoding = dcmst_construction_algorithms.prufer_encode(tree)

        self.assertTrue(np.array_equal(prufer_encoding, np.array([8, 7, 0, 6, 7, 0, 6])))

    def test_prufer_decoding(self):
        """
        Tests whether trees are correctly decoded from Prufer decoding function to adjacency matrix format used by the
        GA that constructs DCMSTs.
        """
        # Example tree taken from original paper on comparison of DCMST construction algorithms (see report for full
        # citation)

        tree = np.array([[0, 0, 0, 1, 0, 0, 1, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [1, 0, 0, 0, 1, 0, 0, 0, 1],
                         [1, 0, 1, 0, 0, 1, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 1, 0, 0]])

        prufer_decoding, _ = dcmst_construction_algorithms.prufer_decode(np.array([8, 7, 0, 6, 7, 0, 6]))

        self.assertTrue(np.array_equal(prufer_decoding, tree))

    def test_prufer_encoding_and_decoding(self):
        """
        Tests whether Prufer encoding/decoding functions used by GA that constructs DCMSTs interact properly.
        """

        # Prufer Number of Tree
        prufer_encoding = np.array([1, 3, 5, 2, 1, 0, 4], dtype=np.int32)

        # Tree as Adjacency Matrix
        tree = np.array([[0, 1, 0, 0, 1, 0, 0, 0, 0],
                         [1, 0, 1, 0, 0, 0, 1, 0, 0],
                         [0, 1, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 1, 0],
                         [1, 0, 0, 0, 0, 0, 0, 0, 1],
                         [0, 0, 1, 1, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0]])

        self.assertTrue(np.array_equal(dcmst_construction_algorithms.prufer_decode(
            dcmst_construction_algorithms.prufer_encode(tree))[0], tree))

        self.assertTrue(np.array_equal(dcmst_construction_algorithms.prufer_encode(
            dcmst_construction_algorithms.prufer_decode(prufer_encoding)[0]), prufer_encoding))

    def test_degree_checker(self):
        """
        Tests whether function that ensures tree encoded as Prufer number correctly determines if degree constraints
        have been violated when DCMST constructed.
        """
        # Prufer Number of Tree
        prufer_encoding = np.array([1, 3, 5, 2, 1, 0, 4], dtype=np.int32)

        self.assertTrue(dcmst_construction_algorithms.check_degree(prufer_encoding, np.array([3 for _ in range(9)]),
                                                                   9))

        self.assertFalse(dcmst_construction_algorithms.check_degree(prufer_encoding, np.array([2 for _ in range(9)]),
                                                                    9))


if __name__ == '__main__':
    unittest.main()

# References
# Example Graph - https://www.geeksforgeeks.org/prims-minimum-spanning-tree-mst-greedy-algo-5/
