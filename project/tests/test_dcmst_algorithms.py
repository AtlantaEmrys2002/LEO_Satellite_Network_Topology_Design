# Libraries
import project.dcmst_construction_algorithms as alg
import numpy as np
import unittest
from scipy.sparse.csgraph import connected_components


class TestDegreeConstrainedMinimumSpanningTreeConstructionAlgorithms(unittest.TestCase):

    # PRIMAL ALGORITHM TESTS #

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
            alg.modified_prims_algorithm(test_costs, np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), 9, 0)[0], answer))

    def test_subtree_builder(self):
        """
        Tests whether function that determines whether exactly two subtrees are created when an edge is deleted from
        the proposed tree, as well as the vertices in each subtree, is correct.
        """

        # Example tree
        tree = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [1, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 1, 0, 0, 0, 0, 1],
                         [0, 0, 1, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0]])

        # Find two subtrees i and j - error would be raised if more than two subtrees were created by edge deletion.
        # i, j = dcmst_construction_algorithms.subtree_builder(tree, np.array([2, 3]))
        i, j = alg.subtree_builder(tree, np.array([2, 3]))

        # Check two subtrees are in correct format and contain the correct vertices
        self.assertTrue((np.array_equal(np.sort(j), np.array([0, 1, 2, 8])) and
                         np.array_equal(np.sort(i), np.array([3, 4, 5, 6, 7]))) or
                        np.array_equal(np.sort(i), np.array([0, 1, 2, 8])) and
                        np.array_equal(np.sort(j), np.array([3, 4, 5, 6, 7])))

    def test_edge_exchange(self):
        """
        Tests whether edges are correctly exchanged (in primal method) to improve total sum cost of DCMST. Example tree
        taken from original DCMST paper by Narula and Ho (see report for full citation).
        """

        cost_matrix = np.array([[0, 2.24, 2.24, 3.61, 6.71, 3.0, 5.39, 8.0, 9.43],
                                [2.24, 0, 2.0, 2.0, 4.47, 2.83, 4.0, 7.28, 7.62],
                                [2.24, 2.0, 0, 4.0, 5.66, 4.47, 6.0, 9.22, 9.49],
                                [3.61, 2.0, 4.0, 0, 4.0, 2.0, 2.0, 5.39, 5.83],
                                [6.71, 4.47, 5.66, 4.0, 0, 6.0, 4.47, 7.81, 5.1],
                                [3.0, 2.83, 4.47, 2.0, 6.0, 0, 2.83, 5.0, 7.07],
                                [5.39, 4.0, 6.0, 2.0, 4.47, 2.83, 0, 3.61, 4.24],
                                [8.0, 7.28, 9.22, 5.39, 7.81, 0, 3.61, 0, 5.0],
                                [9.43, 7.62, 9.49, 5.83, 5.1, 7.07, 0, 5.0, 0], ])

        constraints = np.array([3, 3, 3, 3, 3, 3, 3, 3, 3])

        current_topology = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0],
                                     [1, 0, 1, 1, 0, 0, 0, 0, 0],
                                     [0, 1, 0, 0, 0, 0, 0, 0, 0],
                                     [0, 1, 0, 0, 0, 1, 1, 0, 0],
                                     [0, 0, 0, 0, 0, 0, 0, 0, 1],
                                     [0, 0, 0, 1, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 1, 0, 0, 0, 1, 1],
                                     [0, 0, 0, 0, 0, 0, 1, 0, 0],
                                     [0, 0, 0, 0, 1, 0, 1, 0, 0], ])

        current_degree = np.array([1, 3, 1, 3, 1, 1, 3, 1, 2])

        answer = np.array([[0, 0, 1, 0, 0, 0, 0, 0, 0],
                           [0, 0, 1, 1, 1, 0, 0, 0, 0],
                           [1, 1, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 1, 1, 0, 0],
                           [0, 1, 0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 1, 0, 0, 0, 0, 0],
                           [0, 0, 0, 1, 0, 0, 0, 1, 1],
                           [0, 0, 0, 0, 0, 0, 1, 0, 0],
                           [0, 0, 0, 0, 0, 0, 1, 0, 0], ])

        function_result = alg.edge_exchange(cost_matrix, constraints, 9, current_topology, current_degree)[0]

        self.assertTrue(np.array_equal(function_result, answer))

    # GENETIC ALGORITHM TESTS #

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

        prufer_encoding = alg.prufer_encode(tree)

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

        prufer_decoding, _ = alg.prufer_decode(np.array([8, 7, 0, 6, 7, 0, 6]))

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

        self.assertTrue(np.array_equal(alg.prufer_decode(
            alg.prufer_encode(tree))[0], tree))

        self.assertTrue(np.array_equal(alg.prufer_encode(
            alg.prufer_decode(prufer_encoding)[0]), prufer_encoding))

    def test_degree_checker(self):
        """
        Tests whether function that ensures tree encoded as Prufer number correctly determines if degree constraints
        have been violated when DCMST constructed.
        """
        # Prufer Number of Tree
        prufer_encoding = np.array([1, 3, 5, 2, 1, 0, 4], dtype=np.int32)

        self.assertTrue(alg.check_degree(prufer_encoding, np.array([3 for _ in range(9)]), 9))

        self.assertFalse(alg.check_degree(prufer_encoding, np.array([2 for _ in range(9)]), 9))

    def test_genetic_algorithm(self):
        """
        Checks that the DCMST returned by the genetic algorithm is a valid solution (no degree constraints violated, the
        tree is connected, and there exists a path between all vertex or satellite pairs).
        """
        cost_matrix = np.array([[0, 2.24, 2.24, 3.61, 6.71, 3.0, 5.39, 8.0, 9.43],
                                [2.24, 0, 2.0, 2.0, 4.47, 2.83, 4.0, 7.28, 7.62],
                                [2.24, 2.0, 0, 4.0, 5.66, 4.47, 6.0, 9.22, 9.49],
                                [3.61, 2.0, 4.0, 0, 4.0, 2.0, 2.0, 5.39, 5.83],
                                [6.71, 4.47, 5.66, 4.0, 0, 6.0, 4.47, 7.81, 5.1],
                                [3.0, 2.83, 4.47, 2.0, 6.0, 0, 2.83, 5.0, 7.07],
                                [5.39, 4.0, 6.0, 2.0, 4.47, 2.83, 0, 3.61, 4.24],
                                [8.0, 7.28, 9.22, 5.39, 7.81, 0, 3.61, 0, 5.0],
                                [9.43, 7.62, 9.49, 5.83, 5.1, 7.07, 0, 5.0, 0], ])

        constraints = np.array([3, 3, 3, 3, 3, 3, 3, 3, 3])

        num_sat = 9

        tree, degree = alg.genetic_algorithm(cost_matrix, constraints, num_sat, population_size=7)

        # Check degree of each vertex in resulting tree
        self.assertTrue(False not in (degree <= constraints).tolist())

        # Finds the number of connected components in returned tree (should be 1)
        n, _ = connected_components(tree, directed=False, connection='strong')

        # Check tree is connected
        self.assertTrue(n == 1)

        # Check all vertices are within tree
        self.assertTrue(False not in (np.sum(tree, axis=0) > 0).tolist())

    # ANT COLONY ALGORITHM TESTS #

    def test_ant_colony_algorithm(self):
        """
        Checks that the DCMST returned by the genetic algorithm is a valid solution (no degree constraints violated, the
        tree is connected, and there exists a path between all vertex or satellite pairs).
        """
        cost_matrix = np.array([[0, 2.24, 2.24, 3.61, 6.71, 3.0, 5.39, 8.0, 9.43],
                                [2.24, 0, 2.0, 2.0, 4.47, 2.83, 4.0, 7.28, 7.62],
                                [2.24, 2.0, 0, 4.0, 5.66, 4.47, 6.0, 9.22, 9.49],
                                [3.61, 2.0, 4.0, 0, 4.0, 2.0, 2.0, 5.39, 5.83],
                                [6.71, 4.47, 5.66, 4.0, 0, 6.0, 4.47, 7.81, 5.1],
                                [3.0, 2.83, 4.47, 2.0, 6.0, 0, 2.83, 5.0, 7.07],
                                [5.39, 4.0, 6.0, 2.0, 4.47, 2.83, 0, 3.61, 4.24],
                                [8.0, 7.28, 9.22, 5.39, 7.81, 0, 3.61, 0, 5.0],
                                [9.43, 7.62, 9.49, 5.83, 5.1, 7.07, 0, 5.0, 0], ])

        constraints = np.array([3, 3, 3, 3, 3, 3, 3, 3, 3])

        num_sat = 9

        tree, degree = alg.ant_colony(cost_matrix, constraints, num_sat)

        # Check degree of each vertex in resulting tree
        self.assertTrue(False not in (degree <= constraints).tolist())

        # Finds the number of connected components in returned tree (should be 1)
        n, _ = connected_components(tree, directed=False, connection='strong')

        # Check tree is connected
        self.assertTrue(n == 1)

        # Check all vertices are within tree
        self.assertTrue(False not in (np.sum(tree, axis=0) > 0).tolist())

    # GENERAL FUNCTION TESTS #

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

        current_topology, current_degree = alg.modified_prims_algorithm(test_costs,
                                                                        np.array(
                                                                            [3, 3, 3, 3,
                                                                             3, 3, 3, 3,
                                                                             3]), 9, 0)

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
            alg.increase_connectivity(current_topology, np.array([3, 3, 3, 3, 3, 3, 3, 3, 3]), current_degree,
                                      test_costs, 9), answer))


if __name__ == '__main__':
    unittest.main()

# References
# Example Graph - https://www.geeksforgeeks.org/prims-minimum-spanning-tree-mst-greedy-algo-5/
# Pytest FileNotFound error - https://stackoverflow.com/questions/56755761/filenotfounderror-when-using-python-m-pytest-
# vs-pytest
