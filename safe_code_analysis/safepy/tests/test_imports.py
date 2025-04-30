import unittest
import copy
import networkx as nx
import numpy as np

import safe
import safe_io


class TestImportCys(unittest.TestCase):

    def setUp(self):
        # Load the network
        f = '/Users/zhangxiang/Downloads/Research_CB/GIN/safepy_network_annotation/safe-data/networks/Costanzo_Science_2016.cys'
        self.graph = safe_io.load_network_from_cys(f, verbose=False)
        

    def test_size(self):
        num_nodes = nx.number_of_nodes(self.graph)
        self.assertEqual(num_nodes, 3971, "Should be 3971.")
        num_edges = nx.number_of_edges(self.graph)
        self.assertEqual(num_edges, 28202, "Should be 28202.")


class TestImportAttributes(unittest.TestCase):

    def setUp(self):
        pass

    def test_default(self):

        # Load the default network
        sf = safe.SAFE(verbose=False)
        sf.load_network()
        sf.load_attributes()

    def test_attribute_with_duplicate_values(self):

        # Load the default network
        sf = safe.SAFE(verbose=False)
        sf.load_network()

        f = '/Users/zhangxiang/Downloads/Research_CB/GIN/safepy_network_annotation/safe-data/tests/attribute_file_with_unmatched_duplicated_labels.txt'
        sf.load_attributes(attribute_file=f)


if __name__ == '__main__':
    unittest.main()
