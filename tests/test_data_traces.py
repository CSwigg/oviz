import unittest
import pandas as pd
import numpy as np
from cluster_traces import StarClusterData, StarClusterCollection

class TestStarClusterData(unittest.TestCase):
    def setUp(self):
        # Create a StarClusterData instance with some dummy data
        data = {
            'x': [1, 2, 3],
            'y': [4, 5, 6],
            'z': [7, 8, 9],
            'vx': [10, 11, 12],
            'vy': [13, 14, 15],
            'vz': [16, 17, 18],
            'name': ['a', 'b', 'c'],
            'age_myr': [19, 20, 21]
        }
        df = pd.DataFrame(data)
        self.cluster_data = StarClusterData(df)
        self.time = np.arange(0, 60e6, 1e6)

    def test_init(self):
        self.assertIsInstance(self.cluster_data, StarClusterData)

    def test_create_integrated_dataframe(self):
        # First, call integrate_orbits to populate cluster_int_coords
        self.cluster_data.integrate_orbits(self.time)

        # Then, call create_integrated_dataframe and check if it returns a DataFrame
        result = self.cluster_data.create_integrated_dataframe(self.time)
        self.assertIsInstance(result, pd.DataFrame)

    def test_integrate_orbits(self):
        # Call integrate_orbits and check if it populates cluster_int_coords and df_int
        self.cluster_data.integrate_orbits(self.time)
        self.assertIsNotNone(self.cluster_data.cluster_int_coords)
        self.assertIsNotNone(self.cluster_data.df_int)


class TestStarClusterCollection(unittest.TestCase):
    def setUp(self):
        # Create a StarClusterCollection instance
        self.collection = StarClusterCollection()

        # Create some StarClusterData instances and add them to the collection
        data = {
            'x': [1, 2, 3],
            'y': [4, 5, 6],
            'z': [7, 8, 9],
            'U': [10, 11, 12],
            'V': [13, 14, 15],
            'W': [16, 17, 18],
            'name': ['a', 'b', 'c'],
            'age_myr': [19, 20, 21]
        }
        df = pd.DataFrame(data)
        cluster_data1 = StarClusterData(df)
        cluster_data2 = StarClusterData(df)
        self.collection.add_cluster(cluster_data1)
        self.collection.add_cluster(cluster_data2)

    def test_add_cluster(self):
        # Test that the add_cluster method adds a cluster to the collection
        old_length = len(self.collection.get_all_clusters())
        data = {
            'x': [1, 2, 3],
            'y': [4, 5, 6],
            'z': [7, 8, 9],
            'U': [10, 11, 12],
            'V': [13, 14, 15],
            'W': [16, 17, 18],
            'name': ['a', 'b', 'c'],
            'age_myr': [19, 20, 21]
        }
        df = pd.DataFrame(data)
        cluster_data = StarClusterData(df)
        self.collection.add_cluster(cluster_data)
        new_length = len(self.collection.get_all_clusters())
        self.assertEqual(old_length + 1, new_length)

    def test_get_cluster(self):
        # Test that the get_cluster method returns the correct cluster
        cluster = self.collection.get_cluster(0)
        self.assertIsInstance(cluster, StarClusterData)

    def test_get_all_clusters(self):
        # Test that the get_all_clusters method returns all clusters
        clusters = self.collection.get_all_clusters()
        self.assertEqual(len(clusters), 2)
        for cluster in clusters:
            self.assertIsInstance(cluster, StarClusterData)

if __name__ == '__main__':
    unittest.main()
if __name__ == '__main__':
    unittest.main()