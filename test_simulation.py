import unittest
from simulation_class import SimulationSet
from utils import plot_heatmap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class TestSimulationSet(unittest.TestCase):

    def setUp(self):
        """
        Set up test environment before each test.
        You can initialize necessary objects here.
        """
        self.pdb_files = [
            '../rep1_md_0_500_nolig_nojump_center.pdb',
            '../rep2_md_0_500_no_lig_nojump_center.pdb',
            '../rep3_md_0_500_no_lig_nojump_center.pdb'
        ]
        self.sim_set = SimulationSet(self.pdb_files, resid_range=(50, 70), sub_frames=200)

    def test_initialization(self):
        """
        Test the initialization of the SimulationSet class.
        """
        self.assertEqual(len(self.sim_set.pdb_files), 3)
        self.assertEqual(self.sim_set.resid_range, (50, 70))
        self.assertEqual(self.sim_set.sub_frames, 200)

    def test_load_and_process(self):
        """
        Test the loading and processing of PDB files and generation of pivot tables.
        """
        all_pivot_tables = self.sim_set._load_and_process(n_jobs=3)
        self.assertIsInstance(all_pivot_tables, list)
        self.assertGreater(len(all_pivot_tables), 0)  
        self.assertIsInstance(all_pivot_tables[0], pd.DataFrame)

    def test_plot_heatmap(self):
        """
        Test the plotting of heatmaps for the given pivot tables.
        """
        all_pivot_tables = self.sim_set._load_and_process(n_jobs=3)
        
        try:
            for i, pivot_table in enumerate(all_pivot_tables):
                with self.subTest(i=i):
                    print(f"Plotting Heatmap for Replicate {i+1}:")
                    plot_heatmap(pivot_table, f"Replicate {i+1}")
                    plt.close()  
        except Exception as e:
            self.fail(f"Heatmap plotting failed: {e}")

def test_generate_top_hbonds(self):
    """
    Test the generate_top_hbonds method for plotting hydrogen bonds.
    """
    try:
        self.sim_set.process_data() 
        self.sim_set.generate_top_hbonds(threshold=0.5)
    except Exception as e:
        self.fail(f"Top hydrogen bonds plotting failed: {e}")

if __name__ == '__main__':
    unittest.main()
