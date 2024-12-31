from utils import load_and_select_residues, process_frames, extract_residue_numbers, create_pivot_table, plot_heatmap, plot_top_hbonds

class SimulationSet:
    def __init__(self, pdb_files, resid_range=None, sub_frames=None):
        """
        Initialize the simulation set, managing replicates within the same system.

        Parameters:
        - pdb_files: List of PDB file paths for the replicates.
        - resid_range: Residue range for H-bond analysis (optional).
        - sub_frames: Number of frames to analyze (optional).
        """
        self.pdb_files = pdb_files
        self.resid_range = resid_range
        self.sub_frames = sub_frames
        self.simulation_data = []

    def _load_and_process(self, n_jobs=1):
        """
        Load and process data for each PDB file in the simulation set.

        This method processes the trajectory data, extracts hydrogen bond information, 
        and creates pivot tables for each simulation.

        Parameters:
        - n_jobs: Number of jobs to use for parallel processing (optional, default is 1).

        Returns:
        - all_pivot_tables: A list of pivot tables created from hydrogen bond data.
        """
        all_hbonds_dfs = []  
        all_bond_labels = set()  

        for pdb_file in self.pdb_files:
            trajectory = load_and_select_residues(pdb_file, self.resid_range)
            hbonds_df = process_frames(trajectory, self.sub_frames, n_jobs=n_jobs)
            all_bond_labels.update(hbonds_df['Bond Label'].unique())
            all_hbonds_dfs.append(hbonds_df)

        bond_labels_sorted = sorted(all_bond_labels, key=lambda x: extract_residue_numbers(x))

        all_pivot_tables = []
        for pdb_file, hbonds_df in zip(self.pdb_files, all_hbonds_dfs): 
            hbonds_pivot = create_pivot_table(hbonds_df, bond_labels_sorted)
            all_pivot_tables.append(hbonds_pivot)

        self.simulation_data = [(all_pivot_tables, bond_labels_sorted)]

        return all_pivot_tables

    def generate_heatmap(self):
        """
        Generate heatmaps for all simulations in the set.

        This method creates a heatmap for each simulation based on the pivot tables 
        generated from the hydrogen bond data.

        It uses the `plot_heatmap` function from utils to plot the data.
        """
        if not self.simulation_data:
            raise ValueError("Simulation data is empty. Please process the data first.")
        
        for all_pivot_tables, bond_labels_sorted in self.simulation_data:
            for i, pivot_table in enumerate(all_pivot_tables):
                file_name = self.pdb_files[i] if i < len(self.pdb_files) else f"Simulation_{i+1}"
                plot_heatmap(pivot_table, file_name)

    def generate_top_hbonds(self, threshold=0.5):
        """
        Plot a bar chart showing the top hydrogen bonds with average lifetimes
        above the given threshold.
        """
        if not self.simulation_data:
            raise ValueError("Simulation data is empty. Please process the data first.")
        
        all_pivot_tables = []
        for all_pivot_tables_per_sim, bond_labels_sorted in self.simulation_data:
            print(f"Processing pivot tables for one simulation: {len(all_pivot_tables_per_sim)}")
            all_pivot_tables.extend(all_pivot_tables_per_sim)
        
        print(f"Total pivot tables collected: {len(all_pivot_tables)}")
        
        if not all_pivot_tables:
            print("No pivot tables to plot. Exiting function.")
            return
        
        plot_top_hbonds(all_pivot_tables, threshold, self.sub_frames)



#     def generate_comparative_plot(self):
#         """
#         Generate a comparative plot (e.g., average heatmap) for the simulations.
#         """
#         # Placeholder for generating comparative plots across replicates
#         pass


# class ComparisonSet:
#     def __init__(self, simulation_sets):
#         """
#         Initialize a comparison set with multiple simulation sets for comparing different systems.

#         Parameters:
#         - simulation_sets: List of `SimulationSet` objects for different systems.
#         """
#         self.simulation_sets = simulation_sets

#     def compare_heatmaps(self):
#         """
#         Generate comparative heatmaps comparing the simulation sets.
#         """
#         for sim_set in self.simulation_sets:
#             sim_set.generate_heatmap()
#         # Generate plots comparing heatmaps across systems
#         pass

#     def compare_barplots(self):
#         """
#         Generate comparative bar plots comparing the simulation sets.
#         """
#         for sim_set in self.simulation_sets:
#             sim_set.generate_comparative_plot()
#         # Generate bar plots comparing replicates or systems
#         pass

#     def summarize_comparisons(self):
#         """
#         Summarize the comparison of multiple simulation sets.
#         """
#         # Generate any relevant statistics or summary for comparison
#         pass




# class Simulation:
#    #res range not mand

#    def __init__(self, pdb_files, resid_range, sub_frames):
#         """
#         Initialize the simulation with necessary files and parameters.

#         Parameters:
#         - pdb_files: List of PDB files for the simulation
#         - resid_range: Residue range for H-bond analysis
#         - sub_frames: Number of frames to analyze
#         """
#         self.name = exp_name
#         self.pdb_files = pdb_files
#         if resid_range is None:
#             self.resid_range = 'all'
#         self.resid_range = resid_range
#         self.sub_frames = sub_frames

#     def replicate_heatmaps(self):
#         pass

#     def replicate_barplots(self):
#         pass

#     def comparison_heatmaps(self):
#         pass

#     def comparison_barplots(self):
#         pass

#     def entire_protein_heatmaps(self):
#         pass

#     def region_protein_peptide_map(self):
#         pass

#     def region_protein_peptide_comparison(self):
#         pass

#     def entire_protein_peptide_map(self):
#         pass

#     def entire_protein_peptide_comparison(self):
#         pass

    