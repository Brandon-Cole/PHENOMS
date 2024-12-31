import sys
sys.path.append('../clean_package')
from utils import label_hbond, process_frames, extract_residue_numbers, plot_heatmap
from simulation_class import SimulationSet
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import tqdm
import timeit


pdb_files = [
    '../rep1_md_0_500_nolig_nojump_center.pdb',
    '../rep2_md_0_500_no_lig_nojump_center.pdb',
    '../rep3_md_0_500_no_lig_nojump_center.pdb'
]


sim_set = SimulationSet(pdb_files, resid_range=(50, 70), sub_frames=200)

# this will time the loading and processing step
start_time = timeit.default_timer()
all_pivot_tables = sim_set._load_and_process(n_jobs=3)
loading_time = timeit.default_timer() - start_time
print(f"Time taken to load and process data: {loading_time:.2f} seconds")

# this will time the plotting step
start_time = timeit.default_timer()
for i, pivot_table in enumerate(all_pivot_tables):
    print(f"Plotting Heatmap for Replicate {i+1}:")
    plot_heatmap(pivot_table, f"Replicate {i+1}")
plotting_time = timeit.default_timer() - start_time
print(f"Time taken to plot heatmaps: {plotting_time:.2f} seconds")

# broken: used to time the top hydrogen bonds plotting step
# start_time = timeit.default_timer()
# print("Plotting Top Hydrogen Bonds Above Threshold...")
# sim_set.generate_top_hbonds(threshold=0.5)  # Adjust threshold as needed
# top_hbonds_plotting_time = timeit.default_timer() - start_time
# print(f"Time taken to plot top hydrogen bonds: {top_hbonds_plotting_time:.2f} seconds")