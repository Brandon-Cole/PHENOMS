
import mdtraj as md
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from joblib import Parallel, delayed


def load_and_select_residues(pdb_file, resid_range=None):
    """
    Load the trajectory from a PDB file and select residues based on the provided range.
    If no range is provided, all residues are selected.

    Parameters:
    - pdb_file: Path to the PDB file.
    - resid_range: Tuple containing the start and end residue numbers for selection (optional).

    Returns:
    - A sliced trajectory containing the selected atoms.
    """
    trajectory = md.load(pdb_file)
    
    if resid_range is None:
        selected_atoms = trajectory.top.select('protein')  # selects all protein atoms
    else:
        selected_atoms = trajectory.top.select(f'resid {resid_range[0]} to {resid_range[1]}')
    
    return trajectory.atom_slice(selected_atoms)


def label_hbond(hbond, trajectory):
    """
    Label a hydrogen bond between two atoms in a given trajectory frame based on their names.
    
    Parameters:
    - hbond: A tuple containing the indices of the donor, hydrogen, and acceptor atoms involved in the hydrogen bond.
    - trajectory: The MDTraj trajectory object containing the simulation data.

    Returns:
    - A string representing the label of the hydrogen bond in the format 'residue1 -- residue2' if it's a valid H-bond,
      or None if the bond is not between a nitrogen and oxygen (the typical case for hydrogen bonds).
    """
    atom1 = trajectory.topology.atom(hbond[0])
    atom2 = trajectory.topology.atom(hbond[2])
    atom1_name = atom1.name
    atom2_name = atom2.name
    
    if (atom1_name == 'N' and atom2_name == 'O') or (atom1_name == 'O' and atom2_name == 'N'):
        residue1 = atom1.residue
        residue2 = atom2.residue
        return f'{residue1} -- {residue2}'
    else:
        return None


def process_frame(frame, trajectory):
    """
    Process a single frame to detect hydrogen bonds and label them.
    
    Parameters:
    - frame: The frame index to process.
    - trajectory: The MDTraj trajectory object containing the simulation data.

    Returns:
    - A DataFrame containing the hydrogen bonds detected in the frame, with columns for 'Donor Index', 
      'Hydrogen Index', 'Acceptor Index', 'Bond Label', and 'Frame'.
    """
    hbonds = md.baker_hubbard(trajectory[frame], periodic=False, freq=0)
    frame_df = pd.DataFrame(hbonds, columns=['Donor Index', 'Hydrogen Index', 'Acceptor Index'])
    
    bond_labels = [label_hbond(bond, trajectory[frame]) for i, bond in frame_df.iterrows()]
    frame_df['Bond Label'] = bond_labels
    frame_df['Frame'] = frame
    
    return frame_df


def process_frames(trajectory, sub_frames=None, n_jobs=1):
    """
    Process the trajectory frames to detect hydrogen bonds and label them in parallel.

    Parameters:
    - trajectory: The MDTraj trajectory object containing the simulation data.
    - sub_frames: The number of frames to process. If None, the entire trajectory is used.
    - n_jobs: The number of CPU cores to use for parallel processing. Default is -1 (use all available cores).

    Returns:
    - A DataFrame containing the hydrogen bonds detected in the processed frames, with columns for 'Donor Index', 
      'Hydrogen Index', 'Acceptor Index', 'Bond Label', and 'Frame'.
    """
    total_frames = trajectory.n_frames
    
    if sub_frames is None:
        sub_frames = total_frames

    hbond_frames = Parallel(n_jobs=n_jobs)(
        delayed(process_frame)(frame, trajectory) for frame in tqdm(range(min(sub_frames, total_frames)), desc="Processing Frames", unit="frame")
    )

    hbonds_df = pd.concat(hbond_frames, ignore_index=True)
    hbonds_df = hbonds_df.dropna(subset=['Bond Label'])

    return hbonds_df

def extract_residue_numbers(bond_label):
    """
    Extract residue numbers from a bond label in the format 'Residue1 -- Residue2'.

    Parameters:
    - bond_label (str): A string representing a bond label, e.g., 'Ala1 -- Gly2'.

    Returns:
    - tuple: A tuple containing the residue numbers of the two residues involved in the bond.
      (first_residue_number, second_residue_number)
    """
    residues = bond_label.split(' -- ')
    first_residue_number = int(''.join(filter(str.isdigit, residues[0])))
    second_residue_number = int(''.join(filter(str.isdigit, residues[1])))
    return first_residue_number, second_residue_number


def create_pivot_table(hbonds_df, bond_labels_sorted):
    """
    Create a pivot table from the hydrogen bond data, with consistent bond labels sorted.

    Parameters:
    - hbonds_df (pd.DataFrame): DataFrame containing hydrogen bond information with columns 
      'Bond Label' and 'Frame'.
    - bond_labels_sorted (list): A sorted list of bond labels to use as the index in the pivot table.

    Returns:
    - pd.DataFrame: A pivot table with bond labels as rows and frames as columns, showing the count 
      of occurrences for each bond label per frame.
    """
    hbonds_pivot = hbonds_df.pivot_table(index='Bond Label', columns='Frame', aggfunc='size', fill_value=0)
    hbonds_pivot = hbonds_pivot.reindex(bond_labels_sorted, fill_value=0)
    return hbonds_pivot


def plot_heatmap(pivot_table, file_name):
    """
    Plot a heatmap of hydrogen bond occupation over time.

    Parameters:
    - pivot_table (pd.DataFrame): DataFrame representing the hydrogen bond data, 
      with bond labels as the index and frames as columns.
    - file_name (str): Name of the file used for the plot title.
    
    The heatmap displays hydrogen bond occupancy across frames, with each bond represented 
    as a small square filled with a pastel red color based on the occupation value.
    """
 
    plt.figure(figsize=(20, 12)) 
    sns.heatmap(pivot_table, cmap='bwr', cbar=False, square=False, linewidths=0.5, linecolor='black') 
    plt.title(f'H-Bond Occupation over Time - {file_name}', weight='bold')
    plt.xlabel('Frame (ns)', weight='bold', labelpad=10)
    plt.ylabel('Bond Label', weight='bold', labelpad=10)
    plt.xticks(ticks=np.arange(0, pivot_table.shape[1], 10), labels=np.arange(0, pivot_table.shape[1], 10), rotation=90, fontsize=8)
    plt.yticks(rotation=0, fontsize=6)
    plt.tight_layout()
    plt.show()


def plot_top_hbonds(all_pivot_tables, threshold, sub_frames=None):
    """
    Plot a bar chart showing the top hydrogen bonds with average lifetimes
    above the given threshold.
    """
    if not all_pivot_tables:
        print("Empty pivot tables received.")
        return

    print(f"Received {len(all_pivot_tables)} pivot tables for plotting.")

    combined_pivot = pd.concat(all_pivot_tables, axis=1).fillna(0)
    combined_pivot['Average Lifetime'] = combined_pivot.mean(axis=1)

    total_frames = sub_frames or max([table.shape[1] for table in all_pivot_tables])
    combined_pivot['Presence Fraction'] = (combined_pivot > 0).sum(axis=1) / total_frames


    significant_bonds = combined_pivot[combined_pivot['Presence Fraction'] > threshold]
    print(f"Significant bonds above threshold: {significant_bonds.shape[0]}")
    if significant_bonds.empty:
        print("No significant bonds to plot.")
        return

    significant_bonds = significant_bonds.sort_values(by='Average Lifetime', ascending=False)
    significant_bonds['Average Lifetime'].plot(kind='bar', figsize=(10, 6), color='skyblue')
    plt.title(f"Top Hydrogen Bonds (Threshold: {threshold})")
    plt.xlabel("Hydrogen Bond")
    plt.ylabel("Average Lifetime")
    plt.tight_layout()
    plt.show()



def process_replicates(pdb_files, resid_range, sub_frames):
    all_hbonds_dfs = []
    all_bond_labels = set()

    for pdb_file in pdb_files:
        trajectory = load_and_select_residues(pdb_file, resid_range)
        hbonds_df = process_frames(trajectory, sub_frames)
        all_bond_labels.update(hbonds_df['Bond Label'].unique())
        all_hbonds_dfs.append(hbonds_df)

    bond_labels_sorted = sorted(all_bond_labels, key=lambda x: extract_residue_numbers(x))

    all_pivot_tables = []
    for pdb_file, hbonds_df in zip(pdb_files, all_hbonds_dfs):
        hbonds_pivot = create_pivot_table(hbonds_df, bond_labels_sorted)
        all_pivot_tables.append(hbonds_pivot)
        plot_heatmap(hbonds_pivot, pdb_file)

    return all_pivot_tables, bond_labels_sorted