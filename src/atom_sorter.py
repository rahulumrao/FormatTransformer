import numpy as np
from ase.io import read, write

def traj_sort(traj_file=trajectory_file):

    atoms_trajectory = read(traj_file, index=':')  # ':' reads all frames

    # Sort atoms alphabetically for each frame
    sorted_traj = []
    for atoms in atoms_trajectory:
        # 
        atom_type = np.array(atoms.get_chemical_symbols())
    
        # sorting atoms based on alphabetical order
        sorted_indices = atom_type.argsort()
    
        # sorting atoms
        sorted_atoms = atoms[sorted_indices]
        sorted_traj.append(sorted_atoms)

    return sorted_traj
# Write trajectory to a new XYZ file
# write('sorted_output_trajectory.xyz', sorted_traj)
