#!/usr/bin/env python
import os
import numpy as np
from tqdm import tqdm
import argparse
from ase.io import write
from ase import Atoms
import MDAnalysis as mda
from ase.io import read, write
#============================================

def traj_sort(traj_file='traj.xyz'):
    out_file_name = 'output_trajectory.xyz'
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

    # Write trajectory to a new XYZ file
    write(out_file_name, sorted_traj)
    return out_file_name
#============================================
def xyz_to_poscar(trajectory_file, cell):

#    # Periodic box
#    if len(args.box_size) == 1:
#        # Cubic box
#        box_a = box_b = box_c = args.box_size[0]
#    elif len(args.box_size) == 3:
#        # Orthorombic box
#        box_a, box_b, box_c = args.box_size
#    else:
#        print('Error: parse 1 value for Cubic or 3 values for Orthorohmbic in box argument')
#        exit(1)

    # check in Input file exist
    if not os.path.isfile(trajectory_file):
        print(f'Error: Input file {trajectory_file} not found, exiting...')
        exit(1)
    prefix_path = os.getcwd()
    print("=============================================================")
    print(f' Working directory : {prefix_path}')
    print(f' Reading file : {trajectory_file}')
    print(f' Box dimensions : A = {cell[0,0]}, B = {cell[1,1]}, C = {cell[2,2]}')
    print("=============================================================")

#    # Create cell vectors from A, B, C
#    cell = [[box[0], 0.0, 0.0],
#            [0.0, box[1], 0.0],
#            [0.0, 0.0, box[2]]]
    
    # Sort atoms in trajectory file alphabatically for VASP
    trajectory_file = traj_sort(trajectory_file)
    # Output directory for Quantum Espresso input files
    output_dir = 'poscar_files'

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Read the trajectory
    universe = mda.Universe(trajectory_file)

    # Iterate over each frame
    for ts in tqdm(universe.trajectory):
        # Convert MDAnalysis universe to ASE atoms
        if (ts==0):
            print(f"# Converting xyz to POSCAR using MDAnalysis for structure {ts}")
        atoms = Atoms(symbols=universe.atoms.names,
                      positions=universe.atoms.positions,
                      cell=cell,
                      pbc=True)
    #    print(ts.frame)
        outdir = os.path.join(output_dir, str(ts.frame))
        os.makedirs(outdir, exist_ok=True)
        write(os.path.join(outdir,f"POSCAR"), atoms, format='vasp')

    print(f'POSCAR files written in directory, {output_dir}.')
#==================================== End Program ====================================================#
