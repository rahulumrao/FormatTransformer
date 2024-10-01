import os
import numpy as np
import argparse
from src.xyz_to_POSCAR import xyz_to_poscar  # Import the VASP converter
from src.xyz_to_QE import xyz_to_espresso    # Import the Quantum Espresso converter

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Convert XYZ to VASP POSCAR or Quantum Espresso input files")
    parser.add_argument("-i", "--input", type=str, default="centered_structure.xyz", help="Input XYZ file containing multiple frames")
    parser.add_argument("-b", "--box_size", type=float, nargs='+', default=[12.0], help="Box size in Å (default: Cubic box of 12.0)")
    parser.add_argument("-o", "--output_type", choices=["VASP", "QE"], required=True, help="Specify the output format: VASP or QE")
    parser.add_argument("-pp", "--pseudo", type=str, help="PseudoPotentials directory (required for Espresso)", default="../QE_PP")

    # Parse the arguments
    args = parser.parse_args()

    # Periodic box
    if len(args.box_size) == 1:
        # Cubic box
        box_a = box_b = box_c = args.box_size[0]
    elif len(args.box_size) == 3:
        # Orthorombic box
        box_a, box_b, box_c = args.box_size
    else:
        print('Error: parse 1 value for Cubic or 3 values for Orthorohmbic in box argument')
        exit(1)
    # Create cell vectors from A, B, C
    cell = [[box_a, 0.0, 0.0],
            [0.0, box_b, 0.0],
            [0.0, 0.0, box_c]]
    cell = np.array(cell)

    # Check which output format to use
    if args.output_type == "VASP":
        # Call the function to write POSCAR files
        print("Converting to VASP POSCAR format...")
        xyz_to_poscar(trajectory_file=args.input, cell=cell)
    elif args.output_type == "QE":
        # Call the function to write Quantum Espresso input files
        if not args.pseudo:
            print("Error: PseudoPotential directory is required for Quantum Espresso output.")
            exit(1)
        print("Converting to Quantum Espresso input format...")
        xyz_to_espresso(trajectory_file=args.input, cell=cell, pseudo_dir=args.pseudo)

if __name__ == "__main__":
    main()

