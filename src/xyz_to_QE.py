# Import libraries
import os
import numpy as np
# Import ASE library for read and write
from tqdm import tqdm
from ase.io import read, write
from ase import Atoms
from ase.calculators.espresso import Espresso
# Import MDAnalysis for reading trajectory files
import MDAnalysis as mda

# #====================================================#
# def get_pseudopotential(pseudo_dir, atom_type):
#     '''
#       Find the PseudoPotential for each atom
#       example: pseudopotentials = {'O': 'O.pbe-n-kjpaw_psl.0.1.UPF',
#                                    'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
#                                    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF'}
#     '''
#     pseudopotentials = {}
#     for atom in atom_type:
#         for file in os.listdir(pseudo_dir):
#             if file.startswith(atom) and file.endwith('.UPF'):
#                 pseudopotentials[atoms] = os.path.join(pseudo_dir, file)
#                 break
#         else:
#             raise ValueError(f'Pseudopotential for {atom} not found in {pseudo_dir}')
#     return pseudopotentials

def xyz_to_espresso(trajectory_file, cell, pseudo_dir):
        
    pseudo_dict = {'Ac': 'Ac.us.z_11.ld1.psl.v1.0.0-high.upf', 'Os': 'Os_pbe_v1.2.uspp.F.UPF',
                       'Ag': 'Ag_ONCV_PBE-1.0.oncvpsp.upf', 'Pp': 'P.pbe-n-rrkjus_psl.1.0.0.UPF',
                       'Al': 'Al.pbe-n-kjpaw_psl.1.0.0.UPF','Pa': 'Pa.paw.z_13.ld1.uni-marburg.v0.upf',
                       'Am': 'Am.paw.z_17.ld1.uni-marburg.v0.upf','Pb': 'Pb.pbe-dn-kjpaw_psl.0.2.2.UPF',
                       'Ar': 'Ar_ONCV_PBE-1.1.oncvpsp.upf', 'Pd': 'Pd_ONCV_PBE-1.0.oncvpsp.upf',
                       'As': 'As.pbe-n-rrkjus_psl.0.2.UPF', 'Pm': 'Pm.paw.z_15.atompaw.wentzcovitch.v1.2.upf',
                       'At': 'At.us.z_17.ld1.psl.v1.0.0-high.upf', 'Po': 'Po.pbe-dn-rrkjus_psl.1.0.0.UPF',
                       'Au': 'Au_ONCV_PBE-1.0.oncvpsp.upf', 'Pr': 'Pr.paw.z_13.atompaw.wentzcovitch.v1.2.upf',
                       'Ba': 'Ba.pbe-spn-kjpaw_psl.1.0.0.UPF','Pu': 'Pu.paw.z_16.ld1.uni-marburg.v0.upf',
                       'Bi': 'Bi_pbe_v1.uspp.F.UPF', 'Ra': 'Ra.paw.z_20.ld1.psl.v1.0.0-high.upf',
                       'Bk': 'Bk.paw.z_19.ld1.uni-marburg.v0.upf', 'Rb': 'Rb_ONCV_PBE-1.0.oncvpsp.upf',
                       'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF', 'Re': 'Re_pbe_v1.2.uspp.F.UPF',
                       'Ca': 'Ca_pbe_v1.uspp.F.UPF', 'Rh': 'Rh_ONCV_PBE-1.0.oncvpsp.upf',
                       'Cd': 'Cd.pbe-dn-rrkjus_psl.0.3.1.UPF', 'Rn': 'Rn.pbe-dn-kjpaw_psl.1.0.0.UPF',
                       'Ce': 'Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf','Ru': 'Ru_ONCV_PBE-1.0.oncvpsp.upf',
                       'Cf': 'Cf.paw.z_20.ld1.uni-marburg.v0.upf','Sc': 'Sc_ONCV_PBE-1.0.oncvpsp.upf',
                       'Cm': 'Cm.paw.z_18.ld1.uni-marburg.v0.upf', 'Se': 'Se_pbe_v1.uspp.F.UPF',
                       'Co': 'Co_pbe_v1.2.uspp.F.UPF','Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
                       'Cs': 'Cs_pbe_v1.uspp.F.UPF', 'Sm': 'Sm.paw.z_16.atompaw.wentzcovitch.v1.2.upf',
                       'Cu': 'Cu.paw.z_11.ld1.psl.v1.0.0-low.upf', 'Sn': 'Sn_pbe_v1.uspp.F.UPF',
                       'Dy': 'Dy.paw.z_20.atompaw.wentzcovitch.v1.2.upf', 'Sr': 'Sr_pbe_v1.uspp.F.UPF',
                       'Er': 'Er.paw.z_22.atompaw.wentzcovitch.v1.2.upf', 'Ta': 'Ta_pbe_v1.uspp.F.UPF',
                       'Es': 'Es.paw.z_21.ld1.uni-marburg.v0.upf', 'Tb': 'Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf',
                       'Eu': 'Eu.paw.z_17.atompaw.wentzcovitch.v1.2.upf', 'Tc': 'Tc_ONCV_PBE-1.0.oncvpsp.upf',
                       'Fe': 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF', 'Te': 'Te_pbe_v1.uspp.F.UPF',
                       'Fm': 'Fm.paw.z_22.ld1.uni-marburg.v0.upf', 'Th': 'Th.paw.z_12.ld1.uni-marburg.v0.upf',
                       'Fr': 'Fr.paw.z_19.ld1.psl.v1.0.0-high.upf', 'Tl': 'Tl_pbe_v1.2.uspp.F.UPF',
                       'Ga': 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF', 'Tm': 'Tm.paw.z_23.atompaw.wentzcovitch.v1.2.upf',
                       'Gd': 'Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf', 'U': 'U.paw.z_14.ld1.uni-marburg.v0.upf',
                       'H': 'H.pbe-rrkjus_psl.1.0.0.UPF', 'W': 'W_pbe_v1.2.uspp.F.UPF',
                       'He': 'He_ONCV_PBE-1.0.oncvpsp.upf', 'Xe': 'Xe_ONCV_PBE-1.1.oncvpsp.upf',
                       'Hf': 'Hf-sp.oncvpsp.upf', 'Y': 'Y_pbe_v1.uspp.F.UPF',
                       'Hg': 'Hg_ONCV_PBE-1.0.oncvpsp.upf', 'Yb': 'Yb.paw.z_24.atompaw.wentzcovitch.v1.2.upf',
                       'Ho': 'Ho.paw.z_21.atompaw.wentzcovitch.v1.2.upf', 'Zn': 'Zn_pbe_v1.uspp.F.UPF',
                       'I': 'I.pbe-n-kjpaw_psl.0.2.UPF', 'Zr': 'Zr_pbe_v1.uspp.F.UPF',
                       'In': 'In.pbe-dn-rrkjus_psl.0.2.2.UPF', 'B': 'b_pbe_v1.4.uspp.F.UPF',
                       'Ir': 'Ir_pbe_v1.2.uspp.F.UPF', 'Be': 'be_pbe_v1.4.uspp.F.UPF',
                       'K': 'K.pbe-spn-kjpaw_psl.1.0.0.UPF', 'Br': 'br_pbe_v1.4.uspp.F.UPF',
                       'Kr': 'Kr_ONCV_PBE-1.0.oncvpsp.upf', 'Cl': 'cl_pbe_v1.4.uspp.F.UPF',
                       'La': 'La.paw.z_11.atompaw.wentzcovitch.v1.2.upf', 'Cr': 'cr_pbe_v1.5.uspp.F.UPF',
                       'Lr': 'Lr.paw.z_25.ld1.uni-marburg.v0.upf', 'F': 'f_pbe_v1.4.uspp.F.UPF',
                       'Lu': 'Lu.paw.z_25.atompaw.wentzcovitch.v1.2.upf', 'Ge': 'ge_pbe_v1.4.uspp.F.UPF',
                       'Md': 'Md.paw.z_23.ld1.uni-marburg.v0.upf', 'Li': 'li_pbe_v1.4.uspp.F.UPF',
                       'Mg': 'Mg.pbe-n-kjpaw_psl.0.3.0.UPF', 'Mn': 'mn_pbe_v1.5.uspp.F.UPF',
                       'Mo': 'Mo_ONCV_PBE-1.0.oncvpsp.upf', 'Na': 'na_pbe_v1.5.uspp.F.UPF',
                       'N': 'N.pbe-n-radius_5.UPF', 'Ni': 'ni_pbe_v1.4.uspp.F.UPF',
                       'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF', 'Pt': 'pt_pbe_v1.4.uspp.F.UPF',
                       'Nd': 'Nd.paw.z_14.atompaw.wentzcovitch.v1.2.upf',}

    # check in Input file exist
    if not os.path.isfile(trajectory_file):
        print(f'Error: Input file {trajectory_file} not found, exiting...')
        exit(1)
    prefix_path = os.getcwd()
    print("=============================================================")
    print(f' Working directory : {prefix_path}')
    print(f' Reading file : {trajectory_file}')
    print(f' Box dimenstions : A = {cell[0,0]}, B = {cell[1,1]}, C = {cell[2,2]}')
    print("=============================================================")

#    # Create cell vectors from A, B, C
#    cell = [[box[0], 0.0, 0.0],
#            [0.0, box[1], 0.0],
#            [0.0, 0.0, box[2]]]
#
    # Output directory for Quantum Espresso input files
    output_dir = 'In_file'
    os.makedirs(output_dir, exist_ok=True)

    # Read the trajectory
    universe = mda.Universe(trajectory_file)

    # get atom names
    atom_types = np.unique(universe.atoms.names)
    # get pseudopotentials for the atom types
    # pseudopotentials = get_pseudopotential(pseudo_dir, atom_types)

    # Iterate over each frame
    for ts in tqdm(universe.trajectory):
        # Convert MDAnalysis universe to ASE atoms
        atoms = Atoms(symbols=universe.atoms.names,
                      positions=universe.atoms.positions,
                      cell=cell,
                      pbc=True)
        outdir = os.path.join('scratch', str(ts.frame))
        # Write Quantum Espresso input file for the current frame
        input_qe = {
            'calculation': 'scf',
            'outdir': outdir,
            'pseudo_dir': pseudo_dir,
            'tprnfor': True,
            'tstress': True,
            'disk_io': 'none',
            'system': {
                'ecutwfc': 30,
                'input_dft': 'PBE'
            },
            'electrons': {
                'mixing_beta': 0.5,
                'electron_maxstep': 1000
            }
        }
        if (ts==0):
            print(f"# Converting xyz to QE input using MDAnalysis for structure {ts.frame}")
            print(f"# Be advised: the PseudoPotentials are for PBE level DFT.")

        write(os.path.join(output_dir, f"pw-wt-{ts.frame}.in"), atoms, format='espresso-in',
              input_data=input_qe, pseudopotentials=pseudo_dict)

#===============================================================================================#
