"""Compute energies between all pairs of atoms in a structure"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas
import numpy as np
import math
sys.path.append('/home/haddox/')
import pyrosetta
import pyrosetta.rosetta

flags_f = '/home/haddox/2019/optE_eval/data/HH_run19A_flags_266'
weights_f = '/home/haddox/2019/optE_eval/data/HH_run19A_weights_266.wts'
flags = [
    '-mute all', '-mute core', '-mute protocols',
    '-ignore_unrecognized_res',
    '-read_only_ATOM_entries', f'@{flags_f}'
]
pyrosetta.init(extra_options=' '.join(flags))
sf = pyrosetta.create_score_function(weights_f)
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(sf)

# Define functions
def compute_burial(pdb, probe_radius=6.5, resolution=0.5, poly_x='L'):
    """
    Compute burial of each atom in an input PDB
    
    Args:
        *pdb*: the path to an input PDB file
        *probe_radius*: the radius of the probe used to
            trace the molecular surface
        *resolution*: the resolution used to construct
            the gride for computing depth
    
    Returns: a CSV with per-atom burial
    """
    
    # Read in PDB as pose and score the pose
    pose = pyrosetta.pose_from_pdb(pdb)
    sf(pose)

    # Initiate object for computing SASA
    sasa_calc = \
        pyrosetta.rosetta.protocols.vardist_solaccess.VarSolDistSasaCalculator()
    sasa_map = sasa_calc.calculate(pose)

    # Initiate object for computing atomic depth. If indicated, mutate pose to
    # poly-X before doing so. Otherwise, it will be mutated to poly-leucine
    # using the atomic-depth C++ code
    if poly_x != 'L':
        aa_x = pyrosetta.rosetta.core.conformation.get_residue_from_name1(poly_x)
        poly_x_pose = pyrosetta.pose_from_pdb(pdb)
        for i in range(1, poly_x_pose.size()+1):
            poly_x_pose.replace_residue(i, aa_x, True)
        atomic_depth = pyrosetta.rosetta.core.scoring.atomic_depth.AtomicDepth(
            pose=poly_x_pose,
            probe_radius=probe_radius,
            poly_leu_depth=False,
            resolution=resolution,
            padding=4.0
        )
    else:
        atomic_depth = pyrosetta.rosetta.core.scoring.atomic_depth.AtomicDepth(
            pose=pose,
            probe_radius=probe_radius,
            poly_leu_depth=True,
            resolution=resolution,
            padding=4.0
        )

    # Initiate a dictionary for storing inter-atomic
    # distances
    burial_dict = {
        key : []
        for key in [
            'res_n', 'res_pdb_n', 'res_name', 'res_chain',
            'res_n_nbrs',
            'atom_n', 'atom_name', 'atom_type_name',
            'atom_bb', 'atom_lj_radius', 'atom_sasa',
            'atom_depth', 'atom_x', 'atom_y', 'atom_z',
        ]
    }

    # Loop over all residues and atoms in the protein
    pose_size = pose.size()
    res_ns = list(range(1, pose_size+1))
    for res_n in res_ns:

        # Compute number of C-beta neighbors for given residue
        res = pose.residue(res_n)
        neighbors = pyrosetta.rosetta.core.select.get_tenA_neighbor_residues(
            pose,
            pyrosetta.Vector1([
                i == res_n for i in range(1, pose.size()+1)
            ])
        )
        res_n_nbrs = sum(neighbors) - 1
        
        # Get the atom type set for given residue
        ats = res.type().atom_type_set()
            
        # Loop over all atoms in residue
        for atom_n in list(range(1, res.natoms()+1)):

            # Compute atom depth
            atom = res.atom(atom_n)
            atom_type = res.atom_type(atom_n)
            atom_radius = atom_type.lj_radius()
            xyz = res.xyz(atom_n)
            depth = atomic_depth.calcdepth_unbounded(xyz)

            # Record metadata
            burial_dict['res_n'].append(res_n)
            burial_dict['res_pdb_n'].append(
                pose.pdb_info().number(res_n)
            )
            burial_dict['res_name'].append(res.name3())
            burial_dict['res_chain'].append(
                pose.pdb_info().chain(res_n)
            )
            burial_dict['res_n_nbrs'].append(res_n_nbrs)
            burial_dict['atom_n'].append(atom_n)
            burial_dict['atom_name'].append(
                res.atom_name(atom_n).strip()
            )
            burial_dict['atom_type_name'].append(
                atom_type.name().strip()
            )
            burial_dict['atom_bb'].append(
                res.atom_is_backbone(atom_n)
            )
            burial_dict['atom_lj_radius'].append(atom_radius)
            burial_dict['atom_sasa'].append(
                sasa_map[res_n][atom_n]
            )
            burial_dict['atom_depth'].append(depth)
            burial_dict['atom_x'].append(xyz[0])
            burial_dict['atom_y'].append(xyz[1])
            burial_dict['atom_z'].append(xyz[2])

    burial_df = pandas.DataFrame(burial_dict)
    return burial_df

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", help="the path to a directory with input PDBs")
    parser.add_argument("--output_file", help="the path to an output file")
    parser.add_argument("--probe_radius", help="the radius of the probe used to trace the molecular surface", type=float)
    parser.add_argument("--resolution", help="the resolution of the grid used to compute burial", type=float, default=0.5)
    parser.add_argument("--poly_x", help="will mutate pose to this amino acid for computing depth", type=str, default='L')
    args = parser.parse_args()
    
    # Cycle over PDBs, computing the burial one at a time
    dfs = []
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    pdbs += glob.glob(os.path.join(args.pdb_dir, '*.pdb.gz'))
    for pdb in pdbs:
        df = compute_burial(pdb, args.probe_radius, args.resolution, args.poly_x.upper())
        df['pdb'] = pdb
        dfs.append(df)

    # Write an output file with data from each PDB
    df = pandas.concat(dfs)
    df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
