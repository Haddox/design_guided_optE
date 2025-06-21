"""Compute per-residue energies in Cartesian space for all score terms and for all input PDBs"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas
import pyrosetta
import pyrosetta.rosetta

# Define functions
def compute_pairwise_energies(
    pose, res_i_n, atom_i_n, sf, ij_and_ji=False, max_dist=5,
    max_dist_minus_sigma=0.5
):
    """
    Compute pairwise energies between an input atom and all others
    
    Args:
        *pose*: a PyRosetta pose
        *res_i_n*: the number of the input residue
        *atom_i_n*: the number of the input atom on the input residue
        *sf*: the score function used to compute energies
        *ij_and_ji*: a boolean specifying whether to compute all ij and
            ji interactions. If False (default), then only return
            interactions when the residue number of i is less than j 
        *max_dist*: will only return data for atom pairs within this
            distance of each other
        *max_dist_minus_sigma*: will only return data for atom pairs
            with d-o values less than or equal to this value, where
            o is the sum of the atomic radii
        
    Returns:
        A dataframe where rows are pairwise interactions, and columns
            give metadata on the pair of atoms involved and their
            energies
    """
    
    # Initiate a dictionary to record data
    energies_dict = {
        key : []
        for key in [
            'res_j_n', 'res_j_pdb_n', 'res_j_name', 'res_j_ss',
            'atom_j_n', 'atom_j_name', 'atom_j_type_name',
            'lj_atr', 'lj_rep', 'fa_solv', 'fa_elec',
            'atom_j_lj_radius', 'atom_j_lj_wdepth',
            'd'
        ]
    }
    
    # Get the SS of each site in the protein
    DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
    DSSP.apply(pose)
    ss = pose.secstruct()
    
    # Make an object corresponding to the input residue and get
    # the XYZ coordinates of the input atom
    res_i = pose.residue(res_i_n)
    atom_i_xyz = res_i.xyz(atom_i_n)
    
    # Record pairwise energies between the input atom and all other
    # atoms in the structure
    etable_atom_pair_energies = \
        pyrosetta.toolbox.atom_pair_energy.etable_atom_pair_energies
    for res_j_n in range(1, pose.size()+1):
        if not ij_and_ji:
            if res_j_n < res_i_n:
                continue
        res_j = pose.residue(res_j_n)
        for atom_j_n in range(1, res_j.natoms()+1):
            
            # Skip res_j and atom_j if they are the same as the
            # input residue and atom
            if (res_j_n == res_i_n) and (atom_j_n == atom_i_n):
                continue

            # Compute pairwise energies between atoms i and j
            (lj_atr, lj_rep, fa_solv, fa_elec) = \
                etable_atom_pair_energies(
                    res1=res_i,
                    atom_index_1=atom_i_n,
                    res2=res_j,
                    atom_index_2=atom_j_n,
                    sfxn=sf
                )
            
            # Record data in dictionary
            atom_type_j = res_j.atom_type(atom_j_n)
            atom_j_xyz = res_j.xyz(atom_j_n)
            energies_dict['res_j_n'].append(res_j_n)
            energies_dict['res_j_pdb_n'].append(
                pose.pdb_info().number(res_j_n)
            )
            energies_dict['res_j_name'].append(res_j.name3())
            energies_dict['res_j_ss'].append(ss[res_j_n-1])
            energies_dict['atom_j_n'].append(atom_j_n)
            energies_dict['atom_j_name'].append(
                res_j.atom_name(atom_j_n).strip()
            )
            energies_dict['atom_j_type_name'].append(
                atom_type_j.name().strip()
            )
            energies_dict['lj_atr'].append(lj_atr)
            energies_dict['lj_rep'].append(lj_rep)
            energies_dict['fa_solv'].append(fa_solv)
            energies_dict['fa_elec'].append(fa_elec)
            energies_dict['atom_j_lj_radius'].append(
                atom_type_j.lj_radius()
            )
            energies_dict['atom_j_lj_wdepth'].append(
                atom_type_j.lj_wdepth()
            )
            energies_dict['d'].append(
                (atom_i_xyz - atom_j_xyz).norm()
            )
           
    # Convert the above dictionary into a dataframe and then remove
    # entries where the atoms pairs are separated by a distance
    # that is greater than the `max_dist` variable
    energies_df = pandas.DataFrame(energies_dict)
    energies_df = energies_df[energies_df['d'] <= max_dist]
    
    # Add information on the input residue and atom
    energies_df['res_i_n'] = res_i_n
    energies_df['res_i_pdb_n'] = pose.pdb_info().number(res_i_n)
    energies_df['res_i_name'] = res_i.name3()
    energies_df['res_i_ss'] = ss[res_i_n-1]
    energies_df['atom_i_n'] = atom_i_n
    energies_df['atom_i_name'] = res_i.atom_name(atom_i_n).strip()
    atom_type_i = res_i.atom_type(atom_i_n)
    energies_df['atom_i_type_name'] = atom_type_i.name().strip()
    energies_df['atom_i_lj_radius'] = atom_type_i.lj_radius()
    energies_df['atom_i_lj_wdepth'] = atom_type_i.lj_wdepth()

    # Compute values that characterize VdW interactions, including
    # the sum of atomic radii of atoms i and j (o), and remove entries
    # where d-o is greater than the specified cutoff
    energies_df['o'] = energies_df.apply(
        lambda row: row['atom_i_lj_radius'] + row['atom_j_lj_radius'],
        axis=1
    )
    energies_df['d-o'] = energies_df['d'] - energies_df['o']
    energies_df = energies_df[energies_df['d-o'] <= max_dist_minus_sigma]
    
    #... the well depth, computed as the geometric mean of the wdepth
    # values for atoms i and j...
    energies_df['e'] = energies_df.apply(
        lambda row: math.sqrt(
            row['atom_i_lj_wdepth'] * row['atom_j_lj_wdepth']
        ),
        axis=1
    )
    
    # ... the LJ potential for atoms i and j...
    energies_df['lj_energy'] = energies_df.apply(
        lambda row: compute_lj_energy(row['e'], row['o'], row['d']),
        axis=1
    )
    
    # ... the LJ potential computed using Rosetta...
    energies_df['rosetta_lj_energy'] = energies_df.apply(
        lambda row: compute_rosetta_lj_energy(
            row['e'], row['o'], row['d'], rep_weight=1.0
        ),
        axis=1
    )
    
    # ... the LJ potential computed using Rosetta with a weight of
    # 0.55 on the repulsive term
    energies_df['weighted_rosetta_lj_energy'] = energies_df.apply(
        lambda row: compute_rosetta_lj_energy(
            row['e'], row['o'], row['d'], rep_weight=0.55
        ),
        axis=1
    )
    
    # Return the dataframe
    return energies_df


def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_file", help="the path to a PDB file")
    parser.add_argument("--weights_file", help="the path to a weights file")
    parser.add_argument("--flags_file", help="the path to a flags file")
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    parser.add_argument("--only_analyze_chain_A", help="only analyze chain A", action="store_true")
    args = parser.parse_args()

    # Set up a score function in `PyRosetta`
    pyrosetta.init(f'-mute all -mute core -mute protocols @{args.flags_file}')
    sf = pyrosetta.create_score_function(args.weights_file)
    def fix_scorefxn(sfxn, allow_double_bb=False):
        opts = sfxn.energy_method_options()
        opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
        opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
        sfxn.set_energy_method_options(opts)
    fix_scorefxn(sf)
    
    # Read in the pose and, if specified, subsetting the pose to just
    # chain A
    if args.only_analyze_chain_A:
        initial_pose = pyrosetta.pose_from_pdb(args.pdb_file)
        chains = initial_pose.split_by_chain()
        found_chain_a = False
        for chain in chains:
            chain_name = chain.pdb_info().chain(1)
            if chain_name == 'A':
                pose = chain
                found_chain_a = True
        if not found_chain_a:
            raise ValueError("Could not find chain A")
    else:
        pose = pyrosetta.pose_from_pdb(args.pdb_file)
            
    # Iterate over all residues and compute atom-atom energies for all
    # pairs of atoms from that residue
    dfs = []
    for res_i_n in range(1, pose.size()+1):
        res_i = pose.residue(res_i_n)
        res_i_name = res_i.name3()
        for atom_i_n in list(range(1, res_i.natoms()+1)):
            atom_i_name = res_i.atom_name(atom_i_n)
            df = compute_pairwise_energies(
                pose, res_i_n, atom_i_n, sf
            )
            dfs.append(df)
            
    # Make a single dataframe with all atom pair energies
    energies_df = pandas.concat(dfs)
    energies_df['pdb_file'] = args.pdb_file
    energies_df['pdb'] = os.path.basename(args.pdb_file)
    energies_df['pair_id'] = energies_df.apply(
        lambda row: '_'.join(sorted([
            row['pdb'],
            "{0}_{1}".format(row['res_i_n'], row['atom_i_n']),
            "{0}_{1}".format(row['res_j_n'], row['atom_j_n']),
        ])),
        axis=1
    )
    energies_df.drop_duplicates('pair_id', inplace=True, keep='first')
    energies_df = energies_df[
        (energies_df['res_i_n'] != energies_df['res_j_n'])
    ]
    print("Writing results to the output file: {0}".format(args.output_file))
    energies_df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
