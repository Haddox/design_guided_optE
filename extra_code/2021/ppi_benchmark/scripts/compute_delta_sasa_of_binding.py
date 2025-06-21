"""Identify decoys with chain breaks at an interface between two chains"""

# Import `Python` modules
import os
import argparse
import subprocess
import pandas
import pyrosetta
import pyrosetta.rosetta

# Initialize PyRosetta
flags = [
    '-mute all', '-mute core', '-mute protocols',
    '-ignore_unrecognized_res',
    '-read_only_ATOM_entries',
    '-in:file:silent_struct_type binary',
    '-crystal_refine true',
    '-silent_read_through_errors'
]
pyrosetta.init(extra_options=' '.join(flags))

# Set up a score function in `PyRosetta`
sf = pyrosetta.get_fa_scorefxn()
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(sf)

def pull_apart_chains(pose):

    # Copy pose
    split_pose = pose.clone()

    # Get an arbitrary translation vector that is
    # very large in magnitude
    res = split_pose.residue(1)
    xyz_i = res.xyz(1)
    xyz_j = res.xyz(2)
    dv = xyz_j - xyz_i
    translation_vector = dv.normalize(5000)

    # Translate all atoms that aren't in chain A and then
    # return the split pose
    for res_n in list(range(1, split_pose.size()+1)):
        res = split_pose.residue(res_n)
        chain = split_pose.pdb_info().chain(res_n)
        if chain == 'A':
            continue
        for atom_n in list(range(1, res.natoms()+1)):
            res.set_xyz(
                atom_n,
                res.xyz(atom_n) + translation_vector
            )

    return split_pose

def compute_delta_interface_sasa(input_silent_file, sasa_probe_radius):

    # Loop over each pose and record per-atom SASAs
    sasa_dict = {
        key : []
        for key in [
            'tag',
            'res_n', 'res_name', 'atom_n', 'atom_name', 'atom_type_name',
            'atom_x_bound', 'atom_x_unbound',
            'sasa_bound', 'sasa_unbound',
            'atom_is_donor', 'atom_is_acceptor', 'atom_is_polar_hydrogen',
        ]
    }
    sfd = pyrosetta.rosetta.core.io.silent.SilentFileData(
        pyrosetta.rosetta.core.io.silent.SilentFileOptions()
    )
    sfd.read_file(input_silent_file)
    for tag in sfd.tags():

        # Read in pose and compute SASA for all atoms
        ss = sfd.get_structure(tag)
        pose = pyrosetta.Pose()
        ss.fill_pose(pose)
        sf(pose)
        #bound_sasa_calc = \
            #pyrosetta.rosetta.protocols.vardist_solaccess.VarSolDistSasaCalculator()
        bound_sasa_calc = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
        bound_sasa_calc.set_probe_radius(sasa_probe_radius)
        bound_sasa_calc.calculate(pose)
        bound_sasa_map = bound_sasa_calc.get_atom_sasa()

        # Next, pull apart the chains and compute SASA in the
        # unbound state
        split_pose = pull_apart_chains(pose)
        sf(split_pose)
        #unbound_sasa_calc = \
        #    pyrosetta.rosetta.protocols.vardist_solaccess.VarSolDistSasaCalculator()
        unbound_sasa_calc = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
        unbound_sasa_calc.set_probe_radius(sasa_probe_radius)
        unbound_sasa_calc.calculate(split_pose)
        unbound_sasa_map = unbound_sasa_calc.get_atom_sasa()

        # Record data for each atom
        for res_n in list(range(1, pose.size()+1)):
            res = pose.residue(res_n)
            for atom_n in list(range(1, res.natoms()+1)):

                atom_type = res.atom_type(atom_n)
                atom_name = res.atom_name(atom_n).strip()
                atom_type_name = atom_type.atom_type_name()

                atom_id = pyrosetta.AtomID(atom_n, res_n)
                #atom_hbonds = hbond_set.atom_hbonds(
                #    atom_id, include_only_allowed=False
                #)
                sasa_bound = bound_sasa_map[res_n][atom_n]
                sasa_unbound = unbound_sasa_map[res_n][atom_n]

                sasa_dict['tag'].append(os.path.basename(tag))
                sasa_dict['res_n'].append(res_n)
                sasa_dict['res_name'].append(res.name1())
                sasa_dict['atom_n'].append(atom_n)
                sasa_dict['atom_name'].append(atom_name)
                sasa_dict['atom_type_name'].append(atom_type_name)
                sasa_dict['atom_x_bound'].append(res.xyz(atom_n).x)
                sasa_dict['atom_x_unbound'].append(split_pose.residue(res_n).xyz(atom_n).x)
                sasa_dict['sasa_bound'].append(sasa_bound)
                sasa_dict['sasa_unbound'].append(sasa_unbound)
                sasa_dict['atom_is_donor'].append(atom_type.is_donor())
                sasa_dict['atom_is_acceptor'].append(atom_type.is_acceptor())
                sasa_dict['atom_is_polar_hydrogen'].append(atom_type.is_polar_hydrogen())

    sasa_df = pandas.DataFrame(sasa_dict)
    sasa_df['polar'] = sasa_df[[
        'atom_is_donor', 'atom_is_acceptor',
        'atom_is_polar_hydrogen'
    ]].any(axis=1)
    sasa_df['delta_sasa'] = sasa_df['sasa_bound'] - sasa_df['sasa_unbound']

    # Next, compute the total change in SASA across all atoms
    # upon binding. Do same for just polar or just non-polar
    # atoms.
    delta_sasa_dict = {
        key : []
        for key in ['pdb', 'dsasa_all', 'dsasa_polar', 'dsasa_non_polar']
    }
    tags = list(set(sasa_df['tag']))
    for tag in tags:

        delta_sasa_dict['pdb'].append(tag)

        data = sasa_df[sasa_df['tag'] == tag]
        delta_sasa_dict['dsasa_all'].append(sum(data['delta_sasa']))

        data = sasa_df[
            (sasa_df['tag'] == tag) &
            (sasa_df['polar'] == True)
        ]
        delta_sasa_dict['dsasa_polar'].append(sum(data['delta_sasa']))

        data = sasa_df[
            (sasa_df['tag'] == tag) &
            (sasa_df['polar'] == False)
        ]
        delta_sasa_dict['dsasa_non_polar'].append(sum(data['delta_sasa']))

    delta_sasa_df = pandas.DataFrame(delta_sasa_dict)

    return delta_sasa_df

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_silent_file", help="path to input silent file")
    parser.add_argument("--probe_radius", help="radius used to compute SASA")
    parser.add_argument("--output_file", help="path to output file")
    args = parser.parse_args()

    delta_sasa_df = compute_delta_interface_sasa(args.input_silent_file, float(args.probe_radius))
    delta_sasa_df.to_csv(args.output_file, index=False)
    
if __name__ == '__main__':
    main()
