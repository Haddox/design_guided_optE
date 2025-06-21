"""Identify decoys with chain breaks at an interface between two chains"""

# Import `Python` modules
import os
import argparse
import subprocess
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
    parser.add_argument("--input_silent_file", help="path to input silent file")
    parser.add_argument("--ignore_res_file", help="ignore new-line-delimited residues in this file when recording chain breaks")
    parser.add_argument("--ignore_S_", type=str2bool, default=False, help="Ignore tags with _S_ in the name")
    parser.add_argument("--output_prefix", help="a prefix for output files listing passing and failing tags")
    parser.add_argument("--output_silent_file", help="path to silent file with only passing poses")
    args = parser.parse_args()

    # Read in which residues to ignore when recording chain breaks. I am adding this
    # because some decoys contain artificial chain breaks from stiching together
    # two chains that are different in real life.
    if not os.path.isfile(args.ignore_res_file):
        res_to_ignore = []
    else:
        with open(args.ignore_res_file) as f:
            lines = f.readlines()
            res_to_ignore = [int(line.strip()) for line in lines]

    # Loop over each pose in the input silent file and determine
    # if the pose is problematic
    sfd = pyrosetta.rosetta.core.io.silent.SilentFileData(
        pyrosetta.rosetta.core.io.silent.SilentFileOptions()
    )
    sfd.read_file(args.input_silent_file)
    passing_tags = []
    failing_tags = []
    for tag in sfd.tags():
        if args.ignore_S_:
            if '_S_' in tag:
                continue
        ss = sfd.get_structure(tag)
        pose = pyrosetta.Pose()
        ss.fill_pose(pose)
        
        # Loop over each residue but the last and measure the
        # distance between the CObb of that residue and Nbb of
        # the next residue. Identify residues with distances
        # that are greater than 2A.
        res_at_chain_break = set()
        for res_i_n in range(1, pose.size()):

            # Ignore residue if it is in the list to ignore
            if res_i_n in res_to_ignore:
                continue

            # Read in res_i object and don't consider if it is a
            # virtual residue
            res_i = pose.residue(res_i_n)
            if res_i.name() == 'VRT':
                continue

            # Get index for next residue, and make sure that it
            # is not a virtual residue and that both residues
            # are in the same chain
            res_j_n = res_i_n + 1
            res_j = pose.residue(res_j_n)
            if res_j.name() == 'VRT':
                continue
            chain_i = pose.pdb_info().chain(res_i_n)
            chain_j = pose.pdb_info().chain(res_j_n)
            if chain_i != chain_j:
                continue

            # Measure distance between CObb of res_i and Nbb
            # of res_j
            CObb_xyz = res_i.xyz(3)
            Nbb_xyz = res_j.xyz(1)
            d = (CObb_xyz - Nbb_xyz).norm()
            if d > 2:
                res_at_chain_break = set.union(
                    res_at_chain_break,
                    set([res_i_n, res_j_n])
                )

        # Loop over each chain and get a list of residues that make close
        # contacts with that chain. The list will include residues in
        # the query chain
        residue_selector = pyrosetta.rosetta.core.select.residue_selector
        close_contacts = []
        for chain in ['A', 'B']:
            rs_chain = residue_selector.ChainSelector(chain)
            rs_chain_cc = residue_selector.CloseContactResidueSelector()
            rs_chain_cc.central_residue_group_selector(rs_chain)
            rs_chain_cc.threshold(3.0)
            chain_cc = sorted([
                i for (i, bool_i) in enumerate(rs_chain_cc.apply(pose), 1)
                if bool_i == 1
            ])
            close_contacts.append(set(chain_cc))

        # Define interface residues as the intersection of the above
        # two sets
        interface_res = set.intersection(
            close_contacts[0],
            close_contacts[1]
        )

        # Determine if interface residues overlap with chain breaks
        interface_chain_break_res = set.intersection(
            interface_res,
            res_at_chain_break
        )
        if len(interface_chain_break_res) == 0:
            passing_tags.append(tag)
        else:
            failing_tags.append(tag)
            
    # Write two files, one with passing tags and one with failing tags
    passing_tags_f = f'{args.output_prefix}passing_tags.txt'
    with open(passing_tags_f, 'w') as f:
        for tag in passing_tags:
            f.write(f'{tag}\n')
    failing_tags_f = f'{args.output_prefix}failing_tags.txt'
    with open(failing_tags_f, 'w') as f:
        for tag in failing_tags:
            f.write(f'{tag}\n')

    # Make a silent file that only contains poses that passed
    if len(passing_tags) > 0:
        cmd = ' '.join([
            f'cat {passing_tags_f} |',
            f'/home/haddox/software/silent_tools/silentslice {args.input_silent_file} >',
            args.output_silent_file
        ])
        out = subprocess.check_call(cmd, shell=True)

if __name__ == '__main__':
    main()
