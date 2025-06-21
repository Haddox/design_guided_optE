"""Compute the degree of satisfaction of polar atoms and residues in input structures"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas

# Custom `Python` modules
sys.path.append(os.path.join(os.getcwd(), 'scripts/'))
import hbond_utils

def str2bool(v):
    """
    Function for parsing boolean arguments from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    """
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
    parser.add_argument("--pdb_dir", help="a directory with PDBs for analysis")
    parser.add_argument("--only_analyze_buried_atoms", type=str2bool, help="a bool. if True, the output will only include data for buried atoms and residues where all polar atoms are buried")
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    args = parser.parse_args()

    # Get paths to PDBs to analyze
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    print("Found {0} PDBs to analyze".format(len(pdbs)))

    # Analyze each PDB
    dfs = []
    available_don_acc_dict = {
        key : []
        for key in ['pdb', 'res_n', 'res_name', 'avail_acc', 'avail_don']
    }
    for (pdb_i, pdb) in enumerate(pdbs, 1):

        print("Analyzing PDB {0}".format(pdb_i))

        # Compute the degree of satisfaction of each polar atom
        # in the structure
        degree_sat_df = hbond_utils.compute_degree_satisfaction(pdb)

        # Compute summary statistics on the percent capacity of
        # network residues
        if args.only_analyze_buried_atoms:
            sc_atoms_df = degree_sat_df[
                (degree_sat_df['bb_or_sc'] == 'sc') &
                (degree_sat_df['buried'] == True)
            ].copy()
        else:
            sc_atoms_df = degree_sat_df[
                (degree_sat_df['bb_or_sc'] == 'sc')
            ].copy()
        if len(sc_atoms_df) == 0:
            continue
        df = hbond_utils.compute_percent_capacity_summary_statistics(sc_atoms_df)
        df['pdb'] = pdb

        # Compute summary statistics about the satisfaction at the
        # level of entire residues. If indicated, only do this for residues
        # where all polar side-chain heavy atoms are buried
        res_sat_dict = {}
        for res_n in set(sc_atoms_df['res_n']):

            # Get data for current residue being considered in the loop
            data = degree_sat_df[
                (degree_sat_df['res_n'] == res_n) &
                (degree_sat_df['bb_or_sc'] == 'sc')
            ]
            if args.only_analyze_buried_atoms:
                if not data['buried'].all():
                    continue

            # Determine whether the residue is available as a donor or acceptor
            res_name = list(set(data['res_name']))
            assert len(res_name) == 1
            res_name = res_name[0]
            (available_acc, available_don) = \
                hbond_utils.determine_if_res_has_available_hbond_don_and_acc(data)
            available_don_acc_dict['pdb'].append(pdb)
            available_don_acc_dict['res_n'].append(res_n)
            available_don_acc_dict['res_name'].append(res_name)
            available_don_acc_dict['avail_acc'].append(available_acc)
            available_don_acc_dict['avail_don'].append(available_don)

            # Compute a string giving the residue's overall level of satisfaction
            # as an acceptor and donor
            res_sat = hbond_utils.determine_residue_satisfaction(data)
            if res_sat not in res_sat_dict.keys():
                res_sat_dict[res_sat] = 1
            else:
                res_sat_dict[res_sat] += 1

        # Merge together all atom- and residue-level data from above
        # and append the resulting dataframe to a growing list
        res_sat_df = pandas.DataFrame.from_dict({
            key : [value]
            for (key, value) in res_sat_dict.items()
        })
        df = df.merge(res_sat_df, right_index=True, left_index=True)
        df.set_index('pdb', inplace=True)
        dfs.append(df)

    # Write an output file with the counts of each sat type
    df = pandas.concat(dfs, sort=False)
    print("Writing results to the output file: {0}".format(args.output_file))
    df.to_csv(args.output_file)

    # Write an output file with info on available donors and acceptors
    available_don_acc_df = pandas.DataFrame(available_don_acc_dict)
    available_don_acc_df.to_csv(
        args.output_file.replace('.csv', '_don_acc_availability.csv'),
        index=False
    )

if __name__ == '__main__':
    main()
