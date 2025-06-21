"""Compute the energies and other chemical information of hydrogen bonds"""

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
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    parser.add_argument("--only_return_buried_hbonds", type=str2bool, help="a bool indicating whether to only return data for buried H bonds.")
    args = parser.parse_args()

    print("Read in the following input arguments:")
    print("pdb_dir: {0}".format(args.pdb_dir))
    print("output_file: {0}".format(args.output_file))
    print("only_return_buried_hbonds: {0}".format(args.only_return_buried_hbonds))

    # Get paths to PDBs to analyze
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    print("Found {0} PDBs to analyze".format(len(pdbs)))

    # Analyze each PDB
    dfs = []
    for (pdb_i, pdb) in enumerate(pdbs, 1):
        print("Analyzing PDB {0}".format(pdb_i))
        dfs.append(hbond_utils.compute_hbond_energies(
            pdb,
            only_return_buried_hbonds=args.only_return_buried_hbonds
        ))

    hbond_df = pandas.concat(dfs)
    hbond_df['pdb_basename'] = hbond_df['pdb'].apply(os.path.basename)

    # Write results to output file
    print("Writing results to the output file: {0}".format(args.output_file))
    hbond_df.to_csv(args.output_file)

if __name__ == '__main__':
    main()
