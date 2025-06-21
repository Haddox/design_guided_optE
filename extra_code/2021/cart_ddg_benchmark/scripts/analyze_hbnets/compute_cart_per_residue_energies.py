"""Compute per-residue energies in Cartesian space for all score terms and for all input PDBs"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas

# Custom `Python` modules
sys.path.append(os.path.join(os.getcwd(), 'pymodules/'))
import hbond_cart_utils

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", help="a directory with PDBs for analysis")
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    args = parser.parse_args()

    # Get paths to PDBs to analyze
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    print("Found {0} PDBs to analyze".format(len(pdbs)))

    # For each PDB, compute per-residue energies for each score term
    dfs = []
    for (pdb_i, pdb) in enumerate(pdbs, 1):
        print("Analyzing PDB {0}".format(pdb_i))
        df = hbond_cart_utils.compute_per_residue_energies(pdb)
        df['pdb'] = os.path.basename(pdb).replace('.pdb', '')
        dfs.append(df)

    # Concatenate all scores into a single dataframe, with a column called
    # `pdb` that identifies the relevant PDB
    scores_df = pandas.concat(dfs)

    # Write results to output file
    print("Writing results to the output file: {0}".format(args.output_file))
    scores_df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
