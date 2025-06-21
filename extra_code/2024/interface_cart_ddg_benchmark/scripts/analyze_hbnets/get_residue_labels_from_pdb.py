"""Get residue labels from input PDB files"""

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

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", help="a directory with PDBs for analysis")
    parser.add_argument(
        "--output_file",
        help="the path to an output file with the results of this script"
    )
    args = parser.parse_args()
    
    # Get paths to PDBs to analyze
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    print("Found {0} PDBs to analyze".format(len(pdbs)))

    # Analyze each PDB
    dfs = []
    for (pdb_i, pdb) in enumerate(pdbs, 1):

        print("Analyzing PDB {0}".format(pdb_i))
        df = hbond_utils.get_res_labels_from_pdb(pdb)
        df['pdb'] = pdb
        dfs.append(df)

    # Write an output file with the counts of each sat type
    df = pandas.concat(dfs, sort=False)
    print("Writing results to the output file: {0}".format(args.output_file))
    df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
