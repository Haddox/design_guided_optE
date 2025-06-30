"""Use Needle to compute the identity among all pairs of input sequences"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas

# Custom `Python` modules
sys.path.append(os.path.join(os.getcwd(), 'scripts/'))
import design_utils

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", help="a FASTA file with all input sequences")
    parser.add_argument("--needle_output_file", help="the path to an output file with the results of Needle")
    parser.add_argument("--percent_id_output_file", help="the path to an output file with a dataframe of all pairwise sequence identities")
    args = parser.parse_args()

    # Use needle to compute the sequence identity between all pairs of sequences
    cmd = ' '.join([
        'needleall',
        args.fasta_file,
        args.fasta_file,
        '-gapopen 10',
        '-gapextend 0.5',
        '-outfile {0}'.format(args.needle_output_file),
        '-aformat pair'
    ])
    print("Aligning sequences with Needle using the command:\n{0}".format(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    print("Needle err:")
    print(err)

    # Parse the Needle output, storing all pairwise identities in a dataframe
    print("Parsing the Needle output file")
    df = design_utils.parse_needle_output_file(args.needle_output_file)
    df.to_csv(args.percent_id_output_file, index=False)

if __name__ == '__main__':
    main()
