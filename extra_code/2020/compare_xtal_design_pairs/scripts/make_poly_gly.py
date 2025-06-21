"""Make a polyglycine backbone from an input pose"""

# Import `Python` modules
import os
import sys
import argparse
sys.path.append(os.path.join(os.getcwd(), 'scripts/'))
import design_utils

# Main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pdb", help="the path to an input PDB")
    parser.add_argument("--output_pdb", help="the path to an output PDB with the poly-glycine pose")
    args = parser.parse_args()
    
    # Make poly-glycine backbone from input PDB
    design_utils.make_poly_gly(args.input_pdb, args.output_pdb)

if __name__ == '__main__':
    main()
