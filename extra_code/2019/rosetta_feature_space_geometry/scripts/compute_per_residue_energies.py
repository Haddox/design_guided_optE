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

def get_sf_weights(sf):

    score_types = list(sf.get_nonzero_weighted_scoretypes())
    return {
        str(score_term).replace('ScoreType.', '') : \
            sf.weights().get(score_term)
        for score_term in score_types
    }


def compute_per_residue_energies(pdb, sf):
    """
    Compute per-residue energies for each term in the score function

    Args:
        *pdb*: the path to an input PDB file
        *sf*: PyRosetta score function

    Returns:
        A dataframe with columns giving energies and rows giving
            residues
    """

    # Make a list of all score terms in the energy function
    score_types = list(sf.get_nonzero_weighted_scoretypes())
    score_terms = [
            str(score_term).replace('ScoreType.', '')
            for score_term in score_types
        ] + ['total_score']

    # Read in and score pose
    pose = pyrosetta.pose_from_pdb(pdb)
    sf(pose)

    # Get weights for each score term
    weights = get_sf_weights(sf)
    weights['total_score'] = 1.0

    # Make a dataframe with per-residue weighted scores
    scores_dict = {
        key : []
        for key in ['res_n', 'res_aa'] + score_terms
    }
    for res_n in list(range(1, pose.size()+1)):
        scores_dict['res_n'].append(res_n)
        scores_dict['res_aa'].append(pose.residue(res_n).name1())
        for score_term in score_terms:
            scores_dict[score_term].append(
                weights[score_term] * pose.energies().residue_total_energies(res_n)[
                    pyrosetta.rosetta.core.scoring.score_type_from_name(score_term)
                ]
            )

    scores_df = pandas.DataFrame(scores_dict)
    return scores_df

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", help="a directory with PDBs for analysis")
    parser.add_argument("--flags_f", help="the path to an input flags file")
    parser.add_argument("--weights_f", help="the path to an input weights file")
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    args = parser.parse_args()
    
    
    # Initialize `PyRosetta` with flags from the input flags file
    flags = [
        '-mute all', '-mute core', '-mute protocols',
        f'@{args.flags_f}'
    ]
    pyrosetta.init(extra_options=' '.join(flags))

    # Create a score function with the input weights file
    sf = pyrosetta.create_score_function(args.weights_f)
    def fix_scorefxn(sfxn, allow_double_bb=False):
        opts = sfxn.energy_method_options()
        opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
        opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
        sfxn.set_energy_method_options(opts)
    fix_scorefxn(sf)

    # Get paths to PDBs to analyze
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    print("Found {0} PDBs to analyze".format(len(pdbs)))

    # For each PDB, compute per-residue energies for each score term
    dfs = []
    for (pdb_i, pdb) in enumerate(pdbs, 1):
        print("Analyzing PDB {0}".format(pdb_i))
        df = compute_per_residue_energies(pdb, sf)
        df['pdb'] = os.path.basename(pdb).replace('.pdb', '')
        dfs.append(df)

    # Concatenate all scores into a single dataframe, with a column
    # called `pdb` that identifies the relevant PDB
    scores_df = pandas.concat(dfs)

    # Write results to output file
    print("Writing results to the output file: {0}".format(args.output_file))
    scores_df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
