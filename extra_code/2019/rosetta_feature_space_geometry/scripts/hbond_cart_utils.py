"""
This script contains various modules for analyzing hydrogen bonds using PyRosetta

Hugh Haddox, March 27, 2019
"""

# Import `Python` modules
import os
import pandas
import numpy as np
from Bio.Alphabet import IUPAC

# Import and initialize `PyRosetta`
import pyrosetta
import pyrosetta.rosetta
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
pyrosetta.init('-beta_nov16_cart -corrections:beta_nov16 -mute all -mute core -mute protocols')

# Set up a score function in `PyRosetta`
sf = pyrosetta.create_score_function('beta_nov16_cart.wts')
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(sf)

def get_sf_weights():

    score_types = list(sf.get_nonzero_weighted_scoretypes())
    return {
        str(score_term).replace('ScoreType.', '') : sf.weights().get(score_term)
        for score_term in score_types
    }


def compute_per_residue_energies(pdb):
    """
    Compute per-residue energies for each term in the score function

    Args:
        *pdb*: the path to an input PDB file

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
    weights = get_sf_weights()
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

def compute_total_score(pdb):
    pose = pyrosetta.pose_from_pdb(pdb)
    sf(pose)
    return pyrosetta.rosetta.core.pose.total_energy_from_pose(pose)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
