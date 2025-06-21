"""Compute summary statistics of H-bond energies, distances, and angles"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas
import numpy as np

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--energies_input_file", help="a CSV file with H-bond energies")
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    args = parser.parse_args()

    # Read in per-residue energies
    hbond_df = pandas.read_csv(args.energies_input_file)
    print("Read in energies")

    # Add columns that give the absolute value of the difference between
    # the BAH angles and 120 or 109.5 degrees...
    hbond_df['BAHangle_abs_diff_from_120'] = hbond_df['BAHangle'].apply(
        lambda x: abs(x - 120)
    )
    hbond_df['BAHangle_amnt_below_120'] = hbond_df['BAHangle'].apply(
        lambda x: max([120 - x, 0])
    )
    hbond_df['BAHangle_amnt_above_120'] = hbond_df['BAHangle'].apply(
        lambda x: max([x - 120, 0])
    )

    hbond_df['BAHangle_abs_diff_from_109_5'] = hbond_df['BAHangle'].apply(
        lambda x: abs(x - 109.5)
    )
    hbond_df['BAHangle_amnt_below_109_5'] = hbond_df['BAHangle'].apply(
        lambda x: max([109.5 - x, 0])
    )
    hbond_df['BAHangle_amnt_above_109_5'] = hbond_df['BAHangle'].apply(
        lambda x: max([x - 109.5, 0])
    )

    # ... and the minimum absolute value of the difference between the BAtorsion
    # angle and each of the angles specified below
    hbond_df['BAtorsion_abs_diff_from_planarity'] = hbond_df['BAtorsion'].apply(
        lambda x: min([abs(180 - abs(x)), abs(x)])
    )
    hbond_df['BAtorsion_abs_diff_from_zero'] = hbond_df['BAtorsion'].apply(
        lambda x: abs(x)
    )
    hbond_df['BAtorsion_abs_diff_from_180'] = hbond_df['BAtorsion'].apply(
        lambda x: abs(180 - abs(x))
    )
    hbond_df['BAtorsion_abs_diff_from_135'] = hbond_df['BAtorsion'].apply(
        lambda x: abs(135 - abs(x))
    )
    hbond_df['BAtorsion_abs_diff_from_90'] = hbond_df['BAtorsion'].apply(
        lambda x: abs(90 - abs(x))
    )

    # Make a list of columns to related to H-bond energies to
    # compute stats across all residues
    energy_metrics = [
        'energy', 'weight', 'lj_atr', 'lj_rep', 'fa_solv', 'fa_elec',
        'energy_hbond_plus_elec', 'weighted_energy_hbond_plus_elec',
        'HAdist', 'AHDangle'
    ]
    energy_stats = []
    for metric in energy_metrics:
        energy_stats += [
            '{0}_{1}'.format(stat, metric)
            for stat in ['min', 'max', 'mean', 'total']
        ]

    # Make a list of columns related to specific distances, angles,
    # and hybridization types, for computing stats across a subset
    # of relevant residues
    hybrid_stats = []
    hybrid_metrics = [
        'energy', 'weight', 'lj_atr', 'lj_rep', 'fa_solv', 'fa_elec',
        'energy_hbond_plus_elec', 'weighted_energy_hbond_plus_elec',
        'HAdist', 'AHDangle', 'BAHangle', 'BAtorsion',
        'BAHangle_abs_diff_from_120', 'BAHangle_amnt_below_120',
        'BAHangle_amnt_above_120', 'BAHangle_abs_diff_from_109_5',
        'BAHangle_amnt_below_109_5', 'BAHangle_amnt_above_109_5',
        'BAtorsion_abs_diff_from_planarity', 'BAtorsion_abs_diff_from_zero',
        'BAtorsion_abs_diff_from_180', 'BAtorsion_abs_diff_from_135',
        'BAtorsion_abs_diff_from_90'
    ]
    hybridizations = list(set(hbond_df['acc_hybridization']))
    for h in hybridizations:
        for metric in hybrid_metrics:
            hybrid_stats += [
                '{0}_{1}_{2}'.format(stat, metric, h)
                for stat in ['min', 'max', 'mean']
            ]

    # Compute stats
    summary_stats_dict = {
        key : []
        for key in ['pdb', 'n_hbonds_analyzed'] + energy_stats + hybrid_stats
    }
    pdbs = list(set(hbond_df['pdb']))
    for pdb in pdbs:

        # Get data for the given PDB
        summary_stats_dict['pdb'].append(pdb)
        data = hbond_df[hbond_df['pdb'] == pdb].copy()

        # Record the number of H-bonds being analyzed
        summary_stats_dict['n_hbonds_analyzed'].append(len(data.index.values))

        # Record summary stats for H-bond energies and weights
        for col in energy_metrics:
            vals = data[col]
            summary_stats_dict['min_{0}'.format(col)].append(min(vals))
            summary_stats_dict['max_{0}'.format(col)].append(max(vals))
            summary_stats_dict['mean_{0}'.format(col)].append(np.mean(vals))
            summary_stats_dict['total_{0}'.format(col)].append(sum(vals))

        # Record stats that are specific for given hybridization types
        for col in hybrid_metrics:
            for h in hybridizations:
                data_i = data[
                    data['acc_hybridization'] == h
                ].copy()
                if len(data_i) == 0:
                    summary_stats_dict['min_{0}_{1}'.format(col, h)].append(
                        np.nan
                    )
                    summary_stats_dict['max_{0}_{1}'.format(col, h)].append(
                        np.nan
                    )
                    summary_stats_dict['mean_{0}_{1}'.format(col, h)].append(
                        np.nan
                    )
                else:
                    vals = data_i[col]
                    summary_stats_dict['min_{0}_{1}'.format(col, h)].append(
                        min(vals)
                    )
                    summary_stats_dict['max_{0}_{1}'.format(col, h)].append(
                        max(vals)
                    )
                    summary_stats_dict['mean_{0}_{1}'.format(col, h)].append(
                        np.mean(vals)
                    )

    summary_stats_df = pandas.DataFrame(summary_stats_dict)

    print("Writing results to the output file: {0}".format(args.output_file))
    summary_stats_df.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
