"""Compute summary statistics of per-residue energies for all input PDBs"""

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
    parser.add_argument("--designs_and_res_input_file", help="a CSV file with input PDBs and residues to consider")
    parser.add_argument("--energies_input_file", help="a CSV file with per-residue energies")
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    args = parser.parse_args()

    # Read in list of PDBs and sites to analyze
    data = pandas.read_csv(args.designs_and_res_input_file)
    data.rename(columns={'chip_name':'pdb'}, inplace=True)
    assert len(data['pdb']) == len(set(data['pdb']))
    all_designs = set(data['pdb'])
    print("Read in {0} designs to analyze".format(len(all_designs)))

    # Read in per-residue energies, removing _0001 from file names that got
    # added from relax steps
    energies_df = pandas.read_csv(args.energies_input_file)
    del energies_df['Unnamed: 0']
    energies_df['pdb'] = energies_df['pdb'].apply(
        lambda x: x.replace('_0001', '')
    )
    energies_df = energies_df[
        energies_df['pdb'].isin(all_designs)
    ]
    assert len(set(energies_df['pdb'])) == len(all_designs)
    print("Read in energies")

    # Compute summary statistics for individual score terms
    for score_term in energies_df.columns.values:
        print("Computing summary statistics for {0}".format(score_term))
        if score_term in ['res_n', 'res_aa', 'pdb']:
            continue

        # Cycle over PDBs, computing summary statistics for
        # sites of interest
        min_vals = []
        max_vals = []
        mean_vals = []
        total_vals = []
        for (i, row) in data.iterrows():

            # Get a list of network residues, excluding the backbone
            network_residues = row['network_residues']
            if isinstance(network_residues, int):
                network_residues = [network_residues]
            elif isinstance(network_residues, str):
                network_residues = row['network_residues'].split()
                network_residues = [int(res) for res in network_residues]

            # Compute summary statistics for network residues
            data_i = energies_df[
                (energies_df['pdb'] == row['pdb']) &
                (energies_df['res_n'].isin(network_residues))
            ]
            vals = data_i[score_term]
            if len(vals) == 0:
                min_vals.append(np.nan)
                max_vals.append(np.nan)
                mean_vals.append(np.nan)
                total_vals.append(np.nan)
            else:
                min_vals.append(min(vals))
                max_vals.append(max(vals))
                mean_vals.append(np.mean(vals))
                total_vals.append(sum(vals))

        # Add data to dataframe
        data['min_{0}'.format(score_term)] = min_vals
        data['max_{0}'.format(score_term)] = max_vals
        data['mean_{0}'.format(score_term)] = mean_vals
        data['total_{0}'.format(score_term)] = total_vals

    print("Writing results to the output file: {0}".format(args.output_file))
    data.to_csv(args.output_file, index=False)

if __name__ == '__main__':
    main()
