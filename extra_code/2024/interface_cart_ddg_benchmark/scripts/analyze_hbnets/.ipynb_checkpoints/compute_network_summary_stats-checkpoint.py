"""Identify and characterize H-bond networks for all input PDBs"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas
import numpy as np
import networkx as nx

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", help="a directory with PDBs for analysis")
    parser.add_argument("--energies_input_file", help="a CSV file with H-bond energies")
    parser.add_argument("--output_file", help="the path to an output file with the results of this script")
    args = parser.parse_args()

    # Get paths to PDBs to analyze
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    print("Found {0} PDBs to analyze".format(len(pdbs)))

    # Read in data on energies of all H-bonds in the protein
    hbond_energies_df = pandas.read_csv(args.energies_input_file)

    # For each PDB, identify all H-bond networks in that PDB
    networks_dict = {
        key : []
        for key in [
            'pdb', 'network', 'n_buried_hbonds', #'avg_dist',
            'diameter'
        ]
    }
    for pdb in pdbs:

        # Get data for the PDB of interest
        data = hbond_energies_df[hbond_energies_df['pdb'] == pdb].copy()

        # Find the max residue number
        max_res_n = max([data['acc_res_n'].max(), data['don_res_n'].max()])
        residues = list(range(0, max_res_n + 1))

        # Initialize an empty adjacency matrix starting at with
        # N + 1 columns, where N is the max residue number from above.
        # The zero row/column is just space filler; row/column 1
        # corresponds to res 1
        aM = np.zeros((max_res_n + 1, max_res_n + 1))

        # Go through all H-bonds in the protein and fill out the
        # adjacency matrix with '1's for all residue pairs making
        # H-bonds, indexed by residue number
        for (i, row) in data.iterrows():
            acc_matrix_n = row['acc_res_n']
            don_matrix_n = row['don_res_n']
            aM[acc_matrix_n, don_matrix_n] = 1
            aM[don_matrix_n, acc_matrix_n] = 1

        # Next, use networkx to extract all networks from the
        # adjacency matrix. For each network, count the number
        # of H-bonds in the network that are buried (based on the
        # VSASA of the donor and acceptor atoms)
        adjacency_df = pandas.DataFrame(aM, columns=residues)
        G = nx.from_pandas_adjacency(adjacency_df)
        for network in nx.connected_components(G):
            if len(network) > 1:
                networks_dict['pdb'].append(pdb)
                networks_dict['network'].append(
                    tuple(sorted(list(network)))
                )
                networks_dict['n_buried_hbonds'].append(sum(data[
                    data['acc_res_n'].isin(network) |
                    data['don_res_n'].isin(network)
                ]['buried']))
                subgraph = G.subgraph(network)
                networks_dict['diameter'].append(
                    nx.diameter(subgraph)
                )

    # Make a dataframe with all networks in each PDB, and write
    # this dataframe to an output file
    network_df = pandas.DataFrame(networks_dict)
    print("Writing results to the output file: {0}".format(args.output_file))
    network_df.to_csv(args.output_file)

if __name__ == '__main__':
    main()
