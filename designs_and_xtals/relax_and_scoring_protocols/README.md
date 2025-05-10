* `relax_protocol.txt`: protocol for relaxing structures using a given energy function.
* `scoring_protocol.txt`: protocol for computing Rosetta energies of structures using a given energy function. 
* `make_poly_gly.py`: code for using PyRosetta to convert a structure to a poly-glycine sequence while keeping the backbone atoms fixed in place.
* `all_relaxed_pdbs/`: PDBs of all structures relaxed with each energy function. The subdirectory called `beta16_cart/` contains structures relaxed with beta_nov16, while the subdirectory called `beta16_cart_plus_hpsc_lj_changes/` contains structures relaxed with beta_jan25. For each input structure (from the `../mod/` directory), the relax protocol generated between 6-10 replicate output structures from independent relax runs (replicate number given in the file name's suffix). We only analyzed replicates 1-6. The above subdirectories also contain poly-glycine versions of each relaxed structure. 
* `all_scores.csv`: a file with Rosetta energies of each relax replicate for each design and each crystal structure. Each row in the file corresponds to a single relax replicate of a single structure. Columns report:
	* `orig_design_name`: the name of the design-crystal pair
	* `design_or_xtal`: indicates whether the relax was performed on the design model or crystal structure
	* `e_function`: indicates the energy function used in relax and scoring. `beta16_cart` corresponds to beta_nov16, while `beta16_cart_plus_hpsc_lj_changes` corresponds to beta_jan25.
	* columns giving groups of energy terms plotted in the paper: `total energy`, `electrostatics`, `H-bond`, `sidechain rotamer`, `LJ`, `solvation`, `backbone torsion`, `covalent bonding`
	* other columns include individual energy terms. Columns with a `_poly_gly` suffix indicate the energy was computed from a poly-glycine version of the structure. Columns with a `_corrected` suffix show the difference between the energy computed for the full structure minus the difference computed for the poly-glycine structure.
* `all_egaps.csv`: a file with differences in Rosetta energies between each design-crystal pair.
        * `orig_design_name`: the name of the design-crystal pair
        * `e_function`: indicates the energy function used in relax and scoring. `beta16_cart` corresponds to beta_n
ov16, while `beta16_cart_plus_hpsc_lj_changes` corresponds to beta_jan25. 
	* `score_term`: the name of the energy term or group of energy terms. Names match the scores dataframe from above.
	* `design` and `xtal`: the value of the energy term computed for the design model or crystal structure, repectively.
	* `egap`: the difference in energy between the design model and crystal structure (=crystal-design)
