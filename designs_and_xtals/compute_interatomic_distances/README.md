* `compute_energies_between_atom_pairs.py`: Python script for computing pairwise interatomic distances in a given input structure. It only returns distances of atom pairs below a given distance cutoff.
* `compute_energies_between_atom_pairs.sh`: Command for calling the above python script. I used weights and flags files for the beta_nov16 energy function.
* `interatomic_distances/`: a directory with interatomic distances computed for a given structure
	* `*_design_mod.csv` and `*_xtal_mod.csv`: files computed for design models and crystal structures in `../mod/`
	* `*_design_mod_0001.csv`: files computed for relaxed design models and crystal structures in `../relax_and_scoring_protocols/all_relaxed_pdbs/beta16_cart/`
	* in the paper, when we compare clashing in design models vs. crystal structures, we use interatomic distances of relaxed designs vs. crystal structures that have not been relaxed.
* `distance_distribution_plots/`: plots showing distance distributions for C:C and C:OCbb for each design and crystal pair.
