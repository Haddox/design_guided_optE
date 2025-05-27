* `relax.txt`: template for running the relax protocol
* `cartesian_ddg.txt`: template for running the Cartesian ddG protocol
* `final_scores_250417.csv`: a CSV with ddG estimates for each mutation made by each energy function. Columns report the following values:
	* `pdb`, `chain`, and `mutation`: give the PDB, protein chain, and mutation in question
	* `beta16_cart`: gives the ddG estimated using the beta_nov16 energy function
	* `beta16_cart_plus_hpsc_lj_changes`: gives the ddG estimated using the beta_jan25 energy function
	* other columns give metadata from Frenz et al. (see Methods of manuscript)
* `balanced_ids_final_scores_250417l.csv`: same as above, but only for the balanced set of 735 mutations from Frenz et al.
