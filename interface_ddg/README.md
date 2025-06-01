* We downloaded the following files from the SKEMPI v2.0 database at `https://life.bsc.es/pid/skempi2/database/index`; 08.06.2018 file versions
	* `skempi_v2.csv`: a CSV with experimentally measured ddG values
	* cleaned PDB files
* `relax.txt`: template for running the relax protocol
* `cartesian_ddg.txt`: template for running the Cartesian ddG protocol 
* `computed_and_experimental_ddg_values.csv`: a CSV with computed and experimentally measured ddG values, where each row reports data for a given mutation to a given structure, with columns reporting:
	* `pdb_id`: the four-letter code of the PDB
	* `mut_pdb`: the mutation, with the first letter reporting the wildtype amino acid, the second letter reporting the protein chain of the mutant residue, the number reporting the mutant site, and the final letter reporting the mutant amino acid.
	* `ddg`: the experimentally measured ddG value for a given mutation, computed from Kd and temperature as described in the manuscript.
	* `beta16_cart`: the computationally estimated ddG value for a given mutation, estimated using the beta_nov16 energy function.
	* `beta16_cart_plus_hpsc_lj_changes`: the computationally estimated ddG value for a given mutation, estimated using the beta_jan25 energy function.
* `computed_ddg_values.csv`: similar to the above dataframe, but reporting additional data for computed ddG values. Each row reports data for a given mutation to a given structure with effects measured using a given energy function, with columns reporting:
	* `wt_complex` and `wt_apart`: dG values of the wildtype protein with the two sides of the interface either in complex or pulled apart, with values computed using the Cartesian ddG protocol, averaged across three iterations.
	* `mut_complex` and `mut_apart`: same as above, but for the mutant protein.
