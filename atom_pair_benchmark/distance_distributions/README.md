* `beta_nov16_distributions.csv` and `beta_jan25_distributions.csv`: files with atom-pair distance distributions computed for the 54 high-resolution crystal structures used to validate energy-function performance, where each file reports results for performing the test with a given energy function. Each row reports data for a given atom pair. Columns report:
	* `a1` and `a2`: the names of the two atoms in a given pair
	* `distance`: a distance in Angstroms (distances are binned every 0.05 Angstroms)
	* `ref`: the density of the number of instances of a given atom pair where the atoms are separated by a given distance fall. This column reports the density observed in crystal structures before they are relaxed.
	* `input`: the same as `ref`, but reporting the density observed in structures after they are relaxed.
