This directory contains input files and code for performing the atom-pair distance-distribution benchmark used to train and evaluate the energy function.
* `pdbs_for_training`: 78 high-resolution crystal structures used for training
* `pdbs_for_validation`: 54 high-resolution crystal structures used for validation
* `example_commands_to_run_test.txt`: example commands for computing atom-pair distance distributions from a set of PDBs of high-resolution crystal structures before and after relaxing them. This file also gives an example command for relaxing the structures.
* `distdstr_0.3.py`: Python script for computing atom-pair distance distributions given a set of input PDBs
* `data/`: contains files that encode the atom-typing scheme used to compute distance distributions. The above Python script uses these files.
	* The new atom-typing scheme classifies carbons from sidechains in the following ways:
		* CHR1: carbons from polar sidechains that are directly bonded to a polar atom
		* CHR2: carbons from polar sidechains that are not directly bonded to a polar atom
		* CH1, CH2, CH3, and aroC: carbons from non-polar sidechains
