Each design/crystal pair has one PDB with the suffix `_design.pdb` and another PDB with the suffix `_xtal.pdb`.
The first corresponds to the original design model.
The second corresponds to the corresponding crystal structure.
I downloaded the crystal structure from the PDB and removed any molecules not part of the protein chains in the original design model (e.g., water molecules).
Some crystal structures had multiple copies of the design, in which case I downselected to a single copy (see below for details).
For one crystal structure, I generated a symmetry partner to get the complete design (see below).
* DHD_131: chains AB
* D_3_212
	* fetched PDB (7rmx)
	* created a symmetry partners in pymol using the command `symexp sym, 7rmx, 7rmx, 2` and found the partner that corresponds to the homodimer
	* renamed the chain of the symmetry partner to chain B using the command `alter sym, chain='B'`
	* created a single object with the command `create D_3_212, (7rmx or sym)`
	* saved the PDB with the command  
* DHR76: chain A
* WSHC6: chains DHIJKL
* 5L6HC3_1: chains ABC
* THR2: chain A
* LHD29: chains AB
* MC2_7: chain B
* BB1: chain A
