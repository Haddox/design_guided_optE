* `dms_data_library_1.csv`: this CSV file reports DMS data from the two replicate experiments performed on the library of saturation mutagenesis sequences for 19 of 21 designs. Each row reports data for a given sequence from the library. The column names are the same as the ones used in Rocklin et al. with the following modifications:
	* Column titles with suffixes of `_rep1` and `_rep2` report a given quantity for either the first or second experimental replicate. Removing these suffixes should give column names that are the same as the ones reported in Rocklin et al.
	* `stabilityscore_avg` gives the average stability score between the two replicates. For a given replicate, the stability score is taken to be the minimum between proteases, as in Rocklin et al.
	* `parent_name` is the name of the design being mutagenized from Rocklin et al.
	* `site_n` is the site at which a given variant is mutated
	* `wt_aa` is the original designed amino acid at a given site
	* `mut_aa` is the mutant amino-acid identity in a given variant
	* `protein_sequence` is the full sequence of the protein variant
	* `chip_name` is the name of the oligonucleotide
* `dms_data_library_2.csv`: this CSV file is the same as the one above, except it reports DMS data for the single experiment performed on the library of saturation mutagenesis sequences for the remaining 2 of 21 designs. As there are no replicates, there are no replicate-specific suffixes or averaging of stability scores between replicates.
