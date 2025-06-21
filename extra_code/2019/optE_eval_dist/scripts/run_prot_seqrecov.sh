#!/bin/sh

pdb=$1
flags=$2
weights=$3
output_dir=$4
pdb_bn=$(basename "$pdb" _clean.pdb)

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
        -database ~/Rosetta/main/database \
        -crystal_refine \
        -load_PDB_components false \
        -parser::protocol scripts/prot_fixbb.xml \
	-s $pdb \
	-in::file::native $pdb \
        @$flags \
        -score::weights $weights \
        -parser::script_vars weights_file=$weights \
	-out:prefix $output_dir &> $output_dir/$pdb_bn.prot
        #-out::nooutput \
        #-out::suffix _1 &> opt_$num/$pdb.prot
