#!/bin/sh

pdb=$1
flags=$2
weights=$3
output_dir=$4

~/Rosetta_optE/main/source/bin/rosetta_scripts.default.linuxgccrelease \
        -database ~/Rosetta_optE/main/database \
        -crystal_refine \
        -load_PDB_components false \
        -parser::protocol scripts/prot_fixbb.xml \
	-s $pdb \
	-in::file::native $pdb \
        @$flags \
        -score::weights $weights \
        -parser::script_vars weights_file=$weights \
	-out:prefix $output_dir
        #-out::nooutput \
        #-out::suffix _1 &> opt_$num/$pdb.prot
