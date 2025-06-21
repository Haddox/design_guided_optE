#!/bin/bash

pdb=$1
flags=$2
weights=$3
output_dir=$4

~/Rosetta/main/source/bin/relax.default.linuxgccrelease \
        -database ~/Rosetta/main/database \
	-s $pdb \
        -score:weights $weights \
        @$flags \
        -set_weights cart_bonded 0.5 pro_close 0.0 \
        -dna_move true \
        -exclude_dna_dna false \
        -crystal_refine \
        -nstruct 1 \
        -no_optH false -ignore_unrecognized_res -overwrite \
        -relax:script cartminpack.script \
        -out::prefix $output_dir
