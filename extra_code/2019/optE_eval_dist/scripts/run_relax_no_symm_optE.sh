#!/bin/bash

pdb=$1
flags=$2
weights=$3
output_dir=$4

~/Rosetta_optE/main/source/bin/relax.default.linuxgccrelease \
	@$flags \
	-s $pdb \
	-score:weights $weights \
	-crystal_refine \
	-nstruct 1 \
	-no_optH false -ignore_unrecognized_res -overwrite \
	-default_max_cycles 200 \
	-relax:script cartminpack.script \
	-relax:min_type lbfgs_armijo_nonmonotone \
	-out:prefix $output_dir \
	#-mute all
