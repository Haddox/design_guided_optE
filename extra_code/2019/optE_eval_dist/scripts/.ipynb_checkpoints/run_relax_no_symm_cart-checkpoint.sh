#!/bin/bash

pdb=$1
flags=$2
weights=$3
output_dir=$4

~/Rosetta/main/source/bin/relax.default.linuxgccrelease \
	@$flags \
	-s $pdb \
	-score:weights $weights \
	-crystal_refine \
	-nstruct 1 \
	-no_optH false -ignore_unrecognized_res -overwrite \
	-relax:script MonomerDesign2019 \
    -relax:cartesian true \
	-relax:min_type lbfgs_armijo_nonmonotone \
    -multi_cool_annealer 10 \
    -ex1 -ex2aro -linmem_ig 10 \
	-out:prefix $output_dir \
	#-mute all
