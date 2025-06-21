#!/bin/bash

trg=$1
temp=$2
istart=$3

cd $trg

/software/rosetta/latest/bin/rosetta_scripts \
    @../flags.liquidsim -s restart.pdb \
    -parser:script_vars nmol=4 scfilename=prod.sc outname=prod.out density=1.0 temperature=$temp istart=$istart \
    > prod.log
