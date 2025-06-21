#!/bin/bash

trg=$1
temp=$2

cd $trg

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
    @../flags.liquidsim \
    -parser:protocol ../prod.xml \
    -parser:script_vars nmol=4 scfilename=prod.sc outname=prod.out density=1.0 temperature=$temp istart=1 \
    -s LG_0001.pdb > prod.log
