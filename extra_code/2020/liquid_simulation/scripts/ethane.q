#!/bin/bash
#SBATCH -p medium
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH -o log

cd ~/projects/genpot/optE/liquidsim/test
./run.sh ethane 185
