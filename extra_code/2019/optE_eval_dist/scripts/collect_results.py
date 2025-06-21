#!/usr/bin/env python

import os, sys
import argparse
import numpy
import scipy.special

class SeqrecovParser:
    def __init__(self):
        self.aamap = {
            'ALA':0,'CYS':1,'ASP':2,'GLU':3,'PHE':4,
            'GLY':5,'HIS':6,'ILE':7,'LYS':8,'LEU':9,
            'MET':10,'ASN':11,'PRO':12,'GLN':13,'ARG':14,
            'SER':15,'THR':16,'VAL':17,'TRP':18,'TYR':19
        }

    def parse(self, files):
        Nwrong,Ntotal=0,0
        dist_nat = numpy.zeros(20)
        dist_pred = numpy.zeros(20)
        for fname in files:

            with open(fname) as infile:
                designResBlock = False
                mutationsBlock = False
                for line in infile:
                    if (line.startswith(
                        "protocols.protein_interface_design.filters.DesignableResiduesFilter: Designable residues:")):
                        designResBlock = True
                    elif (line.startswith(
                        "protocols.protein_interface_design.filters.DesignableResiduesFilter: Number of design positions:")):
                        designResBlock = False
                    elif (line.startswith(
                        "protocols.protein_interface_design.filters.SequenceRecoveryFilter: Your design mover mutated")):
                        mutationsBlock = True
                    elif (mutationsBlock and line.startswith(
                        "protocols.protein_interface_design.filters.SequenceRecoveryFilter:")):
                        mutationsBlock = False
                    elif ( designResBlock ):
                        #protocols.protein_interface_design.filters.DesignableResiduesFilter: ARG 134A
                        fields = line.split()
                        Ntotal += 1
                        aaid = self.aamap[fields[1]]
                        dist_nat[aaid] += 1
                        dist_pred[aaid] += 1
                    elif ( mutationsBlock ):
                        #protocols.protein_interface_design.design_utils: THR20 4.77015	ARG20 2.07184
                        fields = line.split()
                        Nwrong += 1
                        srcid = self.aamap[fields[3][:3]]
                        tgtid = self.aamap[fields[1][:3]]
                        dist_pred[srcid] -= 1
                        dist_pred[tgtid] += 1
                    if (line.startswith(
                        "protocols.rosetta_scripts.ParsedProtocol.REPORT: ============End report for timed==================")):
                        break

        dist_nat = (dist_nat+0.01)/(Ntotal+0.2)
        dist_pred = (dist_pred+0.01)/(Ntotal+0.2)

        KLdiv = numpy.sum( scipy.special.rel_entr( dist_nat,dist_pred ) )

        print (1.0-Nwrong/Ntotal,KLdiv,Ntotal)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='collect seqrecov results')
    parser.add_argument('outfiles', help='outfiles', nargs='+')
    args = parser.parse_args()

    print(args.outfiles)

    sp = SeqrecovParser()
    sp.parse(args.outfiles)
