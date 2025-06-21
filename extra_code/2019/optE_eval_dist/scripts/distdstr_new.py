#!/usr/bin/python

import scipy.signal
import numpy
import os
import argparse
import warnings
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROSETTADB = '/home/dimaio/Rosetta/main/database'
ATOMFILE = '%s/chemical/atom_type_sets/fa_standard/atom_properties.txt'%ROSETTADB
DBPATHS = [
    '%s/chemical/residue_type_sets/fa_standard/residue_types/l-caa'%ROSETTADB,
    '%s/chemical/residue_type_sets/fa_standard/residue_types/nucleic/dna'%ROSETTADB
]

MINVAL = 1.5
MAXVAL = 5.0
BINSIZE = 0.05
NBINS = int((MAXVAL-MINVAL)/BINSIZE)
P_PSEUDO = 0.01

def get_bin(d_i):
    bin_i = int((d_i-MINVAL)/BINSIZE)
    if (bin_i<0):
        bin_i = 0
    if (bin_i > NBINS):
        bin_i = NBINS
    return bin_i

def get_bins(ds_i):
    bins_i = ((ds_i-MINVAL)/BINSIZE).astype(numpy.int32)
    bins_i[bins_i<0]=0
    bins_i[bins_i>NBINS]=NBINS
    return bins_i

atom_dtype = numpy.dtype( [
    ("atomname", numpy.unicode_, 4),
    ("atomtype", numpy.int32),
    ("resname", numpy.unicode_, 3),
    ("chnid", numpy.unicode_, 1),
    ("resid", numpy.int32),
    ("X", numpy.float64, 3),
] )

def parse_pdb(pdbfile, aamap, allatms):
    allatoms = []
    with open(pdbfile) as pdbin:
        for line in pdbin:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if ((resname in aamap) and (atomname in aamap[resname])):
                    atmtype = allatms.index(aamap[resname][atomname])
                    split_line = (
                        #line[12:16], line[17:20], line[21], line[22:26], 
                        atomname, atmtype, resname, line[21], line[22:26], 
                        (line[30:38], line[38:46], line[46:54])
                    )
                    allatoms.append(split_line)
    return (numpy.array(allatoms, dtype=atom_dtype))

def read_params(paramfile):
    all_aas = set()
    all_atms = set()
    aamap = dict()
    with open(paramfile) as paramsin:
        aaname = ''
        for line in paramsin:
            if line.startswith('IO_STRING '):
                fields = line.rstrip().split()
                aaname = fields[1]
                all_aas.add(aaname)
                aamap[aaname] = dict()
            elif line.startswith('ATOM '):
                fields = line.rstrip().split()
                atmname = fields[1]
                atype = fields[2]
                all_atms.add(atype)
                aamap[aaname][atmname] = atype

    return aamap,all_aas,all_atms


def load_db():
    all_aas = set()
    all_atms = set()
    aamap = dict()

    for dbpath in DBPATHS:
        for root, dirs, files in os.walk(dbpath):
            for file in files:
                if (file[0] == '.'):
                    continue
                resfile = os.path.join(root, file)
                aamap_i, aas_i, atoms_i = read_params(resfile)
                all_aas.update(aas_i)
                all_atms.update(atoms_i)
                aamap.update(aamap_i)
    return aamap,list(all_aas),list(all_atms)


def load_radii(all_atms):
    radii = numpy.zeros(len(all_atms))
    with open(ATOMFILE) as instream:
        for line in instream:
            words = line[:-1].split()
            atype = words[0]
            if atype in all_atms:
                radius = float(words[2])
                ai = all_atms.index(atype)
                radii[ai] = radius
    return radii

    
def get_neighbors(pdbdata):
    minX = numpy.min(pdbdata['X'],axis=0)
    hashbins = (
        numpy.floor((pdbdata['X'] - minX[None,:])/(MAXVAL+1.0))
    ).astype(numpy.int64)
    FACT = 100
    indices = FACT*FACT*hashbins[:,0]+FACT*hashbins[:,1]+hashbins[:,2]

    nbins = 0
    maxbin = 0
    pthash = dict()
    for i,ind in enumerate(indices):
        if (not ind in pthash):
            nbins += 1
            pthash[ind] = numpy.array([], dtype=numpy.int32)
        pthash[ind] = numpy.concatenate( (pthash[ind], numpy.array([i])) )
        if (len(pthash[ind]) > maxbin):
            maxbin = len(pthash[ind])
        
    halfplane = [
        0, 1,
        FACT-1, FACT, FACT+1,
        FACT*FACT-FACT-1, FACT*FACT-FACT, FACT*FACT-FACT+1,
        FACT*FACT-1, FACT*FACT, FACT*FACT+1,
        FACT*FACT+FACT-1, FACT*FACT+FACT, FACT*FACT+FACT+1,
    ]

    MAXLEN = 14 * nbins * maxbin * maxbin
    dists = numpy.zeros(MAXLEN)
    types1 = numpy.zeros(MAXLEN, dtype=numpy.int32)
    types2 = numpy.zeros(MAXLEN, dtype=numpy.int32)

    outarray_i = 0
    for bin1 in pthash.keys():
        idxs1 = pthash[bin1]

        for offset in halfplane:
            bin2 = bin1+offset
            if (bin2 in pthash):
                idxs2 = pthash[bin2]
                d_ij = numpy.linalg.norm(
                    pdbdata[idxs1]['X'][:,None,:] - pdbdata[idxs2]['X'][None,:,:], axis=-1 )
                mask = numpy.logical_and(
                    d_ij<MAXVAL,
                    numpy.logical_or(
                        pdbdata[idxs1]['chnid'][:,None] != pdbdata[idxs2]['chnid'][None,:],
                        numpy.abs(pdbdata[idxs1]['resid'][:,None] - pdbdata[idxs2]['resid'][None,:]) > 9
                    )
                )
                if (offset == 0):
                    mask = numpy.tril(mask, -1)
                mask = mask.nonzero()
                nmask = len(mask[0])

                dists[outarray_i:(outarray_i+nmask)] = d_ij[mask]
                types1[outarray_i:(outarray_i+nmask)] = pdbdata[idxs1[mask[0]]]['atomtype']
                types2[outarray_i:(outarray_i+nmask)] = pdbdata[idxs2[mask[1]]]['atomtype']

                outarray_i += nmask
    return dists[:outarray_i],types1[:outarray_i],types2[:outarray_i]


def main(files, reffile, radius, genref, plot, verbose, chainids, output_dir, min_counts):
    MINCOUNT = min_counts
    aamap,all_aas,all_atms = load_db()
    radii = load_radii(all_atms)
    natms = len(all_atms)
    hist = numpy.zeros((natms,natms,NBINS+1))

    for pdbfile in files:
        pdbdata = parse_pdb(pdbfile,aamap,all_atms)
        if (len(chainids) > 0):
            pdbdata = numpy.array(
                [ai for ai in pdbdata if ai['chnid'] in chainids], 
                dtype=atom_dtype)
        ds,a1,a2 = get_neighbors(pdbdata)
        vdwsum = radii[a1]+radii[a2]+radius
        ds,a1,a2 = ds[ds<vdwsum],a1[ds<vdwsum],a2[ds<vdwsum]

        # histograms
        indices = get_bins(ds)
        numpy.add.at(hist,(a1,a2,indices),1)
        numpy.add.at(hist,(a2,a1,indices),1)

    counts = numpy.sum(hist[:,:,:-1], axis=2)
    counts[counts==0] = P_PSEUDO
    conv = numpy.array([0.04,0.11,0.21,0.28,0.21,0.11,0.04])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        hcf=scipy.signal.convolve(hist[:,:,:-1], conv[None,None,:], mode='same')
    hcf = hcf / counts[:,:,None]

    if (genref):
        with open(reffile,'w') as outfile:
            for i in range(len(all_atms)):
                for j in range(i,len(all_atms)):
                    if (counts[i,j] > MINCOUNT):
                        outfile.write("{} {} {}".format(all_atms[i],all_atms[j],int(counts[i,j])))
                        for x_i in hcf[i,j,:]:
                            outfile.write(" {:.6f}".format(x_i))
                        outfile.write("\n")
    else:
        # compare to ref
        refhist = dict()
        with open(reffile) as infile:
            for line in infile:
                fields = line.split()
                refhist[ (fields[0],fields[1]) ] = numpy.array([float(x) for x in fields[3:]])

        KLsum = 0.0
        toplot = {}
        for (a1,a2) in refhist.keys():
            i1 = all_atms.index(a1)
            i2 = all_atms.index(a2)
            if (counts[i1,i2] > MINCOUNT):
                kldiv = numpy.sum(scipy.special.rel_entr(hcf[i1,i2]+P_PSEUDO, refhist[a1,a2]+P_PSEUDO))
                if (verbose):
                    print (a1,a2,kldiv)
                KLsum += kldiv * kldiv
                if (plot):
                    name = a1+':'+a2+'\n'+str(kldiv)
                    toplot[name] = (refhist[a1,a2], hcf[i1,i2])

        if (plot):
            nplots = len(toplot.keys())
            fig, axs = plt.subplots(
                int((nplots+9)/10), 10,
                figsize=[40,nplots/2.5],
                frameon=False
            )
            for count,(name,(d1,d2)) in enumerate(toplot.items()):
                i,j = int(count/10),count%10
                axs[i,j].set_title(name)
                axs[i,j].plot(d1,'r')
                axs[i,j].plot(d2,'g')
                axs[i,j].set_yticklabels([])
                axs[i,j].set_xticklabels([])
            plt.savefig(os.path.join(output_dir, f'distrs_{MINCOUNT}.pdf'))
        print('%f %d'%(sqrt(KLsum/count),count))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='collect ssm results')
    parser.add_argument('inpdb', help='outfiles', nargs='+')
    parser.add_argument('--output_dir',  help='output directory in which to store the plot showing the distribution of inter-atomic distances')
    parser.add_argument('--min_counts',  help='the minimum number of observations of an atom pair for it to be plotted and counted in the overall score', type=int)
    parser.add_argument('--gen',  help='Generate reference file', action='store_true')
    parser.add_argument('--ref',  help='Reference file name', default='distr_new.REF')
    parser.add_argument('--plot',  help='Plot distributions', action='store_true')
    parser.add_argument('--verbose',  help='Be verbose', action='store_true')
    parser.add_argument('--chains',  help='Chain IDs', default='')
    parser.add_argument('--radius',  help='Radius of interactions (w.r.t. vdw radii)', default=0.5)
    args = parser.parse_args()

    assert (not (args.gen and args.plot)), "--gen and --plot cannot both be set"

    chains = args.chains
    if (chains != ''):
        chains = chains.split(',')
    else:
        chains = []

    main(
        args.inpdb,
        args.ref,
        args.radius,
        args.gen,
        args.plot,
        args.verbose,
        chains,
        args.output_dir,
        args.min_counts
    )
