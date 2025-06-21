#!/usr/bin/python

import scipy.signal
import scipy.stats
import numpy
import pandas
import os
import glob
import argparse
import warnings
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Define bins for binning inter-atomic distances
BINS = numpy.arange(1.5,5.05,0.05)
P_PSEUDO = 0.01

# Define the minimum number of counts per pair of atoms
# for that pair to be included in the analysis
MINCOUNT = 100


atom_dtype = numpy.dtype( [
    ("atomname", numpy.unicode_, 4),
    ("atomtype", numpy.int32),
    ("resname", numpy.unicode_, 3),
    ("chnid", numpy.unicode_, 1),
    ("resid", numpy.int32),
    ("X", numpy.float64, 3),
] )

def parse_pdb(pdbfile, aamap, allatms):
    """
    Extract residue/atom names and XYZ coordinates of each atom in
    an input PDB file
    
    Args:
        *pdbfile*: the path to an input PDB file
        *aamap*: a dictionary with the same format as the one returned
            by the *read_params* function
        *allatms*: a list with all Rosetta atom types to be considered
    
    Returns:
        A numpy array with one row for each atom in the
            input PDB and with the following columns:
                *atomname*: the atom name from the PDB file (e.g., N)
                *atmtype*: the index of the corresponding Rosetta atom
                    type in the list *allatms*
                *resname*: the residue name from the PDB file
                *chnid*: the chain ID from the PDB file
                *resid*: the residue ID from the PDB file
                *X*: a tuple with the XYZ coordiates of the atom
    """
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
    """
    Extract info from a Rosetta `.params` file
    
    Args:
        *paramfile*: the path the an input params file
        
    Returns:
        *aamap*: A dictionary keyed first by the name of an
            amino acid (e.g., ARG) and next by the standard
            name of an atom (e.g., N), returning the
            corresponding atom type in Rosetta (e.g., Nbb),
            omitting atom types that are in the list called
            `SKIPATOMS`
        *all_aas*: A set of all amino-acid names seen in
            the input params file
        *all_atms*: A set of all Rosetta atom types see in
            the input params file
    """
    
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
                if (atype not in SKIPATOMS):
                    all_atms.add(atype)
                    aamap[aaname][atmname] = atype

    return aamap,all_aas,all_atms


def load_db():
    """
    Read in metadata that allows mapping between atoms in a
    residue and atom types in Rosetta
    
    Returns:
        The same things as the `read_params` function, but with
            metadata for all residues/atoms from all params files,
            and returning lists instead of the sets
    """
    
    # Iterate over the pre-defined list of folders with params files
    # and record metadat from each file
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
    """
    Read in atomic radii for each atom type in Rosetta
    
    Args:
        *all_atms*: a list of all Rosetta atom types to
            consider
    Returns:
        *radii*: a numpy array with atomic radii listed
            in the same order as their corresponding in
            *all_atms*
    """
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

    # 
    MAXVAL=numpy.max(BINS)
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
                        numpy.abs(pdbdata[idxs1]['resid'][:,None] - pdbdata[idxs2]['resid'][None,:]) > SEQSEP #9
                    )
                )
                if (offset == 0):
                    mask = numpy.tril(mask, -1)
                mask = mask.nonzero()
                nmask = len(mask[0])

                dists[outarray_i:(outarray_i+nmask)] = d_ij[mask]
                types1[outarray_i:(outarray_i+nmask)] = \
                    pdbdata[idxs1[mask[0]]]['atomtype']
                types2[outarray_i:(outarray_i+nmask)] = \
                    pdbdata[idxs2[mask[1]]]['atomtype']

                outarray_i += nmask
    return dists[:outarray_i],types1[:outarray_i],types2[:outarray_i]


def main(
        files, reffile, radius, genref, plot, verbose, chainids,
        output_dir, allow_missing_atom_pairs, analyze_all_atoms,
        use_new_atom_typing
    ):
    """
    The main code of the script
    """
    
    # Read in metadata that allows mapping between each atom in each
    # residue and the corresponding atom type in Rosetta
    aamap,all_aas,all_atms = load_db()
    
    # Read in radii for atom type
    radii = load_radii(all_atms)

    # Iterate over all input PDB files and compute distances
    # between all pairs of atoms
    natms = len(all_atms)
    dlist = {}
    for pdbfile in files:
        
        # Read in a matrix with one row for each atom in the input
        # PDB and columns giving atom/residue name and coordinates
        pdbdata = parse_pdb(pdbfile,aamap,all_atms)
        
        # If specified, subset to atoms from indicated chains
        if (len(chainids) > 0):
            pdbdata = numpy.array(
                [ai for ai in pdbdata if ai['chnid'] in chainids], 
                dtype=atom_dtype)
            
        # Compute all inter-atomic distances between nearby atoms 
        ds,a1,a2 = get_neighbors(pdbdata)

        # Get data for atom pairs within a distance cutoff defined
        # by the sum of the atomic radii plus a constant
        vdwsum = radii[a1]+radii[a2]+radius
        ds,a1,a2 = ds[ds<vdwsum],a1[ds<vdwsum],a2[ds<vdwsum]
        
        for x in range(len(ds)):
            ai,aj = int(a1[x]),int(a2[x])
            if (ai<aj):
                tag = str(ai)+' '+str(aj)
            else:
                tag = str(aj)+' '+str(ai)
            if (tag not in dlist):
                dlist[tag] = []
            dlist[tag].append(ds[x])

    # Compute counts for each atom pair, and smooth the
    # distribution using a Gaussian KDE function
    hist = {}
    counts = {}
    for x,y in dlist.items():
        counts[x] = len(y)
        if (counts[x] > 1) and (len(set(y)) > 1):
            kde = scipy.stats.gaussian_kde(y, 'scott')
            hist[x] = kde(BINS)

    if (genref):
        toplot = {}
        with open(reffile,'w') as outfile:
            for x,y in hist.items():
                cij = int(counts[x])
                i,j = [int(k) for k in x.split()]
                outfile.write(
                    "{} {} {}".format(all_atms[i],all_atms[j],cij)
                )
                for y_i in y:
                    outfile.write(" {:.6f}".format(y_i))
                outfile.write("\n")

                if (plot):
                    name = all_atms[i]+':'+all_atms[j]+'\n'+str(cij)
                    toplot[name] = y

        if (plot):
            nplots = len(toplot.keys())
            fig, axs = plt.subplots(
                int((nplots+9)/10), 10, figsize=[40,40], frameon=False
            )
            for count,(name,ref) in enumerate(toplot.items()):
                i,j = int(count/10),count%10
                axs[i,j].set_title(name)
                axs[i,j].plot(ref,'r')
                axs[i,j].set_yticklabels([])
                axs[i,j].set_xticklabels([])
            plt.savefig(os.path.join(output_dir, 'ref.pdf'))  
            
    else:
        # compare to ref
        refhist, refcounts = dict(), dict()
        with open(reffile) as infile:
            for line in infile:
                fields = line.split()
                if (int(fields[2]) >= MINCOUNT):
                    refhist[ (fields[0],fields[1]) ] = numpy.array([float(x) for x in fields[3:]])
                    refcounts[ (fields[0],fields[1]) ] = int(fields[2])

        KLsum, KLwtsum = 0.0, 0.0
        toplot = {}
        verbose_output = []
        dist_dfs = []
        for (a1,a2) in refhist.keys():
            i1 = all_atms.index(a1)
            i2 = all_atms.index(a2)
            if (i1<i2):
                tag = str(i1)+' '+str(i2)
            else:
                tag = str(i2)+' '+str(i1)
            
            # If specified, skip over atom pairs that are present in
            # the reference proteins, but not in the input proteins.
            # Otherwise, raise an error
            if tag not in list(hist.keys()):
                if allow_missing_atom_pairs:
                    continue
                else:
                    raise ValueError(f'Atom pair {tag} not in input proteins')

            # Compute the KL divergence of the atom pair, weighing it
            # by the number of counts of the atom pair
            kldiv = (BINS[1]-BINS[0]) * numpy.sum(scipy.special.rel_entr(hist[tag]+P_PSEUDO, refhist[a1,a2]+P_PSEUDO))
            weight = numpy.sqrt( refcounts[a1,a2] )
            verbose_output.append('{0} {1} {2} {3}'.format(a1, a2, kldiv, weight))
            KLwtsum += weight
            KLsum += weight * kldiv * kldiv

            if (plot):

                # Add data to dictionary for plotting
                name = a1+':'+a2+'('+str(refcounts[a1,a2])+')\n'+str(kldiv)
                toplot[name] = (refhist[a1,a2], hist[tag])
                
                # Add data to dataframe for output CSV
                dist_df = pandas.DataFrame({
                    'distance' : BINS,
                    'ref' : refhist[a1,a2],
                    'input' : hist[tag]
                })
                dist_df['a1'] = a1
                dist_df['a2'] = a2
                dist_df['atoms'] = a1+':'+a2
                dist_dfs.append(dist_df)
        
        dist_df = pandas.concat(dist_dfs, sort=False)
        output_dist_file = os.path.join(output_dir, 'distrs.csv')
        if use_new_atom_typing:
            output_dist_file = output_dist_file.replace(
                '.csv', '_new_atom_typing.csv'
            )
        if analyze_all_atoms:
            output_dist_file = output_dist_file.replace(
                '.csv', '_all_atoms.csv'
            )
        dist_df.to_csv(output_dist_file, index=False)

        if (plot):
            nplots = len(toplot.keys())
            fig, axs = plt.subplots(int((nplots+9)/10), 10, figsize=[40,40], frameon=False)
            for count,(name,(d1,d2)) in enumerate(toplot.items()):
                i,j = int(count/10),count%10
                axs[i,j].set_title(name)
                axs[i,j].plot(d1,'r')
                axs[i,j].plot(d2,'g')
                axs[i,j].set_yticklabels([])
                axs[i,j].set_xticklabels([])
            output_plot_file = os.path.join(output_dir, 'distrs.pdf')
            if use_new_atom_typing:
                output_plot_file = output_plot_file.replace(
                    '.pdf', '_new_atom_typing.pdf'
                )
            if analyze_all_atoms:
                output_plot_file = output_plot_file.replace(
                    '.pdf', '_all_atoms.pdf'
                )
            plt.savefig(output_plot_file)

        ndstrs = len(refhist.keys()) 
        print('%f %d'%(sqrt(KLsum/KLwtsum), ndstrs)) # HH RMS divergence?
        if verbose:
            print('a1 a2 kldiv weight')
            for output in verbose_output:
                print(output)


if __name__ == "__main__":
    
    # Read in command-line arguments
    parser = argparse.ArgumentParser(description='Compute interatomic distance distributions')
    parser.add_argument('inpdb', help='outfiles', nargs='+')
    parser.add_argument('--output_dir',  help='output directory in which to store the plot, if it is generated', default='./')
    parser.add_argument('--gen',  help='Generate reference file', action='store_true')
    parser.add_argument('--ref',  help='Reference file name', default='distr_new.REF')
    parser.add_argument('--plot',  help='Plot distributions', action='store_true')
    parser.add_argument('--verbose',  help='Be verbose', action='store_true')
    parser.add_argument('--chains',  help='Chain IDs', default='')
    parser.add_argument('--radius',  help='Radius of interactions (w.r.t. vdw radii)', default=0.5)
    parser.add_argument('--min_seq_sep',  help='Minimum separation between residues in primary sequence for data between them to be counted', type=float, default=10)
    parser.add_argument('--allow_missing_atom_pairs',  help='Allow atom pairs that are present in the reference file to be absent in the set of input proteins', action='store_true')
    parser.add_argument('--inpdb_arg_is_dir',  help='The argument giving input PDBs is a directory', action='store_true')
    parser.add_argument('--use_all_files_in_dir',  help='If the inpdb_arg_is_dir argument is provided, then use all files in that directory, not just ones with .pdb extensions ', action='store_true')
    parser.add_argument('--analyze_all_atoms',  help='Analyze all atom types in Rosetta. Otherwise, the default is to exclude atoms in the list named SKIPATOMS', action='store_true')
    parser.add_argument('--use_new_atom_typing',  help='Use the new atom-typing scheme with CHR1/CHR2/HapR/HAbb', action='store_true')
    args = parser.parse_args()

    # If an inputdirectory was provided instead of PDB paths,
    # then make a list of PDB paths
    if args.inpdb_arg_is_dir:
        if args.use_all_files_in_dir:
            pdbs = glob.glob(os.path.join(args.inpdb[0], '*'))
        else:
            pdbs = glob.glob(os.path.join(args.inpdb[0], '*.pdb'))
        assert len(pdbs) > 0, len(pdbs)
    else:
        pdbs = args.inpdb
        
    # Define the variable related to minimum sequence separation
    # between residues for data to be counted
    global SEQSEP
    SEQSEP = args.min_seq_sep - 1
        
    # Determine which chains to analyze
    chains = args.chains
    if (chains != ''):
        chains = chains.split(',')
    else:
        chains = []
        
    # Define paths to input files with atom types and radii
    global ATOMFILE
    global DBPATHS
    ROSETTADB = '/home/dimaio/Rosetta/main/database'
    if args.use_new_atom_typing:
        ATOMFILE = \
            'data/new_atom_typing/atom_properties.txt'
        DBPATHS = [
            'data/new_atom_typing/params_files/',
            '%s/chemical/residue_type_sets/fa_standard/residue_types/nucleic/dna'%ROSETTADB
        ]
    else:
        ATOMFILE = '%s/chemical/atom_type_sets/fa_standard/atom_properties.txt'%ROSETTADB
        DBPATHS = [
            '%s/chemical/residue_type_sets/fa_standard/residue_types/l-caa'%ROSETTADB,
            '%s/chemical/residue_type_sets/fa_standard/residue_types/nucleic/dna'%ROSETTADB
        ]
        
    # Define atoms to skip over in the analysis
    if args.analyze_all_atoms:
        SKIPATOMS = []
    else:
        SKIPATOMS = ['Hapo','Hpol','Haro','HS', 'HNbb', 'HAbb']

    # Run the main code
    main(
        pdbs,
        args.ref,
        args.radius,
        args.gen,
        args.plot,
        args.verbose,
        chains,
        args.output_dir,
        args.allow_missing_atom_pairs,
        args.analyze_all_atoms,
        args.use_new_atom_typing
    )
