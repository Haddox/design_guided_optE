"""
Functions that are useful for putting together an order of a chip of oligos of protein designs
"""

import os
import random
import doctest
import pandas

def pad_sequences_on_end(seq, pattern, length):
    """
    Pad an input sequence with a repeating pattern (e.g., `GSS`) until it reaches a target length.
    
    The pattern will be cut off at the end of the sequence if the target length is reached in the
    middle of the pattern.
    
    Args:
        `seq`: the input protein sequence (string)
        `pattern`: a sequence of amino acids that will be used to pad the sequence in a
            repeated fashion until it reaches the target amino-acid lentgh (string)
        `length`: the target amino-acid length (int)
    
    Returns:
        The padded sequence of the target length (string)
        
    Code for `doctest`:
    
    >>> seq = 'AQNP'
    >>> pattern = 'GSS'
    >>> length = 4
    >>> pad_sequences_on_end(seq, pattern, length)
    'AQNP'
    
    >>> seq = 'AQNP'
    >>> pattern = 'GSS'
    >>> length = 10
    >>> pad_sequences_on_end(seq, pattern, length)
    'AQNPGSSGSS'
    
    >>> seq = 'AQNP'
    >>> pattern = 'GSS'
    >>> length = 9
    >>> pad_sequences_on_end(seq, pattern, length)
    'AQNPGSSGS'
    
    >>> seq = 'AQNP'
    >>> pattern = 'GSS'
    >>> length = 8
    >>> pad_sequences_on_end(seq, pattern, length)
    'AQNPGSSG'
    """
    # Make sure the input sequence is less than or equal to the target length
    assert len(seq) <= length, "The length of the sequence is greater than the target length"
    
    # Pad the sequence using the input pattern until it reaches the desired length.
    i = 0
    while len(seq) < length:
        aa = pattern[i]
        seq += aa
        if i == len(pattern)-1:
            i = 0
        else:
            i += 1
    return seq

def make_patterned_scrambled_sequences(seq):
    """
    Scramble the input sequence while preserving both the overall amino-acid composition and the hydrophobic or polar
    character at each position
    
    Args:
        `seq`: the input amino-acid sequence that will be used to generate the patterned, scrambled output sequence
            (string)
    
    Returns:
        A scrambled output sequence (string)
    """
    
    # Set the random seed
    random.seed(a=3)
    
    # Define which amino acid is in each group
    aa_pattern = {
        'R': 'P',
        'N': 'P',
        'D': 'P',
        'Q': 'P',
        'E': 'P',
        'H': 'P',
        'K': 'P',
        'S': 'P',
        'T': 'P',
        'C': 'P',
        'I': 'H',
        'L': 'H',
        'M': 'H',
        'F': 'H',
        'W': 'H',
        'Y': 'H',
        'V': 'H',
        'A': 'H',
        'P': 'X',
        'G': 'Z'
    }
    assert len(aa_pattern) == 20
    patterns_not_to_shuffle = ['X', 'Z']
    
    # Read in the input sequence, convert it to a dataframe with each amino-acid as a row
    df = pandas.DataFrame.from_dict({'amino_acid':list(seq)})
    
    # Add a column indicating the group of each amino acid and make lists of the indices
    # of either all polar amino acids or all hydrophobic amino acids
    df['hydrophilicity'] = df['amino_acid'].apply(lambda x: aa_pattern[x])
    starting_hydrophilicity_string = list(df['hydrophilicity'])
    polar_indices = list(df[df['hydrophilicity'] == 'P'].index.values)
    hydrophobic_indices = list(df[df['hydrophilicity'] == 'H'].index.values)
    indices_to_keep_in_place = list(df[df['hydrophilicity'].isin(patterns_not_to_shuffle)].index.values)
    
    # Scramble the indices for the different groups and add these scrambled indices to the
    # dataframe in a column called `new_index`
    scrambled_polar_indices = random.sample(polar_indices, len(polar_indices))
    scrambled_hydrophobic_indices = random.sample(hydrophobic_indices, len(hydrophobic_indices))
    assert len(scrambled_polar_indices) == len(polar_indices)
    assert len(scrambled_hydrophobic_indices) == len(hydrophobic_indices)
    old_indices = polar_indices + hydrophobic_indices + indices_to_keep_in_place
    new_indices = scrambled_polar_indices + scrambled_hydrophobic_indices + indices_to_keep_in_place
    assert len(old_indices) == len(new_indices)
    scramble_indices_dict = {
        old_i : new_i
        for (old_i, new_i) in zip(old_indices, new_indices)
    }
    df['scramble_index'] = df.apply(lambda row: scramble_indices_dict[row.name], axis=1)
    
    # Sort the dataframe by the scramble index and get the scrambled amino-acid sequence
    df.sort_values(by='scramble_index', inplace=True)
    ending_hydrophilicity_string = list(df['hydrophilicity'])
    assert starting_hydrophilicity_string == ending_hydrophilicity_string, "The hydrophilicity has changed"
    scrambled_seq = ''.join(list(df['amino_acid']))

    return scrambled_seq

def WriteSbatchFile(sbatch_file_name, command_file_name=None, command=None, queue_type='short', memory='2g', array=None):
    """
    Write a file to submit a job via `sbatch`

    Args:
        `sbatch_file_name` : The name of the sbatch file. This will also
            serve as a prefix for the name of the output and error files.

        `command_file_name` : The name of a file with a list of commands
            to execute (string). The default is `None`, but if one is given,
            then this function will write an `sbatch` file that is specific
            for carying out an array of commands. If this is not given,
            the `command` argument must be specified. Only one can be
            specified at once.

        `command` : A command-line argument to be carried out (string).
            The default is `None`. If this is not given, the
            `command_file_name` argument must be specified. Only one can
            be specified at once.

        `queue_type` : The queue type ('short', 'medium', 'long')

        `array`: The size of the array (default: None)

	`memory` : The amount of memory in megabytes, unless other unit is specified. Default is '2g', i.e., 2 GB.

    Retruns:
        A file for submitting a job via `sbatch`
    """

    # If a file with an array of commands is provided, write an `sbatch`
    # file that is suited for this task
    if command_file_name:
        # Make sure the `command` argument hasn't been called
        assert command==None

        # Write the file
        with open(sbatch_file_name, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -p {0}\n'.format(queue_type))
            if array:
                f.write('#SBATCH -a {0}\n'.format(array))
            f.write('#SBATCH --mem={0}\n'.format(memory))
            f.write('#SBATCH -o {0}.out\n'.format(sbatch_file_name))
            f.write('#SBATCH -e {0}.err\n'.format(sbatch_file_name))
            f.write('CMD=$(head -n $SLURM_ARRAY_TASK_ID {0} | tail -1)\n'.format(command_file_name))
            f.write('exec ${CMD}')

    # Write an `sbatch` file to carry out the specified `command`
    else:
        with open(sbatch_file_name, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -p {0}\n'.format(queue_type))
            if array:
                f.write('#SBATCH -a {0}\n'.format(array))
            f.write('#SBATCH --mem={0}\n'.format(memory))
            f.write('#SBATCH -o {0}.out\n'.format(sbatch_file_name))
            f.write('#SBATCH -e {0}.err\n'.format(sbatch_file_name))
            f.write(command)


if __name__ == '__main__':
    doctest.testmod()