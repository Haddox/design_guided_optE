from __future__ import division

# This program accepts arguments like this:

#./remove_superfluous_trp.py pdb1.pdb pdb2.pdb pdb3.pdb
# or
#./remove_superfluous_trp.py -in:file:silent my.silent

import os
import sys
import math

from pyrosetta import *
from pyrosetta.rosetta import *

import numpy as np
from collections import defaultdict
import time
import argparse
import itertools
import subprocess
import time
import scipy.spatial.transform
import tqdm

# Append the directory of the current script to sys.path
scripts_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(scripts_dir, 'npose'))
import npose_util as nu
import npose_util_pyrosetta as nup


init("-corrections:beta_nov16  -in:file:silent_struct_type binary -keep_input_scores false")


parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("file_listing_pdbs", type=str)
parser.add_argument("-max_translate", type=float, default=4, help="Max translation amount")
parser.add_argument("-max_angle", type=float, default=15, help="Max rotation angle")
parser.add_argument("-cluster_rmsd", type=float, default=1.2, help="RMSD to cluster before outputting")
parser.add_argument("-internal_samples", type=int, default=100000, help="Internally how many xforms to generate")
parser.add_argument("-silent_out", action="store_true", help="force output as silent file")
parser.add_argument("-output_dir", help="output directory where all results from script will be stored")
parser.add_argument("-output_pdb_prefix", help="prefix that will be added to the basename of output PDB files")

args = parser.parse_args(sys.argv[1:])

file_listing_pdbs = args.file_listing_pdbs
with open(file_listing_pdbs) as f:
    lines = f.readlines()
    pdbs = [line.strip() for line in lines]

print(pdbs)

silent = args.__getattribute__("in:file:silent")


scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16")


def my_rstrip(string, strip):
    if (string.endswith(strip)):
        return string[:-len(strip)]
    return string



the_locals = None

def jitter_binder(pose, name_no_suffix):
    """
    Jitter a binder in 3D space relative to a target
    """

    # Make objects that correspond to the chains of either the binder or the target
    chains = pose.split_by_chain()
    split_pose = pose.split_by_chain()
    if len(split_pose) != 2:
        return
    target = split_pose[1]
    binder = split_pose[2]

    N = args.internal_samples

    # we need unit vectors for translation and rotation
    # using a gaussian distribution for x, y, and z results in uniform random vectors
    # make 2*N so we can split it later
    random_units = np.random.normal(size=(N*2, 3))

    # they weren't normalized, so calculate the norm and divide
    random_units /= np.linalg.norm( random_units, axis=-1 )[:,None]

    # split the list in 2
    rotation_units = random_units[:N]
    translation_units = random_units[N:]

    # randomly select translation amount
    translation_scale = np.random.random(N) * args.max_translate

    # randomly select rotation amount
    rotation_scale = np.random.random(N) * args.max_angle


    # here comes the nasty part. transformations are given as 4 x 4 matrices that look like this
    #
    #   R  R  R  X
    #   R  R  R  Y
    #   R  R  R  Z
    #   0  0  0  1
    #
    # where R represents a 3x3 rotation matrix and X Y Z is the translation vector

    # first make the rotation matrices
    rotation_matrices = scipy.spatial.transform.Rotation.from_rotvec(rotation_units * np.radians( rotation_scale )[:,None] ).as_dcm()

    # next make the translation vectors
    translation_vectors = translation_units * translation_scale[:,None]

    # Now we make the transformations
    xforms = np.zeros((N, 4, 4), np.float)

    xforms[:,:3,:3] = rotation_matrices
    xforms[:,:3,3] = translation_vectors
    xforms[:,3,3] = 1


    # When you rotate something, its center of mass needs to be at the origin.
    # Let's get the xforms for that set up right now

    # get CA locations
    binder_cas = nup.get_pose_cas( binder )

    center_of_mass = nu.center_of_mass( binder_cas )

    move_com_to_orig = np.identity(4)
    move_orig_to_com = np.identity(4)

    move_com_to_orig[:3,3] = -center_of_mass
    move_orig_to_com[:3,3] = center_of_mass


    # here comes the sketchy part
    # we use Will's method to calculate RMSD from two different transforms

    # first calculate the radius of gyration of CAs for the binder

    radius_of_gyration = nu.radius_of_gyration( binder_cas )

    # Will's method to calculate RMSD from two different transforms
    #
    #   We want the RMSD between transform A and transform B
    #
    # 1. RelativeTransform = inv( A ) * B
    # 2. TranslationDist = RelativeTransform.translation
    # 3. AngularDist = RadiusOfGyration * sin( RelativeTransform.angle )
    # 4. RMSD = sqrt( TranslationDist^2 + AngularDist^2 )

    # Next we cluster by Will's RMSD using a function that brian wrote

    cluster_centers, _ = nu.cluster_xforms(args.cluster_rmsd, radius_of_gyration, xforms, info_every=100 )

    print("")
    print("Found %i clusters. Dumping poses"%len(cluster_centers))

    # We have the jitters! Let's start jittering

    for i, ixform in enumerate(tqdm.tqdm(cluster_centers)):

        # new pose to translate
        pose = binder.clone()

        # get all the information about this xform
        xform = xforms[ixform]
        rotation_angle = rotation_scale[ixform]
        translation = translation_scale[ixform]

        # the actual xform we're going to use needs to move the binder to the origin
        #   apply the xform, and move it back. Turns out this is just matrix multiplication

        actual_xform = move_orig_to_com @ xform @ move_com_to_orig

        # pull out the rotation matrix and the translation vector in rosetta data types
        rotation_matrix = nup.xform_to_matrix( actual_xform )
        translation_vector = nup.xform_to_vector( actual_xform )

        pose.apply_transform_Rx_plus_v( rotation_matrix, translation_vector )

        target_pose = target.clone()
        target_pose.append_pose_by_jump( pose, 1 )
        #pose.append_pose_by_jump( target, 1 )

        # prepare the score file
        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()

        score_map['translation'] = translation
        score_map['rotation_angle'] = rotation_angle

        # give it a new name
        name = os.path.join(
            args.output_dir,
            args.output_pdb_prefix + name_no_suffix + "_chain_{0}".format(binder.pdb_info().chain(1)) + "_%i"%i
        )

        # and then write it to score file
        sfd = core.io.raw_data.ScoreFileData(os.path.join(args.output_dir, "jitter_metadata.sc"))
        sfd.write_pose( target_pose, score_map, name, string_map)

        # and then either dump pdb or write to silent
        if ( silent != "" or args.silent_out ):
            silent_name = os.path.join(args.output_dir, "out.silent")
            sfd_out = core.io.silent.SilentFileData( silent_name, False, False, "binary", core.io.silent.SilentFileOptions())
            struct = sfd_out.create_SilentStructOP()
            struct.fill_struct(target_pose, name)
            sfd_out.add_structure(struct)
            sfd_out.write_all(silent_name, False)
        else:
            target_pose.dump_pdb(name + ".pdb")





############### BEGIN MAIN FUNCTION ###########################

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    pdbs = list(sfd_in.tags())


num = -1
for pdb in pdbs:
    t0 = time.time()
    print("Attempting pose: " + pdb)

    # try:
    for k in [1]:
        if ( silent == "" ):
            pose = pose_from_file(pdb)
        else:
            pose = Pose()
            sfd_in.get_structure(pdb).fill_pose(pose)

        name_no_suffix = my_rstrip(my_rstrip(os.path.basename(pdb), ".gz"), ".pdb")

        out_stuff = jitter_binder(pose, name_no_suffix )

        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)
