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

sys.path.append("/home/bcov/sc/scaffold_comparison/tools")
from elife_write_PDBcontacts import get_pi_pi_contacts

import distutils.spawn
sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
import silent_tools

init("-corrections:beta_nov16  -in:file:silent_struct_type binary " +
    "-holes:dalphaball /home/bcov/dev_rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc")


parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("-overwrite", action='store_true', default=False)


args = parser.parse_args(sys.argv[1:])

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")


# script_dir = os.path.dirname(os.path.realpath(__file__))
# xml = script_dir + "/py_xml/remove_superfluous_nonpolar.xml"


# objs = protocols.rosetta_scripts.XmlObjects.create_from_file(xml)

# Set up score functions
scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16")
blank_scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("none")
hb_scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("none")
hb_scorefxn.set_weight(core.scoring.hbond_sc, 1)
hb_scorefxn.set_weight(core.scoring.hbond_sr_bb, 1)
hb_scorefxn.set_weight(core.scoring.hbond_lr_bb, 1)
hb_scorefxn.set_weight(core.scoring.hbond_bb_sc, 1)

salt_scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("none")
salt_scorefxn.set_weight(core.scoring.hbond_sc, 1)

def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)

fix_scorefxn(hb_scorefxn, True)
fix_scorefxn(salt_scorefxn)



def my_rstrip(string, strip):
    if (string.endswith(strip)):
        return string[:-len(strip)]
    return string


def add_to_score_map(og_map, to_add, prefix, suffix=""):
    for name, score in list(to_add.items()):    # this iterator is broken. use list()
        og_map[prefix + name + suffix] = score

def move_chainA_far_away(pose):
    pose = pose.clone()
    sel = core.select.residue_selector.ChainSelector("A")
    subset = sel.apply(pose)

    x_unit = numeric.xyzVector_double_t(1, 0, 0)
    far_away = numeric.xyzVector_double_t(10000, 0, 0)

    protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, far_away, pose, subset)

    return pose


def which_chain(pose, resi):
    for i in range(1, pose.num_chains()+1):
        if ( pose.conformation().chain_begin(i) <= resi and pose.conformation().chain_end(i) >= resi):
            return i
    assert(False)

def get_monomer_score(pose, scorefxn):
    pose = pose.split_by_chain()[1]
    return scorefxn(pose)


def get_filter_by_name(filtername):
    the_filter = objs.get_filter(filtername)

    # Get rid of stochastic filter
    if ( isinstance(the_filter, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
        the_filter = the_filter.subfilter()

    return the_filter

def add_filter_to_results(pose, filtername, out_score_map):
    filter = get_filter_by_name(filtername)
    print("protocols.rosetta_scripts.ParsedProtocol.REPORT: ============Begin report for " + filtername + "=======================")
    if (isinstance(filter, protocols.simple_filters.ShapeComplementarityFilter)):
        value = filter.compute(pose)
        out_score_map[filtername] = value.sc
        out_score_map[filtername+"_median_dist"] = value.distance
    else:
        value = filter.report_sm(pose)
        out_score_map[filtername] = value
    print("============End report for " + filtername + "=======================")

def score_with_these_filters(pose, filters, out_score_map):
    for filtername in filters:
        add_filter_to_results(pose, filtername, out_score_map)

def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return str(the_stuff[0]) + str(the_stuff[1])

# types of burial:
# very-buried, sasa-buried

# types of satisfied
# rosetta-hbond, unsatisfied

# groups:
# arg, amide, bb-amide, carboxyl, his, ser, tyr

# atom types:
# arg, lys, bb-n, bb-o, amide-n, amide-o, carboxyl, tyr-o, tyr-h, ser-o, ser-h, trp-n, his-n


# modifiers
# amide-stack, carb-stack, ring-stack, arg-stack, lys-stack, salt-bridged

def atid(resnum, atno):
    return core.id.AtomID( atno, resnum )

the_locals = None

def types_of_buried_polars(pose, out_score_map, out_string_map, suffix):

    global the_locals

    # Setup the atom maps:
    #   The one to tell dalphaball where to look
    #   The one to identify if buried by sasa
    #   The one to identify if very buried

    atoms = core.id.AtomID_Map_bool_t()
    atoms.resize(pose.size())
    is_sasa_buried = core.id.AtomID_Map_bool_t()
    is_sasa_buried.resize(pose.size())
    is_very_buried = core.id.AtomID_Map_bool_t()
    is_very_buried.resize(pose.size())
    is_both_buried = core.id.AtomID_Map_bool_t()
    is_both_buried.resize(pose.size())
    bad_atom_sasa = core.id.AtomID_Map_double_t()
    bad_atom_sasa.resize(pose.size())
    for i in range(1, pose.size()+1):
        atoms.resize( i, pose.residue(i).natoms(), True)
        is_sasa_buried.resize( i, pose.residue(i).natoms(), False)
        is_very_buried.resize( i, pose.residue(i).natoms(), False)
        is_both_buried.resize( i, pose.residue(i).natoms(), False)
        is_both_buried.resize( i, pose.residue(i).natoms(), False)
        bad_atom_sasa.resize(i, pose.residue(i).natoms(), 0)


    # Do the burial calculations
    surf_vol = core.scoring.packing.get_surf_vol( pose, atoms, 0.9)

    bad_sasa = None
    for i in range(1, pose.size()+1):
        for j in range(1, atoms.n_atom(i)+1):
            if ( utility.isnan( surf_vol.surf( i, j ) ) ):
                print("BAD SASA")

                residue_sasa = utility.vector1_double()
                core.scoring.calc_per_atom_sasa( pose, bad_atom_sasa, residue_sasa, 1.4)
                bad_sasa = bad_atom_sasa
                return None
        # if ( not bad_sasa is None ):
        #     break



    atomic_depth = core.scoring.atomic_depth.AtomicDepth( pose, 2.3, True, 0.5 )

    type_set = pose.residue(1).type().atom_type_set()

    # Figure out which atoms are buried

    for resnum in range(1, pose.size()+1):
        res = pose.residue(resnum)

        for at in range(1, res.nheavyatoms()+1):

            depth = atomic_depth.calcdepth(res.atom(at), type_set)

            # These are the stupid surface residues and I don't think using total residue sasa is a good idea
            if ( depth < 2.0 ):
                continue

            # for the sasa method, be sure to include the attached hydrogens
            cursasa = surf_vol.surf(resnum,at)
            for hindex in range(res.attached_H_begin( at ), res.attached_H_end( at) + 1):
                cursasa += surf_vol.surf(resnum,hindex)

            if ( cursasa < 0.01 ):
                is_sasa_buried.set(atid(resnum, at), True)

            # for very buried, we'll use 4.5 to allow more

            if ( depth >= 5.5 ):
                is_very_buried.set(atid(resnum, at), True)

            if ( is_sasa_buried.get(atid( resnum, at)) and is_very_buried.get(atid(resnum, at )) ):
                is_both_buried.set(atid(resnum, at), True)


    # get the hbond set
    hb_scorefxn(pose)
    hbset = core.scoring.hbonds.HBondSet()
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)
    hbset.hbond_options().bb_donor_acceptor_check(False)
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)


    # build salt bridges

    positive = ["ARG", "LYS"]
    negative = ["ASP", "GLU"]

    minimum_hb_cut = -1.0
    salt_scorefxn( pose )

    is_salt_bridged = core.id.AtomID_Map_bool_t()
    is_salt_bridged.resize(pose.size())
    for i in range(1, pose.size()+1):
        is_salt_bridged.resize( i, pose.residue(i).natoms(), False)


    graph = pose.energies().energy_graph()
    for i in range(1, pose.size()+1):
        node = graph.get_node(i)
        iterr = node.const_upper_edge_list_begin()
        while (iterr.valid()):
            edge = iterr.__mul__()
            iterr.plus_plus()
            j = edge.get_other_ind(i)
        # for j in range(i+1, pose.size()+1):
            # if ( i >= j ):
            #     continue
            # edge = graph.find_edge(i, j)
            if ( edge is None ):
                continue

            name31 = pose.residue(i).name3()
            name32 = pose.residue(j).name3()
            valid = False
            if ( name31 in positive and name32 in negative):
                valid = True
            if ( name31 in negative and name32 in positive):
                valid = True
            if ( not valid ):
                continue

            edge = graph.find_edge(i, j)
            if ( edge is None ):
                continue
            emap = edge.fill_energy_map()
            score = emap[core.scoring.hbond_sc]
            if ( score < minimum_hb_cut ):
                for atno in range(5, pose.residue(i).nheavyatoms()+1):
                    is_salt_bridged.set(atid(i, atno), True)
                for atno in range(5, pose.residue(j).nheavyatoms()+1):
                    is_salt_bridged.set(atid(j, atno), True)

    # use the elife program to get the pi-pi interactions
    dump_pose = pose.clone()
    dump_pose.pdb_info( core.pose.PDBInfo( dump_pose ) )
    dump_pose.dump_pdb("tmp.pdb")
    # results = cmd("OPENBLAS_NUM_THREADS=1 /home/bcov/sc/scaffold_comparison/tools/elife_write_PDBcontacts.py tmp.pdb")
    # print(results)
    # results = results.split("\n")
    results = get_pi_pi_contacts("tmp.pdb")


    # Setup the atom maps:

    is_amide_stack = core.id.AtomID_Map_bool_t()
    is_amide_stack.resize(pose.size())
    is_ring_stack = core.id.AtomID_Map_bool_t()
    is_ring_stack.resize(pose.size())
    is_arg_stack = core.id.AtomID_Map_bool_t()
    is_arg_stack.resize(pose.size())
    is_carb_stack = core.id.AtomID_Map_bool_t()
    is_carb_stack.resize(pose.size())
    for i in range(1, pose.size()+1):
        is_amide_stack.resize( i, pose.residue(i).natoms(), False)
        is_ring_stack.resize( i, pose.residue(i).natoms(), False)
        is_arg_stack.resize( i, pose.residue(i).natoms(), False)
        is_carb_stack.resize( i, pose.residue(i).natoms(), False)

    # We are going to store part2's type into part 1
    def store_interaction_oneway( part1, part2, pose, is_amide_stack, is_ring_stack, is_arg_stack, is_carb_stack ):
        res1 = int(part1[1])
        res2 = int(part2[1])

        store_map = None
        if (len(part2[2]) == 6):
            store_map = is_amide_stack
        else:
            arg = ["ARG"]
            amides = ["ASN", "GLN"]
            carbs = ["ASP", "GLU"]
            rings = ["TRP", "TYR", "PHE", "HIS"]
            if ( part2[2] in arg ):
                store_map = is_arg_stack
            if ( part2[2] in amides ):
                store_map = is_amide_stack
            if ( part2[2] in carbs ):
                store_map = is_carb_stack
            if ( part2[2] in rings ):
                store_map = is_ring_stack
        if ( store_map is None ):
            print("Error! Unrecognized interaction type: " + ".".join(part2))
            sys.exit(1)

        if ( len(part1[2]) == 6 ):
            store_map.set(atid(res1, 3), True)
            store_map.set(atid(res1, 4), True)
            store_map.set(atid(res1+1, 1), True)
        else:
            for atno in range(5, pose.residue(res1).nheavyatoms()+1):
                store_map.set(atid(res1, atno), True)


    for result in results:
        if ("*" not in result):
            continue
        parts = result.split("*")

        parts[0] = parts[0].split(".")
        parts[1] = parts[1].split(".")

        store_interaction_oneway( parts[0], parts[1], pose, is_amide_stack, is_ring_stack, is_arg_stack, is_carb_stack )
        store_interaction_oneway( parts[1], parts[0], pose, is_amide_stack, is_ring_stack, is_arg_stack, is_carb_stack )



    # It's time to start classifying
    modifiers = {
        "amide-stack":is_amide_stack,
        "ring-stack":is_ring_stack,
        "arg-stack":is_arg_stack,
        "carb-stack":is_carb_stack,
        "salt-bridged":is_salt_bridged,
        "-":None
    }

    burials = {
        "sasa-buried":is_sasa_buried,
        "very-buried":is_very_buried,
        "both-buried":is_both_buried
    }

    atom_types = {
        "ARG-N":["NH1", "NH2"],
        "ARG-NE":["NE"],
        "LYS-N":["NZ"],
        "BB-N":["N"],
        "BB-O":["O"],
        "AMIDE-N":["ND2", "NE2"],
        "AMIDE-O":["OD1", "OE1"],
        "CARB-O":["OD1", "OD2", "OE1", "OE2"],
        "TYR-O":["OH"],
        "TYR-H":["HH"],
        "SER-O":["OG", "OG1"],
        "SER-H":["HG", "HG1"],
        "TRP-N":["NE1"],
        "HIS-N":["ND1", "NE2"]
    }

    groups = {
        "G-ARG":["ARG"],
        "G-AMIDE":["ASN", "GLN"],
        "G-BB":[],
        "G-CARB":["ASP", "GLU"],
        "G-HIS":["HIS"],
        "G-SER":["SER", "THR"],
        "G-TYR":["TYR"]
    }


    def get_max_atom_hb(atom_type):
        if (atom_type == "LYS-N"):
            return 5
        else:
            return 4

    def get_max_group_hb(group):
        if (group == "G-ARG" or group == "G-CARB"):
            return 6
        else:
            return 5

    # Prepare the scorefile ahead of time so that they're all the same

    salt_bridge_terms = ["ARG", "LYS", "CARB"]

    for modifier in modifiers:
        for burial in burials:

            for atom_type in atom_types:
                name = atom_type + "_" + modifier + "_" + burial

                if ( modifier == "salt-bridged" ):
                    valid = False
                    for term in salt_bridge_terms:
                        if ( term in atom_type ):
                            valid = True
                    if ( not valid ):
                        continue

                if ( "SER" in atom_type and modifier != "-"):
                    continue

                for i in range(0, get_max_atom_hb(atom_type) + 1):
                    out_score_map[name + "_hb%i"%i] = 0
                    out_string_map[name + "_hb%i_"%i] = ","


            for group in groups:
                name = group + "_" + modifier + "_" + burial

                if ( modifier == "salt-bridged" ):
                    valid = False
                    for term in salt_bridge_terms:
                        if ( term in group ):
                            valid = True
                    if ( not valid ):
                        continue

                if ( "SER" in group and modifier != "-"):
                    continue

                for i in range(0, get_max_group_hb(group) + 1):
                    out_score_map[name + "_hb%i"%i] = 0
                    out_string_map[name + "_hb%i_"%i] = ","


    def get_num_hbonds(seqpos, atno):
        num_hbonds = 0
        for hb in hbset.atom_hbonds(atid(seqpos, atno), False):
            if ( hb.energy() < -0.25 ):
                num_hbonds += 1
        return num_hbonds

    def register_atom_hbonds( seqpos, atno, atom_type):
        res = pose.residue(seqpos)
        if ( res.type().atom_type(atno).element() == "N" and res.heavyatom_has_polar_hydrogens(atno)):
            num_hbonds = 0
            for hat in range(res.attached_H_begin(atno), res.attached_H_end(atno)+1):
                num_hbonds += get_num_hbonds(seqpos, hat)
        else:
            num_hbonds = get_num_hbonds(seqpos, atno)

        num_hbonds = min(num_hbonds, get_max_atom_hb(atom_type))

        use_atno = atno
        if ( pose.residue(seqpos).type().atom_is_hydrogen( atno ) ):
            use_atno = pose.residue(seqpos).type().atom_base( atno )

        any_modifier = False
        for modifier in reversed(sorted(modifiers.keys())):
            if ( modifier == "-"):
                if ( any_modifier ):
                    continue
            else:
                is_modifier = modifiers[modifier]
                if ( not is_modifier.get(atid(seqpos, use_atno )) ):
                    continue
                any_modifier = True

            for burial in burials:
                is_buried = burials[burial]
                if ( is_buried.get(atid( seqpos, use_atno )) ):
                    name = atom_type + "_" + modifier + "_" + burial + "_hb%i"%num_hbonds
                    out_score_map[name] += 1
                    out_string_map[name + "_"] += "%i/%i-%s,"%(seqpos, pose.pdb_info().number(seqpos),
                                                                pose.residue(seqpos).atom_name(atno).strip())


    def register_group_hbonds( atom_pairs, group_type):
        num_hbonds = 0
        for seqpos, atno in atom_pairs:
            num_hbonds += get_num_hbonds(seqpos, atno)
        num_hbonds = min(num_hbonds, get_max_group_hb(group_type))

        any_modifier = False
        for modifier in reversed(sorted(modifiers.keys())):
            if ( modifier == "-"):
                if ( any_modifier ):
                    continue
            else:
                is_modifier = modifiers[modifier]
                matches = False
                for seqpos, atno in atom_pairs:
                    use_atno = atno
                    if ( pose.residue(seqpos).type().atom_is_hydrogen( atno ) ):
                        use_atno = pose.residue(seqpos).type().atom_base( atno )

                    if ( is_modifier.get(atid(seqpos, use_atno )) ):
                        matches = True
                        break
                if ( not matches ):
                    continue
                any_modifier = True

            for burial in burials:
                is_buried = burials[burial]

                all_buried = True
                for seqpos, atno in atom_pairs:
                    use_atno = atno
                    if ( pose.residue(seqpos).type().atom_is_hydrogen( atno ) ):
                        use_atno = pose.residue(seqpos).type().atom_base( atno )

                    if ( not is_buried.get(atid( seqpos, use_atno )) ):
                        all_buried = False
                        break

                if ( all_buried ):
                    name = group_type + "_" + modifier + "_" + burial + "_hb%i"%num_hbonds
                    out_score_map[name] += 1
                    out_string_map[name+"_"] += "%i/%i-%s,"%(seqpos, pose.pdb_info().number(seqpos),
                                                            pose.residue(seqpos).atom_name(atom_pairs[0][1]).strip())


    def get_sc_polars(seqpos):
        pairs = []
        res = pose.residue(seqpos)
        for at in range(5, res.nheavyatoms()+1):
            if ( res.heavyatom_is_an_acceptor(at) ):
                pairs.append((seqpos, at))
            if ( res.heavyatom_has_polar_hydrogens(at) ):
                for hat in range(res.attached_H_begin(at), res.attached_H_end(at)+1):
                    pairs.append((seqpos, hat))
        return pairs

    # This is where all the atoms are classified

    for seqpos in range(1, pose.size() + 1):
        res = pose.residue(seqpos)
        name3 = res.name3()
        tp = res.type()

        # Deal with bb amides
        if ( seqpos > 1 and res.has("H") ):
            register_atom_hbonds( seqpos, res.atom_index("H"), "BB-N")

        if ( seqpos < pose.size() and res.has("O") ):
            register_atom_hbonds( seqpos, res.atom_index("O"), "BB-O")

            if ( pose.residue(seqpos+1).has("H") ):
                register_group_hbonds( [(seqpos, res.atom_index("O")),
                                        (seqpos+1, pose.residue(seqpos+1).atom_index("H"))], "G-BB" )


        # Look through all the atoms
        for at in range(5, res.natoms()+1):
            if ( not ( tp.atom_is_polar_hydrogen(at) or res.heavyatom_is_an_acceptor(at) or res.heavyatom_has_polar_hydrogens(at) ) ):
                continue
            atname = res.atom_name(at).strip()
            for atom_type in atom_types:

                # believe it or not, there's only three collisions doing this
                if ( atname in atom_types[atom_type] ):
                    if ( atname == "NE2" ):
                        if ( name3 == "HIS" and atom_type == "AMIDE-N"):
                            continue
                        if ( name3 == "GLN" and atom_type == "HIS-N" ):
                            continue
                    if ( atname == "OD1" ):
                        if ( name3 == "ASN" and atom_type == "CARB-O"):
                            continue
                        if ( name3 == "ASP" and atom_type == "AMIDE-O" ):
                            continue
                    if ( atname == "OE1" ):
                        if ( name3 == "GLN" and atom_type == "CARB-O"):
                            continue
                        if ( name3 == "GLU" and atom_type == "AMIDE-O" ):
                            continue

                    register_atom_hbonds( seqpos, at, atom_type )

        # Identify relevant groups
        for group in groups:
            if ( name3 in groups[group] ):
                atom_pairs = get_sc_polars(seqpos)
                register_group_hbonds(atom_pairs, group)


    # debugging
    the_locals = locals()

    return None




############### BEGIN MAIN FUNCTION ###########################

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    pdbs = silent_tools.get_silent_index(silent)["tags"]

    sfd_out = core.io.silent.SilentFileData( "out.silent", False, False, "binary", core.io.silent.SilentFileOptions())

starting_tags = []
if ( os.path.exists("score.sc") ):
    sfd = core.io.raw_data.ScoreFileData("score.sc")
    starting_tags = list(sfd.read_tags_fast("score.sc"))

num = -1
for pdb in pdbs:
    t0 = time.time()
    print("Attempting pose: " + pdb)

    # try:
    for k in [1]:

        name_no_suffix = my_rstrip(my_rstrip(os.path.basename(pdb), ".gz"), ".pdb")
        if ( not args.overwrite and name_no_suffix in starting_tags ):
            print("Already done! Skipping")
            continue

        if ( silent == "" ):
            pose = pose_from_file(pdb)
        else:
            pose = Pose()
            sfd_in.get_structure(pdb).fill_pose(pose)


        sfd = core.io.raw_data.ScoreFileData("score.sc")

        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()


        out_pose = types_of_buried_polars(pose, score_map, string_map, "")

        print("Buried polars found:")
        for name, value in list(string_map.items()):
            if ( len(value) > 1):
                print("%40s: %s"%(name[:-1], value[1:]))

        core.io.raw_data.ScoreMap.add_arbitrary_score_data_from_pose( pose, score_map)
        core.io.raw_data.ScoreMap.add_arbitrary_string_data_from_pose( pose, string_map)

        sfd.write_pose( pose, score_map, name_no_suffix, string_map)
        if (out_pose != None):


            # pdb_info = core.pose.PDBInfo(pose)
            # pose.pdb_info(pdb_info)
            if ( silent == "" ):
                out_pose.dump_pdb(name_no_suffix + ".pdb")
            else:
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(out_pose, pdb)
                sfd_out.add_structure(struct)


        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



if ( silent != "" ):
    sfd_out.write_all("out.silent", False)
