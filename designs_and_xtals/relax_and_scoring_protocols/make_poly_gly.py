import pyrosetta
import pyrosetta.rosetta
init_flags = ' '.join([
    '-mute all', '-mute core', '-mute protocols',
    '-detect_disulf False', '-read_only_ATOM_entries',
])
pyrosetta.init(init_flags)

def make_poly_gly(input_pdb, output_pdb):
    """
    Use PyRosetta to convert an input PDB to a poly-Gly backbone
    """

    # Read in pose
    pose = pyrosetta.pose_from_pdb(input_pdb)

    # Mutate each residue to glycine
    for res_i in range(1, pose.size()+1):
        pyrosetta.toolbox.mutants.mutate_residue(
            pose, res_i, 'G'
        )

    # Save new pose to output file
    if not os.path.isfile(output_pdb):
        pose.dump_pdb(output_pdb)
