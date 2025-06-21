import pyrosetta
import sys
weights_file = sys.argv[1]
flags_file = sys.argv[2]
pdb = sys.argv[3]

init_flags = ' '.join([
    '-mute all', '-mute core', '-mute protocols',
    '-detect_disulf False', '-read_only_ATOM_entries',
    f'@{flags_file}'
])
pyrosetta.init(init_flags)

sf = pyrosetta.create_score_function(args.weights_f)
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(sf)

pose = pyrosetta.pose_from_pdb(pdb)
print(sf(pose))