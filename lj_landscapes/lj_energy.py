import argparse
import math
import os
import pandas as pd

# Parse the score function name from the command line. PyRosetta must be
# initialized with the correct flag before any Rosetta objects are created,
# so argument parsing happens before the import/init block below.
parser = argparse.ArgumentParser(description="Compute LJ energy landscape for an atom pair.")
parser.add_argument("--score_fn", default="beta_jan25", help="Rosetta score function name (e.g. beta_jan25, beta_nov16)")
args = parser.parse_args()

SCORE_FN_NAME = args.score_fn

import pyrosetta
# Pass the score function name as an init flag (e.g. -beta_jan25) to load
# the matching energy method options and etable parameters.
pyrosetta.init(f"-{SCORE_FN_NAME} -mute all")

from pyrosetta.rosetta.core.scoring import ScoringManager, fa_atr, fa_rep
from pyrosetta.rosetta.core.chemical import ChemicalManager

print(SCORE_FN_NAME)

# Create the score function and extract the fa_atr and fa_rep weights.
# These weights scale the raw LJ energies when computing the total score.
sfxn = pyrosetta.create_score_function(SCORE_FN_NAME)
w_atr = sfxn.get_weight(fa_atr)
w_rep = sfxn.get_weight(fa_rep)
print(f"fa_atr weight: {w_atr}, fa_rep weight: {w_rep}")

# Retrieve the etable type string from the score function's energy method options,
# then fetch the corresponding Etable from the ScoringManager. This ensures the
# etable always matches the chosen score function. The ScoringManager returns a
# weak_ptr (CAP), so we call .lock() to get the underlying Etable object.
# The Etable stores precomputed LJ and solvation parameters for all atom type pairs.
etable_type = sfxn.energy_method_options().etable_type()
etable = ScoringManager.get_instance().etable(etable_type).lock()

# Get the fa_standard atom type set, which maps atom type names (e.g. "CH3")
# to integer indices used by the etable.
atom_type_set = ChemicalManager.get_instance().atom_type_set("fa_standard")


def lj_atr_rep(etable, atype1, atype2, d2):
    """
    Compute unweighted LJ attractive and repulsive energies for an atom pair
    at distance squared d2, using Rosetta's analytic etable parameters.

    The split mirrors Rosetta's smooth etable convention:
      - lj_atr: capped at the well minimum for d < d_min, tapers to 0 beyond
      - lj_rep: excess above the well minimum cap for d < d_min, else 0
    """
    # Look up the analytic LJ parameters for this atom type pair.
    # These include the r12/r6 coefficients, well depth/position, and
    # the parameters for the linear ramp and cubic polynomial regions.
    p = etable.analytic_params_for_pair(atype1, atype2)

    # Beyond the interaction cutoff (~6 Å), both terms are zero.
    if d2 > p.maxd2:
        return 0.0, 0.0

    d = math.sqrt(d2)

    # Region 1: very short range — linear ramp.
    # At very close distances the standard LJ diverges, so Rosetta replaces
    # it with a linear ramp to keep energies finite. The attraction is capped
    # at the well minimum; all excess energy is repulsion.
    if d2 < p.ljrep_linear_ramp_d2_cutoff:
        lj_val = p.lj_switch_slope * d + p.lj_switch_intercept
        return p.lj_val_at_minimum, lj_val - p.lj_val_at_minimum

    # Region 2: cubic polynomial taper for the attractive term.
    # Between xlo (~4.5 Å) and xhi (6.0 Å = cutoff), the attractive energy
    # is smoothly tapered to zero using a cubic polynomial, avoiding a
    # discontinuity at the cutoff. No repulsion in this range.
    if p.ljatr_cubic_poly_xlo <= d <= p.ljatr_cubic_poly_xhi:
        cp = p.ljatr_cubic_poly_parameters
        lj_atr = cp.c0 + cp.c1 * d + cp.c2 * d**2 + cp.c3 * d**3
        return lj_atr, 0.0

    # Region 3: standard 12-6 LJ potential.
    # E = lj_r12_coeff / r^12 + lj_r6_coeff / r^6
    # (lj_r6_coeff is negative, giving the attractive well)
    inv_d2  = 1.0 / d2
    inv_d6  = inv_d2 ** 3
    inv_d12 = inv_d6 ** 2
    lj_val = p.lj_r12_coeff * inv_d12 + p.lj_r6_coeff * inv_d6

    if d < p.lj_minimum:
        # Inside the well: attraction is capped at the well minimum energy,
        # and repulsion carries the excess above that cap.
        return p.lj_val_at_minimum, lj_val - p.lj_val_at_minimum
    else:
        # Beyond the well minimum: purely attractive, no repulsion.
        return lj_val, 0.0

# Atom type pairs to evaluate. Each entry is (type1, type2); the distance
# window is centered on the vdW contact distance (sum of LJ radii) for that pair.
PAIRS = [
    ("CH3", "CH3"),
    ("CH3", "CH2"),
    ("CH3", "CH1"),
    ("CH3", "OCbb"),
    # ("CH3", "OH"),
]

n_points = 40
rows = []

for name1, name2 in PAIRS:
    atype1 = atom_type_set.atom_type_index(name1)
    atype2 = atom_type_set.atom_type_index(name2)

    # Center the distance window on the vdW contact distance for this pair.
    r1 = atom_type_set[atype1].lj_radius()
    r2 = atom_type_set[atype2].lj_radius()
    d_center = r1 + r2

    d_lo = d_center - 1.2
    d_hi = d_center + 1.2
    step = (d_hi - d_lo) / (n_points - 1)

    for i in range(n_points):
        dist = d_lo + i * step
        atr, rep = lj_atr_rep(etable, atype1, atype2, dist**2)
        rows.append({
            "pair": f"{name1}:{name2}",
            "d": dist,
            "lj_atr": atr,
            "lj_rep": rep,
            "fa_atr": w_atr * atr,
            "fa_rep": w_rep * rep,
            "total": w_atr * atr + w_rep * rep,
        })

df = pd.DataFrame(rows)
print(df.to_string(index=False))

out_path = f"results/lj_energies_{SCORE_FN_NAME}.csv"
os.makedirs("results", exist_ok=True)
df.to_csv(out_path, index=False)
print(f"\nSaved to {out_path}")
