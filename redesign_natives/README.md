# Redesign Natives

Snakemake pipeline that redesigns native protein structures using Rosetta, then evaluates how well each energy function recovers native-like interatomic distance distributions.

## Overview

For each of 20 native crystal structures and two Rosetta energy functions (`beta_nov16_cart`, `beta_jan25_cart`), the pipeline:

1. **Relaxes** each native structure (Rosetta `relax` with cartesian refinement).
2. **Redesigns** each native structure (Rosetta `FastDesign` with surface amino acid composition constraints).
3. **Generates a reference** distribution of interatomic distances from the original crystal structures.
4. **Computes distance distributions** for both relaxed and designed structures, comparing them against the crystal reference.

The goal is to assess whether the energy functions preserve native-like atom-pair distance distributions after relaxation and design.

## Requirements

- Rosetta (relax and rosetta_scripts binaries)
- Python with numpy, scipy, pandas, matplotlib
- Snakemake

## Configuration

All paths and parameters are set in `config.yaml`:

- `rosetta_bin` / `rosetta_scripts_bin`: Paths to Rosetta executables
- `rosetta_db`: Path to Rosetta database
- `pdb_dir`: Directory containing input crystal structure PDBs
- `distdstr_script`: Path to the distance distribution analysis script
- `design_script`: RosettaScripts XML protocol for FastDesign
- `surface_aa_comp_file`: Amino acid composition constraint file for surface residues
- `energy_functions`: List of Rosetta energy functions to compare
- `stems`: List of PDB identifiers to process

## Running

```bash
snakemake --cores <N>
```

## Output

Results are organized under `results/`:

```
results/
  crystal_ref/          # Reference distance distributions from crystal structures
  relax/{ef}/           # Relaxed PDBs per energy function
  design/{ef}/          # Designed PDBs per energy function
  test_distributions/{ef}/    # Distance distributions for relaxed structures
  design_distributions/{ef}/  # Distance distributions for designed structures
```

## Files

- `Snakefile`: Pipeline rules
- `config.yaml`: Configuration (paths, PDB list, energy functions)
- `scripts/distdstr_0.3.py`: Distance distribution analysis script (local copy, patched for NumPy 2.x compatibility)
- `scripts/fast_design_with_surface_and_ss_aa_comp_3.xml`: RosettaScripts XML for FastDesign
- `scripts/total_hydrophobicity.comp`: Surface amino acid composition constraints
- `data/`: Symlink to `../atom_pair_benchmark/data` (atom typing parameters used by distdstr)