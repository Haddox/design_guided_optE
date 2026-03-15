# CLAUDE.md

## Project summary

Snakemake pipeline that relaxes and redesigns 20 native protein structures under two Rosetta energy functions (`beta_nov16_cart`, `beta_jan25_cart`), then compares interatomic distance distributions against crystal structure references.

## Pipeline rules (in dependency order)

1. `relax` - Rosetta cartesian relaxation of each native PDB (per energy function)
2. `design` - Rosetta FastDesign of each native PDB (per energy function)
3. `gen_ref` - Generate reference atom-pair distance distributions from crystal structures
4. `compute_test_distributions` - Compare relaxed structures against crystal reference
5. `compute_design_distributions` - Compare designed structures against crystal reference

## Key files

- `Snakefile` - Pipeline definition
- `config.yaml` - All configurable paths, PDB stems, energy functions
- `scripts/distdstr_0.3.py` - Local copy of distance distribution script, patched for NumPy 2.x (`numpy.str_` instead of removed `numpy.unicode_`)
- `scripts/fast_design_with_surface_and_ss_aa_comp_3.xml` - RosettaScripts design protocol
- `data/` - Symlink to `../atom_pair_benchmark/data` (required by distdstr at runtime)

## Running the pipeline

```bash
snakemake --cores <N>
```

Relax and design jobs are complete. Remaining jobs: `gen_ref` (done), `compute_test_distributions` (x2), `compute_design_distributions` (x2).

## Notes

- The `data/` symlink is required because `distdstr_0.3.py` expects `data/new_atom_typing/atom_properties.txt` as a relative path.
- Input PDBs come from `../atom_pair_benchmark/pdbs_for_validation/`.