# pred_rescue: Cartesian ddG pipeline for rescue mutations

Runs Rosetta's `cartesian_ddg` protocol on nonpolar rescue mutations identified by
deep mutational scanning (DMS). For each mutation classified as `large-to-small` or
`small-to-large` (nonpolar size changes), ddG is estimated under two energy functions
(`beta_nov16` and `beta_jan25`).

## Overview

```
DMS libraries (lib1, lib2)
        |
  analyze_dms.py   →   results/rescue_mutations.csv   (32 mutations)
        |
  cartesian_ddg        ×2 energy functions  =  64 ddG jobs
        |
  results/ddg/{energy_func}/{design}/{mutation}/{design}_{mutation}.ddg
```

Each ddG job runs 5 iterations of WT and mutant, writing scores to a `.ddg` file
(10 lines: 5 `WT_` + 5 `MUT_` rounds).

## Requirements

- `pyrosetta` conda environment with `snakemake` and `pandas` installed
- Rosetta `cartesian_ddg` binary (macOS or Linux release)
- PDB files for each design

## Configuration

Edit `config.yaml` to set paths for your system:

```yaml
rosetta_bin:  # path to cartesian_ddg binary
rosetta_db:   # path to Rosetta database
pdb_dir:      # root directory containing rd1/, rd2/, rd4/ subdirs with design PDBs
dms_lib1:     # path to DMS library 1 CSV
dms_lib2:     # path to DMS library 2 CSV
energy_functions:
  - beta_nov16
  - beta_jan25
```

PDB files are expected at `{pdb_dir}/{rd}/{design}.pdb`, where `rd` (e.g. `rd4`) is
parsed from the design name.

## Running

Due to a Snakemake 9 limitation, the pipeline must be invoked in two passes.

**Pass 1** — identify rescue mutations:
```bash
conda run -n pyrosetta snakemake run_analyze_dms --cores 1
```

**Pass 2** — run all ddG calculations:
```bash
conda run -n pyrosetta snakemake "results/ddg/.done" --cores N
```

Replace `N` with the number of parallel jobs desired. Each job takes ~2 minutes on a
modern Mac; with `--cores 4` the full set of 64 jobs takes ~30 minutes.

If restarting after an interruption, add `--rerun-incomplete`:
```bash
conda run -n pyrosetta snakemake "results/ddg/.done" --cores N --rerun-incomplete
```

## Output

```
results/
  rescue_mutations.csv              — all rescue mutations with size-change annotation
  ddg/
    {energy_func}/
      {design}/
        {mutation}/
          {design}_{mutation}.mut   — Rosetta mutfile
          {design}_{mutation}.ddg   — ddG scores (5 WT + 5 MUT rounds)
          {design}_{mutation}.sc    — per-round score breakdown
    .done                           — marker file created when all jobs complete
```

## Notes

- The `cartesian_ddg` app in non-legacy mode does not write a `.sc` scorefile;
  the `.ddg` file contains all energy information needed for downstream analysis.
- Each job `cd`s into its output directory before running Rosetta so that crash
  logs and intermediate files are isolated per mutation.
- Weights `cart_bonded 0.5 pro_close 0` are required for Cartesian minimization.
