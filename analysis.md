# `ana_chikaterasu.sh` — Automated MD Analysis for MD runs initiated by Chikaterasu

This bash script automates post-processing of molecular dynamics (MD) simulation results. 
It is meant to simplify the typical workflow after running simulations with GROMACS.

## What it does

The script performs the following steps automatically:

1. **Calculate RMSD** (Root Mean Square Deviation)  
   ➤ Output: `rmsd/rmsd.xvg`

2. **Calculate RMSF** (Root Mean Square Fluctuation)  
   ➤ Output: `rmsf/rmsf.xvg`

3. **Calculate custom distances** (e.g. between specific atoms like C-alpha atoms of M1 and V5)  
   ➤ Output: `distance/distance.xvg`

4. **Prepare trajectory file** for visualization  
   ➤ Output: `traj/traj.pdb`

These output files are automatically saved into appropriately named folders.

## Viewing the results

### On Linux/macOS

Most output is in `.xvg` format (used by GROMACS). You can view them using `xmgrace`, a common plotting tool:

```bash
xmgrace rmsd.xvg
xmgrace rmsf.xvg
xmgrace distance.xvg
```

You must first `cd` into the directory containing the `.xvg` files.

> If you're using macOS, install xmgrace with Homebrew:
> ```
> brew install grace
> ```

### To view the trajectory

Use PyMOL to view the MD trajectory:

```bash
pymol traj/traj.pdb
```

This lets you see how the protein moved during the simulation.

## For further analysis

If you have GROMACS installed, you can do more advanced analysis using `.xtc` trajectory files.



This script is meant to save time and help you quickly extract useful plots and data from your simulations.
