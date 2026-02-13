# Using AI2BMD Trajectories (with ASE, MDAnalysis, and PyMOL)

This guide documents how to:

* Install required tools (`ase`, `MDAnalysis`)
* Convert `AI2BMD` trajectory files (`.traj`) into `.dcd`
* Visualize trajectories in PyMOL
* Use `.dcd` files for downstream analysis with MDAnalysis

---

## Install Required Python Packages (Conda)

Create or activate your conda environment, then install:

```bash
# Install core packages
conda install -c conda-forge ase mdanalysis -y
```

Optional (recommended for visualization workflows):

```bash
conda install -c conda-forge pymol-open-source -y
```

---

## Convert `.traj` to `.dcd` (for PyMOL, MDanalysis)

AI2BMD outputs trajectories in ASE `.traj` format.
For easier visualization in PyMOL and compatibility with many MD tools, convert to `.dcd`.

Example:

```bash
python traj2dcd.py --input ../../Logs-ub/ub-traj.traj --output test.dcd --pdb ../../ub_preprocessed/ub-preeq.pdb --stride 1
```

### Argument Explanation

* `--input` → ASE trajectory file (`.traj`)
* `--output` → Output trajectory file in DCD format
* `--pdb` → Static structure file used as topology reference
* `--stride` → Frame skipping (1 = use every frame)

---

## Load into PyMOL

To correctly visualize:

1. Load the **static PDB first**
2. Then load the **DCD file**

Order matters because:

* The PDB provides atom topology
* The DCD provides coordinates only

In PyMOL:

```python
load chig-preeq.pdb
load_traj test.dcd
```

---

## Using DCD with MDAnalysis

DCD is a standard trajectory format originally from CHARMM/NAMD and is widely supported.

Example usage in MDAnalysis:

```python
import MDAnalysis as mda

u = mda.Universe("chig-preeq.pdb", "test.dcd")

print(len(u.trajectory))   # number of frames
print(len(u.atoms))        # number of atoms
```

---

| Format  | Best For               | Notes                |
| ------- | ---------------------- | -------------------- |
| `.traj` | ASE workflows          | Native AI2BMD output |
| `.dcd`  | PyMOL, MDAnalysis, VMD | More interoperable   |
| `.pdb`  | Static topology        | Required for DCD     |

---

## Recommended Workflow Summary

AI2BMD simulation → `.traj`
⬇
Convert using `traj2dcd.py`
⬇
Visualize in PyMOL
⬇
Analyze in MDAnalysis


