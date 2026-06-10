# ChikaterasuGUI

A Qt C++ GUI front-end for the [Chikaterasu](https://github.com/yurayura-nmr/chikaterasu) GROMACS automation script.  
Designed to lower the barrier for students setting up MD simulations — no manual script editing required.

---

## Requirements

- GROMACS built from source, installed to `/usr/local/gromacs`
- Tested currently only for gromacs 2025 and not all functions tested
- Qt 5.15+ or Qt 6.x (Widgets module)
- CMake 3.16+
- `awk` (used for nsteps calculation — standard on macOS/Linux)

---

## Build

```bash
cd gui
mkdir build && cd build
cmake ..
make -j
```

Install Qt if needed:

```bash
# macOS (Homebrew)
brew install qt

# Ubuntu / Debian
sudo apt install qt6-base-dev cmake build-essential
```

---

## Changes to chikaterasu.sh due to GUI development

Because the `source` line comes after the manual defaults, GUI values silently override them.  
The script remains fully usable from the terminal without any other changes.

The GUI sets a `CHIKA_GUI=1` environment variable. The script uses this to select
non-interactive GROMACS command paths (see below).

---

## Parameters exposed

### Simulation

| Widget | `chika_gui.conf` variable | Notes |
|--------|--------------------------|-------|
| Molecule / PDB name | `protein_name` | Must match a `.pdb` file in `gromacs/coord/` |
| Force field | `forcefield` | Passed as `-ff` to `pdb2gmx`; skips interactive menu |
| Temperature | `ref_t` | Patched into `ref_t` in `nvt.mdp`, `npt.mdp`, and `md.mdp` |
| Simulation time | `sim_time_ns` | Converted to `nsteps` in `md.mdp` via `awk` |
| Number of runs | `nruns` | 1 for testing, up to 20 for high-reproducibility requiring studies |
| Debug / stop after | `debug_level` | 0 = full run; 1–6 = stop after a specific stage |

### Topology

| Widget | `chika_gui.conf` variable | Notes |
|--------|--------------------------|-------|
| Histidine protonation | `his_manual` | `true` adds `-his` flag to `pdb2gmx` |
| Disulfide bonds | `disulfide` | `true` adds `-ss` (auto-detect); interactive mode CLI-only |

### Ions

| Widget | `chika_gui.conf` variable | Notes |
|--------|--------------------------|-------|
| Ion mode | `specify_salt_concentration` | Toggle between concentration and manual count |
| Salt concentration | `salt_concentration` | Molar; used with `gmx genion -conc` |
| Positive ions | `pos_ions` | Manual mode only |
| Negative ions | `neg_ions` | Manual mode only |
| Mg²⁺ instead of Na⁺ | `magnesium` | Applies to both ion modes. Useful if ATP etc. expect specific counter ions. |

### MDP file patching

The GUI values for temperature and simulation time are patched directly into the
`chika_mdp/` files at runtime by `chikaterasu.sh` using `sed`:

```bash
# In chikaterasu.sh, after sourcing chika_gui.conf:
REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
MDP_DIR="$REPO_ROOT/chika_mdp"

sed -i '' "s/^ref_t[[:space:]]*=.*/ref_t   = $ref_t   $ref_t/" "$MDP_DIR/nvt.mdp"
sed -i '' "s/^ref_t[[:space:]]*=.*/ref_t   = $ref_t   $ref_t/" "$MDP_DIR/md.mdp"

nsteps=$(awk "BEGIN {printf \"%d\", $sim_time_ns * 500000}")
sed -i '' "s/^nsteps[[:space:]]*=.*/nsteps  = $nsteps/"        "$MDP_DIR/md.mdp"
```

`dt = 0.002 ps` is assumed (standard AMBER-ILDN setup), giving 500,000 steps per ns.

---

## Interactive vs GUI mode

Several GROMACS commands (`pdb2gmx -chainsep interactive`, `-merge interactive`, `-ter`)
require stdin and cannot run inside the GUI's output log. The script detects GUI mode
via the `CHIKA_GUI` environment variable and switches to a non-interactive path:

```bash
if [ -n "$CHIKA_GUI" ]; then
    # non-interactive: force field pre-selected via -ff, no -ter, no -merge
    gmx pdb2gmx ... -ff $forcefield -ignh
else
    # original interactive CLI path unchanged
    gmx pdb2gmx ... -chainsep interactive -merge interactive -ter
fi
```

This means **the script behaves identically from the terminal** — the GUI path is purely additive.

---

## Persistent settings

The last used script path is remembered across sessions via `QSettings`.  

---

## Extending the GUI

Every new parameter follows the same four-step pattern:

1. Add a member widget pointer in `mainwindow.h`
2. Construct and lay it out inside `setupUI()` in `mainwindow.cpp`
3. Write its value to `chika_gui.conf` inside `buildAndWriteConfig()`
4. Read that variable in `chikaterasu.sh` (it is already sourced)

### Parameters not yet in the GUI (planned)

- Box shape and size (`cell_shape`, `box_dim`, `box_manual`)
- Water model (`water` — tip3p / tip4p / tip5p / spce)
- Small molecule insertion (`insert_small_molecules`, `protein_added_small_molecule_name`)
- Distance restraints for Zn²⁺ (`distance_restraints`)
- AMBER force field toggle (`amber`)
- Progress bar tied to debug level stages
- Plot button invoking `gmx energy` on completed output
