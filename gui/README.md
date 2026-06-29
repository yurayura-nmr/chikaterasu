# ChikaterasuGUI

A Qt C++ GUI front-end for the [Chikaterasu](https://github.com/yurayura-nmr/chikaterasu) GROMACS automation script.  
Designed to lower the barrier for students setting up MD simulations — no manual script editing required.

---

## Requirements

- GROMACS built from source, installed to `/usr/local/gromacs`
- Tested currently on Linux with GROMACS 2025/2026; not all functions tested on macOS
- Qt 5.15+ or Qt 6.x (Widgets module)
- CMake 3.16+
- `awk` (used for nsteps calculation)

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

Because the `source` line for `chika_gui.conf` comes after the manual header defaults,
GUI values silently override them. The script remains fully usable from the terminal
without any other changes.

The GUI sets a `CHIKA_GUI=1` environment variable. The script uses this to select
non-interactive GROMACS command paths (see below).

---

## Parameters exposed

### Simulation

| Widget | `chika_gui.conf` variable | Notes |
|--------|--------------------------|-------|
| Molecule / PDB name | `protein_name` | Must match a `.pdb` file in `gromacs/coord/` |
| Force field | `forcefield` | Passed as `-ff` to `pdb2gmx`; skips interactive menu |
| Water model | `water` | Passed as `-water` to `pdb2gmx`; also derives `water_file` (see below) |
| Temperature | `ref_t` | Patched into `ref_t` in `nvt.mdp`, `npt.mdp`, and `md.mdp` |
| Simulation time | `sim_time_ns` | Converted to `nsteps` in `md.mdp` via `awk` |
| Number of runs | `nruns` | 1 for testing, up to 20 for high-reproducibility studies |
| Debug / stop after | `debug_level` | 0 = full run; 1–6 = stop after a specific stage |

### Topology

| Widget | `chika_gui.conf` variable | Notes |
|--------|--------------------------|-------|
| Disulfide bonds | `disulfide` | `true` adds `-ss` (auto-detect); interactive chain/merge selection remains CLI-only |

> **Removed:** manual histidine protonation (`-his`). This GROMACS flag prompts
> interactively per-residue and cannot run inside the GUI's output log — there is
> no safe non-interactive equivalent, so the option was dropped from the GUI.
> Histidine protonation now always uses `pdb2gmx`'s automatic default.

### Box configuration

| Widget | `chika_gui.conf` variable | Notes |
|--------|--------------------------|-------|
| Manually specify box dimensions | `box_manual` | Off by default; reveals dimension/shape fields when checked |
| Box dimensions (nm) | `box_dim` | Free-text `x y z`, only used if `box_manual=true` |
| Cell shape | `cell_shape` | triclinic / cubic / dodecahedron / octahedron |

### Ions

| Widget | `chika_gui.conf` variable | Notes |
|--------|--------------------------|-------|
| Ion mode | `specify_salt_concentration` | Toggle between concentration and manual count |
| Salt concentration | `salt_concentration` | Molar; used with `gmx genion -conc` |
| Positive ions | `pos_ions` | Manual mode only |
| Negative ions | `neg_ions` | Manual mode only |
| Mg²⁺ instead of Na⁺ | `magnesium` | Applies to both ion modes. Useful if ATP etc. expect specific counter ions |

---

## Water model → water file mapping

`chikaterasu.sh` derives `water_file` from the GUI's `water` selection right after sourcing
`chika_gui.conf`.

---

## MDP file patching

The GUI values for temperature and simulation time are patched directly into the
`chika_mdp/` files at runtime by `chikaterasu.sh` using `sed`, **after** sourcing
`chika_gui.conf`.

`dt = 0.002 ps` is assumed (standard AMBER-ILDN setup), giving 500,000 steps per ns.
Note: `sed -i` here uses GNU sed syntax (Linux).

---

## Interactive vs GUI mode

Several GROMACS commands (`pdb2gmx -chainsep interactive`, `-merge interactive`, `-ter`)
require stdin and cannot run inside the GUI's output log. The script detects GUI mode
via the `CHIKA_GUI` environment variable and switches to a non-interactive path:

```bash
if [ -n "$CHIKA_GUI" ]; then
    # non-interactive: force field and water model pre-selected, no -ter, no -merge
    gmx pdb2gmx ... -ff $forcefield -water $water -ignh
else
    # original interactive CLI path unchanged
    gmx pdb2gmx ... -chainsep interactive -merge interactive -ter
fi
```

This means **the script behaves identically from the terminal** — the GUI path is purely additive.

---

## Live progress bar

`mdrun` is launched with `-v`, which writes step progress to stderr in the form:

```
step 15800, will finish Tue ...
```

The GUI parses this directly from the running process's stderr stream (no log-file
polling needed) and updates a progress bar as `current_step / total_steps`:

```cpp
static QRegularExpression reStep(R"(^step\s+(\d+),)");
```

`total_steps` is computed client-side from the simulation time using the same
`* 500000` formula as the bash script's `nsteps` calculation, so the two stay in sync.

---

## Stop button

The Run button is paired with a Stop button that terminates the simulation, including
all child processes (`bash` → `gmx mdrun`). 

This is mainly useful for quick iteration while testing GUI changes, but is also handy
for students who want to abort a misconfigured run early.

---

## Persistent settings

The last used script path is remembered across sessions via `QSettings`.

---

### Parameters not yet in the GUI (planned)

- Shear-flow settings
- Small molecule insertion (`insert_small_molecules`, `protein_added_small_molecule_name`)
- Distance restraints for Zn²⁺ (`distance_restraints`)
