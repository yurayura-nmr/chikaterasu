# Launching AI2BMD Simulations

After preparing your input PDB (`ub.pdb`) and generating the preprocessing directory (`ub_preprocessed/`), you are ready to launch molecular dynamics.

---

## Expected Directory Structure

Your working directory should look like this:

```
project_folder/
│
├── ub.pdb
├── ub_preprocessed/
│   ├── ub-preeq.pdb
│   └── ub-preeq-nowat.pdb
└── ai2bmd
```

### File Descriptions

| File                 | Purpose                                           |
| -------------------- | ------------------------------------------------- |
| `ub.pdb`             | Original input structure                          |
| `ub_preprocessed/`   | Generated preprocessing directory                 |
| `ub-preeq.pdb`       | Pre-equilibrated structure (with solvent if used) |
| `ub-preeq-nowat.pdb` | Pre-equilibrated structure (no water)             |

---

## Short Test Run (≈ 1 ps)

Use this for quick validation that everything runs correctly.

```bash
./ai2bmd \
  --prot-file ub.pdb \
  --preprocess-dir ub_preprocessed \
  --preeq-steps 0 \
  --sim-steps 1000 \
  --record-per-steps 1
```

### What This Does

* `--prot-file` → Input protein structure
* `--preprocess-dir` → Directory generated during preprocessing
* `--preeq-steps 0` → Skip additional pre-equilibration
* `--sim-steps 1000` → Run 1000 MD steps (~1 ps depending on timestep)
* `--record-per-steps 1` → Save every step (high-resolution output)

Use this to:

* Confirm no crashes
* Check energy behavior
* Validate trajectory output

---

## Production Run (≈ 10 ns)

Once the short test run succeeds, launch a longer simulation:

```bash
./ai2bmd \
  --prot-file ub.pdb \
  --preprocess-dir ub_preprocessed \
  --preeq-steps 0 \
  --sim-steps 10000000 \
  --record-per-steps 1000
```

### Key Differences

* `--sim-steps 10000000` → Longer production trajectory (~10 ns depending on timestep)
* `--record-per-steps 1000` → Save every 1000 steps (reduces file size)

