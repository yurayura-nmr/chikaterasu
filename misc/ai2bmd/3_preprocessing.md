# Preprocessing with AI2BMD

This step prepares the system for simulation using a classical force field before running AI-driven MD.

⚠️ **Important:**
Preprocessing can be **slow** and may not use the GPU. 
It typically runs Amber's `sander` on CPU, which makes it significantly slower than production runs. 
Avoid repeating this step unless necessary.

---

## Command

```bash
./ai2bmd \
  --prot-file from_step_before.pdb \
  --preprocess-method FF19SB
```

---

## What This Does

* `--prot-file`
  Input structure prepared in the previous step (hydrogens added, termini capped, cleaned).

* `--preprocess-method FF19SB`
  Uses the **AMBER FF19SB force field** for:

  * Topology generation
  * Parameter assignment
  * Pre-equilibration
  * System preparation

This generates a preprocessing directory containing equilibrated structures used for AI2BMD simulations.

---

## Expected Output

After completion, a directory will be created:

```
<protein_name>_preprocessed/
```

Example:

```
ub_preprocessed/
├── ub-preeq.pdb
├── ub-preeq-nowat.pdb
├── topology files
├── coordinate files
└── logs
```

### Important Files

| File                | Purpose                                    |
| ------------------- | ------------------------------------------ |
| `*-preeq.pdb`       | Pre-equilibrated structure                 |
| `*-preeq-nowat.pdb` | Pre-equilibrated structure without solvent |
| Log files           | Useful for diagnosing failures             |

---

## Performance Notes

* This stage often runs **Amber `sander`**
* Typically **CPU-only**
* Can take a long time depending on:

  * System size
  * Solvation
  * Number of equilibration steps
* GPU acceleration is generally **not used** in this phase

Because of this:

✅ Ensure your input PDB is clean before running
✅ Avoid repeating preprocessing unnecessarily
✅ Back up the generated `_preprocessed` directory
