Here is your updated, clean, GitHub-ready section for **Step 1: Preparing the Input PDB in PyMOL**.

---

# Step 1 — Preparing the Input PDB (PyMOL Preprocessing)

Before running **AI2BMD**, the input structure should:

* Have all hydrogens added
* Have properly capped N- and C-termini (if required)
* Be saved as a clean, simulation-ready PDB

This preprocessing is done **inside PyMOL**, not standard Python.

---

## 🔧 Install PyMOL (Conda)

If not already installed:

```bash
conda install -c conda-forge pymol-open-source
```

Launch PyMOL:

```bash
pymol
```

---

## 🧬 Example: Preparing Ubiquitin (1UBQ)

Structure source:
RCSB PDB → [https://www.rcsb.org/structure/1UBQ](https://www.rcsb.org/structure/1UBQ)

---

## 📌 Important

Run the following commands **inside the PyMOL command line**, not in a Python interpreter.

---

## 📝 PyMOL Commands

```python
# Load structure
cmd.load("1UBQ.pdb", "molecule")

# Add all hydrogens
cmd.h_add("molecule")

# Start mutagenesis wizard (used for terminal capping)
cmd.wizard("mutagenesis")

# Set N-terminal cap to acetyl (ACE)
cmd.get_wizard().set_n_cap("acet")
# Manually click the N-terminal residue in the PyMOL viewer
cmd.get_wizard().apply()

# Set C-terminal cap to N-methylamide (NME)
cmd.get_wizard().set_c_cap("nmet")
# Manually click the C-terminal residue in the PyMOL viewer
cmd.get_wizard().apply()

# Exit wizard
cmd.set_wizard()
```

---

## 🔎 What This Does

| Step                | Purpose                                  |
| ------------------- | ---------------------------------------- |
| `h_add()`           | Adds missing hydrogen atoms              |
| `set_n_cap("acet")` | Caps N-terminus with acetyl (ACE)        |
| `set_c_cap("nmet")` | Caps C-terminus with N-methylamide (NME) |

### Why Cap Termini?

Many simulation methods (including AI2BMD) assume:

* Neutral termini
* No artificial terminal charges

Capping prevents:

* Spurious electrostatics
* Edge artifacts
* Unphysical terminal interactions

---

## 💾 Save the Processed Structure

After modifications:

```python
cmd.save("1UBQ_preprocessed.pdb", "molecule")
```

This file will be used as input for AI2BMD.

---

## ⚠️ Notes

* Ensure there are **no alternate locations (ALTLOC)** remaining.
* Remove crystallographic waters unless explicitly needed.
* Verify no missing residues.
* Check protonation states if working at specific pH.

Optional cleanup:

```python
remove solvent
remove resn HOH
```

---

## ✅ Final Output

```
1UBQ_preprocessed.pdb
```


