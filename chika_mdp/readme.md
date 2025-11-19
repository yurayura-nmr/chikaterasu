# MD Parameter Files for Reproducible GROMACS Simulations

This subfolder contains the **GROMACS `.mdp` files** used for a standard and reproducible molecular dynamics workflow:

1. **Energy minimization** (`minim.mdp`)
2. **NVT equilibration** (`nvt.mdp`)
3. **NPT equilibration** (`npt.mdp`)
4. **Production MD under NPT** (`md.mdp`, `mdrheoE.mdp`)

---

## 📁 File Overview

### **`minim.mdp`**

Standard steepest-descent energy minimization. Used to relax any steric clashes prior to equilibration.

### **`nvt.mdp`**

Performs **constant-volume (NVT)** equilibration using temperature coupling. 
This step stabilizes system temperature before applying pressure control.

### **`npt.mdp`**

Performs **constant-pressure (NPT)** equilibration. 
This equilibrates system density and volume before production dynamics.

### **`md.mdp`**

Main **production MD** file for long simulations. 
Contains output intervals and integration settings. 
Often modified for specific simulation lengths.

### **`mdrheoE.mdp`**

A variant of the production file, used for rheology-related simulations.

### **`ions.mdp`**

Used for ion addition via `gmx genion`.

---

## 📘 Notes on `md.mdp`

Below is an excerpt explaining key settings:

```ini
; Simulation length
nsteps  = 50000000    ; 100 ns
dt      = 0.002       ; 2 fs (2022/3/6 note: dt=0.0025 may also work)

; Output
nstxout          = 10000   ; coordinates every 20 ps
nstvout          = 10000   ; velocities every 20 ps
nstenergy        = 10000   ; energies every 20 ps
nstlog           = 10000   ; log update every 20 ps
nstxout-compressed = 10000 ; compressed coordinates every 20 ps
compressed-x-grps  = System
```

These values are suitable for **100–200 ns** simulations where trajectory detail matters.

---

## 🧠 Recommendation for Very Long MD Runs (> 1 μs)

For very long trajectories (e.g., **5 μs simulations**), 
the default output frequency can generate massive files and wear down storage devices.

**Recommended high-frequency output settings:**

```ini
nstxout               = 5000000
nstvout               = 5000000
nstenergy             = 5000000
nstlog                = 5000000
nstxout-compressed    = 5000000
```

At 2 fs time step, this saves frames every **~10 ns**, drastically reducing output size:

* A small system typically produces **≈2.5 GB per 1 us** at high frequency.
* With coarse output (5M-step intervals), even **multi-microsecond simulations** remain practical in size.
* It still produces 10 data points every 100 ns, so large-scale domain-domain motion dynamics should be somewhat covered.
* It will produce poor RMSD/RMSF/distance etc. statistics though.

---

The purpose of this collection is to:

* Provide a **consistent, reproducible workflow** for GROMACS MD.
* Allow fast swapping between systems while keeping parameters stable.
* Enable both **short exploratory runs** and **long production simulations**.

Feel free to adapt these files for your own systems while keeping reproducibility in mind.
