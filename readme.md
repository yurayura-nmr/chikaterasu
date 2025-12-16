# Chikaterasu

**C**ommand-line **H**elper for **I**nstalling **K**ernel-**A**ccelerated **T**hermodynamic **E**xploration **R**uns **A**nd **S**etting **U**p

<img src="logo.png" alt="Chikaterasu logo" align="right" />

Chikaterasu is a Bash script designed to automate the setup of multiple identical molecular dynamics (MD) simulations. 
Current parameters are taken for the amber-type force fields like amber99sb-ildn and amber03ws.

## Requirements

Chikaterasu requires GROMACS, and any recent version should be compatible. 
To install GROMACS, download, extract, and build from source with GPU support.
Verify installation by running:

```bash
gmx --version  # e.g., after installing running the commands in misc/fresh_setup_linux_pc.sh
```

Ensure the GPU and CUDA are recognized correctly.

## Usage

Clone the repository:

```bash
git clone https://github.com/yurayura-nmr/chikaterasu.git
```

### Setup Steps

1. Place your system's PDB file (e.g., 1UBQ.pdb) in `./gromacs/coord/`. [e.g. 1UBQ.pdb]
2. Edit the filename in the top section of `chikaterasu.sh` (e.g., 1UBQ). [e.g. 1UBQ]
3. Configure MD parameters in `chika_mdp/*.mdp` files, setting values for NVT, NPT, and production MD stages (e.g., time, temperature, frame frequency).
4. For new systems, test each level incrementally (1-6) to confirm successful MD preparation and execution. For familiar systems, a full run can be initiated at level 0.

### Execution Levels

Use the following commands to run Chikaterasu at various levels:

**Level 1**: Checks PDB to GROMACS conversion:

```bash
sh chikaterasu.sh 1
```

**Level 2**: Tests protein solvation:

```bash
sh chikaterasu.sh 2
```

**Level 3 - Level 6**: Sequential tests for preparation, ending with equilibration stages.

**Production (Level 0)**: Start production MD:

```bash
sh chikaterasu.sh 0
```

---

## Special bonds

### Metal coordination

Metal coordination involves forming partially covalent interactions between Lewis acids such as Zn²⁺ and Lewis bases such as the thiolate groups of cysteine residues. 
The lone pair on sulfur interacts strongly with the metal center, and in practice this interaction is stiff enough that it behaves like a short bond during MD.

The example below describes how to introduce a tetrahedrally coordinated Zn²⁺ ion in the protein. 
Although this example uses Zn²⁺ and cysteines, the same approach applies to other divalent ions (e.g., Ca²⁺, Mg²⁺) and their coordinating residues.

First, choose the cysteine residues that will coordinate the Zn²⁺ ion and define the deprotonated ones as **CYM** in the PDB file. 
A typical Zn²⁺ site uses two CYM residues (anionic thiolates) for charge balance, together with two neutral CYS residues.
For this example, let us assume:

* CYM: Cys197, Cys214
* CYS: Cys200, Cys211

Run:

```bash
./chikaterasu.sh 1
```

The generated topology is placed in `gromacs/top`. From there, identify:

* the atom index of the Zn²⁺ ion (e.g., 3358), and
* the Sγ atoms of the four coordinating cysteine residues (e.g., 2564, 2595, 2768, 2811).

These numbers are used to define the metal–ligand bonds.

### Coordination bonds

In this workflow, **the coordination geometry is enforced by explicit bonds**, not distance restraints. 
Our tests show that Zn²⁺ stays stably coordinated even if the `DISRES` flag is absent, because these covalent-style bonds dominate the interaction.

Create a file called `zn_bond.top` containing the Zn–Sγ bonds:

```
2564  3358     6 0.24  4000    ; CYM197 SG
2595  3358     6 0.24  4000    ; CYS200 SG
2768  3358     6 0.24  4000    ; CYS211 SG
2811  3358     6 0.24  4000    ; CYM214 SG
```

These use a stiff harmonic bond (type 6) at 0.24 nm, which reliably holds the Zn²⁺ coordination sphere together.

Chikaterasu automatically inserts `zn_bond.top` into the final `topol.top` right before the `[ pairs ]` section, so no manual editing is needed.

### Optional distance restraints

Distance restraints can still be added for experiments or for enforcing a particular coordination geometry, but they are **not required** for stability when these bonds are defined. 
If you choose to use them, place them in `distance_restraints.itp` and enable them with `define = -DDISRES` in your `md.mdp` file. 
Otherwise, they are simply ignored.

Chikaterasu prints reminder messages indicating whether:

* `-DDISRES` is present in the MDP (if you intend to use restraints), and
* the Zn²⁺–protein bonds are correctly included in the final topology.

Future updates will expand these checks, but the core metal-coordination topology is now fully automated.
