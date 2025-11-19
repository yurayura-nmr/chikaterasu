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

## Special bonds

### Metal coordination

Coordination bonds are formed between Lewis acids (e.g., Zn2+) and Lewis bases (e.g., thiol groups R-S-H or R-S(-)). Lone pair electrons on S atom to Zn2+ d-orbital.

The example here is written for the case of Zn2+ ions bound to a protein in the input PDB file. Similar steps should work for other ions (Ca2+, Mg2+, etc.).

Choose 2 coordinating Cys residues and define them as CYM (instead of CYS) in the PDB file to ensure local charge neutrality (Zn2+ + Cym- + Cym-). For this example, let us assume that CYS residues Cys197 and Cys214 were chosen as CYM residues. Cys200 and Cys211 are chosen as regular CYS residues.

Now we need the topology:

```bash
./chikaterasu.sh 1   # this should run normally
```

Now the topology is already accessible in `gromacs/top`. We can see the atom numbers as defined by GROMACS (not the PDB file). Find the atom number of the Zn2+ ion and note it down (let's say, it is 3358). This will be used soon.
Also find the atom numbers of all the coordinating (e.g., S-gamma atoms of Cys residues. Let us say for this example, they are 2564, 2595, 2768, and 2811.

In the chikaterasu base folder, provide an additional file containing the distance restraints for the Zn2+ ion (e.g., to coordinating His or Cys residues). As an example `distance_restraints.itp` file:

```
#ifdef DISRES
[ distance_restraints ]
; ai aj type index type' low up1 up2 fac
  2564 3358   1   1     1     0.0 0.23 0.3 1      ; CYM197 SG <> Zn2+
  2595 3358   1   2     1     0.0 0.23 0.3 1      ; CYS200 SG <> Zn2+  
  2768 3358   1   3     1     0.0 0.23 0.3 1      ; CYS211 SG <> Zn2+
  2811 3358   1   4     1     0.0 0.23 0.3 1      ; CYM214 SG <> Zn2+
#endif
```

Next, in `chikaterasu.sh`, set the distance restraints flag to true:

```bash
distance_restratints=true
```

Finally, we need some Zn2+-protein bonds. Make another file called "metal_protein_bonds.top" or something like that:

```
2564  3358     6 0.24  4000    ; 16 CYS-SG G-C-T
2595  3358     6 0.24  4000    ; 13 CYS-SG Q-C-P
2768  3358     6 0.24  4000    ; 27 CYS-SG G-C-E
2811  3358     6 0.24  4000    ; 30 CYS-SG G-C-C
```

Now, unfortunately, this file (containing the coordination bonds between the Zn2+ atom and Cys S-gamma atoms) is not yet automatically integrated into final topology file `topol.top`.

At present, Chikaterasu will pause once to give the user the chance to ninja-edit the `topol.top` by adding these special bonds, right before going to energy minimization. (i.e., just copy paste these 4 lines into `topol.top`).

In the future, this will be automated so that the `metal_protein_bonds.top` file is read automatically into `topol.top` and this user-invention is not required.

## Modelling missing loops in the structure

For loop modeling, use UCSF's ModLoop service with your modified PDB file to define loop segments (e.g., 70:A:71:A).
Alternatives are SwissModel, AlphaFold, PyMOL etc.
