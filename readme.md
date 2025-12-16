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


## Documentation

More detailed documentation, tutorials, and workflow explanations can be found in our GitHub Wiki:

👉 https://github.com/yurayura-nmr/chikaterasu/wiki

