#!/bin/sh

: '
*************************************************************
Chikaterasu         Version: dev
gmx                 Version: 2022.3

For recent changes or issues, please refer to the GitHub repository.

Author: Erik Walinda
Affiliation: Kyoto University, Graduate School of Medicine
*************************************************************
'

: '
*************************************************************
Manual Parameter Setup for Current Run
Note: Most parameters are configured via .mdp files.
The starting PDB file must be placed in the gromacs/coord directory.

Debug Level Options:
0   Full production MD run (NPT); no debug output
1   Topology generation (pdb2gmx)
2   Solvation step
3   Addition of counterions and distance restraints
4   Energy minimization
5   NVT equilibration completed
6   NPT equilibration completed
*************************************************************
'

# === Simulation Configuration ===
protein_name="1UBQ"         # Protein or molecule to simulate. For small molecule simulations, use the appropriate file (e.g., ATP.pdb).
nruns=1                     # Number of runs. Use 1 for testing and 10 for production-level simulations.

# === Debugging Options ===
debug_level=0               # Set the debug level for output verbosity. Can be passed as an argument (e.g., ./chikaterasu 0).

# === Histidine Protonation and Zn2+ Options ===
his_manual=false            # Manually specify histidine protonation states (set to true if required).
distance_restraints=false   # Enable distance restraints, typically used for Zn2+ interactions.
disulfide=false             # Specify whether to account for disulfide bridges during topology generation.


# === Ion Configuration ===
specify_salt_concentration=true  # Specify salt concentration in molar units (true), or manually count ions (false).
salt_concentration=0.050         # Desired salt concentration (in mol/L) if specify_salt_concentration=true.

pos_ions=0                  # Number of positive ions (only applicable when manually specifying ion count).
neg_ions=2                  # Number of negative ions (only applicable when manually specifying ion count).

magnesium=false             # If true, use Mg2+ as the positive ion; otherwise, Na+ is used by default.

# === Box Shape and Simulation Setup ===
# Choose one of the following simulation scenarios:
# a) Protein-only simulation:
#    Set insert_small_molecules=false
# b) Small molecule-only simulation:
#    Set insert_small_molecules=true
# c) Protein with small molecules:
#    Set insert_small_molecules=true and protein_with_small_molecules=true

# Define simulation box dimensions for molecule insertion 
# Must set box_dim when working with inserting molecules
insert_small_molecules=false       # If we want to simulate small molecules like ATP or sialic acid
insert_small_molecules_number=2

# Define whether to include a protein with small molecules
protein_with_small_molecules=false
protein_added_small_molecule_name="ATP"  # Name of the small molecule to be inserted alongside the protein (file must be in the gromacs/coord directory)

# Specify box configuration options
amber=false                 # For N-/C-terminus. False for CHARMM36m etc.
box_manual=false            # Enable manual specification of the box size
box_empty=false             # Set to true for a water-only simulation (no solute)

# Manually specified box dimensions (used if box_manual=true). 
# Box dimensions (in nm) should accommodate the solute, especially for rheological MD or hydrodynamics simulations. 
# Ensure no periodic boundary condition (PBC) artifacts.
box_dim="    7.845   6.497   6.363 "  
cell_shape="triclinic"      # Available shapes: triclinic, cubic, dodecahedron, octahedron

# --- Water Model Configuration ---
# Specify the water model to be used in the simulation.
# Available options: spce, spc, tip3p, tip4p, tip5p, ...
# (For TIP4P-2005, set `water="tip4p"` and manually replace the corresponding .itp file.)
# Note: For rheological MD (rheoMD) simulations, use the SPCE water model ("spce" / "spc216.gro").
# TIP4P-based models are known to exhibit instability and may lead to simulation crashes under shear conditions.
water="tip4p"      
water_file="tip4p.gro"

: '
*************************************************************
No further modifications should be needed beyond this point.

Setting up directories for the run.
*************************************************************
'

if [ -d "./runs" ]; then
    echo "[Chikaterasu] Warning: './runs' directory already exists."
    echo "[Chikaterasu] This suggests a previous run was initiated here."
    echo "[Chikaterasu] While this is usually harmless, it could lead to potential artifacts."
    echo "[Chikaterasu] For example, errors that should have caused the simulation to crash"
    echo "[Chikaterasu] might be bypassed if existing data from a previous run is used."
    echo "[Chikaterasu] This could happen if an entire simulation was copied or moved"
    echo "[Chikaterasu] to this folder for a fresh run."
    #exit 1 # Optional: uncomment to force an exit if folder exists (untested)
else
    echo "[Chikaterasu] No existing './runs' directory found. Starting a fresh run."
fi

if [ -z "$1" ]
then
    echo "[Chikaterasu] Command line arguments are empty. Using manually set debug level $debug_level."
else
    echo "[Chikaterasu] Command line arugments provided. Using first argument as debug level $1."
    debug_level=$1
fi

mkdir -p gromacs
rm -rf gromacs/top
rm -rf gromacs/solvation
rm -rf gromacs/addions
rm -rf gromacs/emin
#rm -rf gromacs/mdp   # unused?
mkdir -p gromacs/top
mkdir -p gromacs/solvation
mkdir -p gromacs/addions
mkdir -p gromacs/emin
#mkdir -p gromacs/mdp
mkdir -p gromacs/coord
mkdir -p runs
mkdir -p runs/nvt
mkdir -p runs/npt
mkdir -p runs/md
mkdir -p custom_analysis

: '
*************************************************************
[pdb2gmx]
Generate topology and position restraints
*************************************************************
'

cd gromacs/coord

# Check if user-defined file exists. Exit if file not found.
if [ ! -e $protein_name.pdb ]; then
    echo "[Chikaterasu] PDB input file does not exist in ./gromacs/coord/ folder. Exiting!"
    exit 2
  else
    echo "[Chikaterasu] Found PDB input file. Will now convert to gromacs format."
fi

# Case a: normal protein simulation, no added molecules
# insert_small_molecules=false

if [ "$insert_small_molecules" = false ] ; then
    if [ "$his_manual" = true ] ; then
        gmx pdb2gmx -f $protein_name.pdb -o ../top/$protein_name.pdb_processed.gro -p ../top/topol.top -water $water -chainsep interactive -ignh -rtpres -merge interactive -his
    fi

    if [ "$his_manual" = false ] ; then
        if [ "$disulfide" = false ] ; then
            if [ "$amber" = true ] ; then
                gmx pdb2gmx -f $protein_name.pdb -o ../top/$protein_name.pdb_processed.gro -p ../top/topol.top -water $water -chainsep interactive -ignh -rtpres -merge interactive
            fi
            if [ "$amber" = false ] ; then
                gmx pdb2gmx -f $protein_name.pdb -o ../top/$protein_name.pdb_processed.gro -p ../top/topol.top -water $water -chainsep interactive -ignh -rtpres -merge interactive -ter
            fi
        fi
        if [ "$disulfide" = true ] ; then
          gmx pdb2gmx -f $protein_name.pdb -o ../top/$protein_name.pdb_processed.gro -p ../top/topol.top -water $water -chainsep interactive -ignh -rtpres -merge interactive -ss
        fi
    fi
fi

# Case b: only small molecule simulation, no protein
# insert_small_molecules=true

if [ "$insert_small_molecules" = true ] ; then
    if [ "$protein_with_small_molecules" = false ]; then
        echo "" > empty.pdb
        gmx insert-molecules -f empty.pdb -ci $protein_name.pdb -o teibunshibox.pdb -nmol $insert_small_molecules_number -box $box_dim
        gmx pdb2gmx -f teibunshibox.pdb -o ../top/$protein_name.pdb_processed.gro -p ../top/topol.top -water $water -chainsep interactive -ignh -rtpres -merge interactive
    fi

    # New: Case c: small molecule and protein. Add molecule to a copy for the initial protein file.
    if [ "$protein_with_small_molecules" = true ]; then
        gmx insert-molecules -f $protein_name.pdb -ci $protein_added_small_molecule_name.pdb -o protein_with_added_molecules.pdb -nmol $insert_small_molecules_number -box $box_dim
        gmx pdb2gmx -f protein_with_added_molecules.pdb -o ../top/$protein_name.pdb_processed.gro -p ../top/topol.top -water $water -chainsep interactive -ignh -rtpres -merge interactive
    fi
fi

mv *.itp ../top/
rm \#*
rm ../top/\#*

cd ../..

if [ "$debug_level" = 1 ] ; then
    echo "[Chikaterasu] Debug level 1 set. Exiting after initial topology generation [pdb2gmx]"
    exit 1
fi


: '
*************************************************************
Solvate the protein 
If defined set the box dimensions manually

Gromacs will print a warning starting with:
[WARNING: Masses and atomic (Van der Waals) radii will be guessed ...]
... but that is fine.
See: https://www.researchgate.net/post/MD-simulation-using-gromacs-using-the-command-gmx-solvate
*************************************************************
'

cd gromacs/solvation

gmx editconf -f ../top/$protein_name.pdb_processed.gro -o ./$protein_name.pdb_newbox.gro -c -d 1.0 -bt $cell_shape

if [ "$box_manual" = true ] ; then
    sed -i '$ d' $protein_name.pdb_newbox.gro 
    echo "$box_dim" >> $protein_name.pdb_newbox.gro
fi

if [ "$box_empty" = true ] ; then
    echo "WATER" > $protein_name.pdb_newbox.gro
    echo " 0" >> $protein_name.pdb_newbox.gro
    echo "$box_dim" >> $protein_name.pdb_newbox.gro
fi

gmx solvate -cp ./$protein_name.pdb_newbox.gro -cs $water_file -o ./$protein_name.pdb_solv.gro -p ../top/topol.top

cd ../..

if [ "$debug_level" = 2 ] ; then
    echo "[Chikaterasu] Debug level 2 set. Exiting after solvation."
    exit 1
fi

: '
*************************************************************
Add ions

On some builds GROMACS will complain that charge is not neutral.
Although the charge is like 0.0000001. So it might be a precision
issue. But of course we want to neutralize it.
The maxwarn flag helps in such a case. Otherwise, recompile
and hope for better results. Just make sure the topology
of the system is neutral at the end.
*************************************************************
'

cd gromacs/addions
cp ../../chika_mdp/ions.mdp ./ions.mdp

gmx grompp -f ./ions.mdp -c ../solvation/$protein_name.pdb_solv.gro -p ../top/topol.top -o ./ions.tpr -maxwarn 1
#exit 1

if [ "$specify_salt_concentration" = true ] ; then
    if [ "$magnesium" = false ] ; then
        printf "SOL" | gmx genion -s ions.tpr -o $protein_name.pdb_solv_ions.gro -p ../top/topol.top -conc $salt_concentration -neutral
    fi
    if [ "$magnesium" = true ] ; then
        printf "SOL" | gmx genion -s ions.tpr -o $protein_name.pdb_solv_ions.gro -p ../top/topol.top -conc $salt_concentration -neutral -pname MG -nname CL
        #exit 1
    fi
fi

if [ "$specify_salt_concentration" = false ] ; then
    if [ "$magnesium" = false ] ; then
        printf "SOL" | gmx genion -s ions.tpr -o $protein_name.pdb_solv_ions.gro -p ../top/topol.top -pname NA -nname CL -np $pos_ions -nn $neg_ions
    fi
    if [ "$magnesium" = true ] ; then
        printf "SOL" | gmx genion -s ions.tpr -o $protein_name.pdb_solv_ions.gro -p ../top/topol.top -pname MG -nname CL -np $pos_ions -nn $neg_ions
    fi
fi

cd ../..

if [ "$debug_level" = 3 ] ; then
    echo "[Chikaterasu] Debug level 3 set. Exiting after adding counterions."
    exit 1
fi

: '
*************************************************************
Add distance restraints
e.g. Zn2+ and other ligands or to fix one domain on a
     receptor etc.

Requires that the file "distance_restraints.itp" is prepared
and in the main directory of this MD simulation

Make sure that chika_mdp files are correctly set up:

NVT
define  = -DPOSRES -DDISRES     ; position restrain the protein
disre   = simple                ; Enable Distance Restraints

NPT
define  = -DPOSRES -DDISRES     ; position restrain the protein
disre   = simple                ; Enable Distance Restraints

MD

define  = -DDISRES              ; only distance restrains left
disre   = simple                ; Enable Distance Restraints

Make sure that numbers in distance restraint file match
the topology numbers even after hydrogens are added and 
CYS is converted to CYM etc.
*************************************************************
'

if [ "$distance_restraints" = true ] ; then
    cd gromacs/top
    cp ../../distance_restraints.itp ./

    # untested addition
    if [[ ! -f ./distance_restraints.itp ]] ; then
      echo "[Chikaterasu-dev] File distance_restraints.itp is not there, aborting."
      exit 1
    fi

    # Add the include statement to the toplogy file [topol.top]
    #; Include distance restraints
    #ifdef DISRES
    #include "./distance_restraints.itp"
    #endif

    sed -i '/Include Position restraint file/a ; Include distance restraints\n#ifdef DISRES\n#include "./distance_restraints.itp"\n#endif\n; Include Position restraint file' topol.top

    cd ../..

    #exit 1

    echo "[Chikaterasu-dev] Added distance restraints to topology file."
    echo "[Chikaterasu-dev] Make sure mdp file md.mdp also contains define -DDISRES . Otherwise the distance restraints are ignored even though we included them in the topology."
    echo "[Chikaterasu-dev] Make sure the topology file is correct."
    echo "[Chikaterasu-dev] e.g., if using Zn2+ ion, confirm the topol.top now to make sure the correct protein-Zn2+ bonds have been added to the [bonds] and before the [pairs] section."

    sed -i "/\[ pairs ]/i $(sed ':a;N;$!ba;s/\n/\\n/g' zn_bond.top)" gromacs/top/topol.top

    read -p "[Chikaterasu-dev] Distance restraints are ON. zn_bond.top has been appended to the topology. Continue?" dummy


fi

: '
*************************************************************
Energy minimization
*************************************************************
'

cd gromacs/emin
cp ../../chika_mdp/minim.mdp ./minim.mdp

gmx grompp -f ./minim.mdp -c ../addions/$protein_name.pdb_solv_ions.gro -p ../top/topol.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em

printf "Potential" | gmx energy -f em.edr -o ./potential.xvg

cd ../..

if [ "$debug_level" = 4 ] ; then
    echo "[Chikaterasu] Debug level 4 set. Exiting after energy minimization."
    exit 1
fi

: '
*************************************************************
Start the MD loop
*************************************************************
'

for i in `seq 1 $nruns`;

do
        # 1. NVT
        # NVT also generates a warning in the newer gromacs versions like this:
        #     Removing center of mass motion in the presence of position restraints
        #     might cause artifacts. When you are using position restraints to
        #     equilibrate a macro-molecule, the artifacts are usually negligible.
        # I never saw any problems here and gromacs also thinks it is fine like this.
        cd runs/nvt
        cp ../../chika_mdp/nvt.mdp ./nvt.mdp

        gmx grompp -f ./nvt.mdp -c ../../gromacs/emin/em.gro -r ../../gromacs/emin/em.gro -p ../../gromacs/top/topol.top -o nvt.tpr -maxwarn 1
        gmx mdrun -deffnm nvt -nb gpu -v

        if [ "$debug_level" = 5 ] ; then
            echo "[Chikaterasu] Debug level 5 set. Exiting after NVT equilibration."
            exit 1
        fi

        printf "Temperature" | gmx energy -f nvt.edr -o ./temperature.xvg
        cd ../..

        # 2. NPT
        cd runs/npt
        cp ../../chika_mdp/npt.mdp ./npt.mdp

        gmx grompp -f ./npt.mdp -c ../nvt/nvt.gro -r ../nvt/nvt.gro -t ../nvt/nvt.cpt -p ../../gromacs/top/topol.top -o npt.tpr -maxwarn 1
        gmx mdrun -deffnm npt -nb gpu -v

        if [ "$debug_level" = 6 ] ; then
            echo "[Chikaterasu] Debug level 6 set. Exiting after NPT equilibration."
            exit 1
        fi

        # Diagnostics: The following commands are disabled by default and can be uncommented to
        # extract additional thermodynamic and box-condition metrics when needed (for example, if
        # variation in temperature, pressure or box dimensions must be verified).
        #printf "Temperature" | gmx energy -f npt.edr -o ./temperature.xvg
        #printf "Density" | gmx energy -f npt.edr -o ./density.xvg
        #printf "Pres-XX" | gmx energy -f npt.edr -o ./pressure_XX.xvg
        #printf "Pres-YY" | gmx energy -f npt.edr -o ./pressure_YY.xvg
        #printf "Pres-ZZ" | gmx energy -f npt.edr -o ./pressure_ZZ.xvg
        #printf "Box-XX" | gmx energy -f npt.edr -o ./box_XX.xvg
        #printf "Box-YY" | gmx energy -f npt.edr -o ./box_YY.xvg
        #printf "Box-ZZ" | gmx energy -f npt.edr -o ./box_ZZ.xvg
        #printf "Volume" | gmx energy -f npt.edr -o ./volume.xvg
        #printf "Pres-XY" | gmx energy -f npt.edr -o ./pressure_XY.xvg
        #printf "Pres-XZ" | gmx energy -f npt.edr -o ./pressure_XZ.xvg
        #printf "Pres-YX" | gmx energy -f npt.edr -o ./pressure_YX.xvg
        #printf "Pres-YZ" | gmx energy -f npt.edr -o ./pressure_YZ.xvg
        #printf "Pres-ZX" | gmx energy -f npt.edr -o ./pressure_ZX.xvg
        #printf "Pres-ZY" | gmx energy -f npt.edr -o ./pressure_ZY.xvg

        cd ../..

        # 3. MD
        cd runs/md

        cp ../../chika_mdp/md.mdp ./md.mdp

        gmx grompp -f ./md.mdp -c ../npt/npt.gro -t ../npt/npt.cpt -p ../../gromacs/top/topol.top -o md.tpr -maxwarn 1

        # Optional: Use the following commands only if you need to manually inspect the simulation box or setup
        # before starting the production run. These steps are not required under normal conditions.
        # (Temporarily disabled, since routine runs do not require manual verification.)
        # gmx editconf -f md.tpr -o system_before_md.pdb
        # gmx editconf -f md.tpr -o md.gro

        # Note (2022-05-09): When running rheological (shear-flow) MD simulations using the `mdrun -deform` option,
        # certain systems may exhibit issues with Thread-MPI optimization. In such cases, using the fallback command
        # below resolves the problem, though it disables multi-thread optimization:
        #   gmx mdrun -deffnm md -ntmpi 1 -v
        # For most modern hardware, the standard command works as expected.
        gmx mdrun -deffnm md -nb gpu -v

        # Diagnostics: The following commands are disabled by default and can be uncommented to
        # extract additional thermodynamic and box-condition metrics when needed (for example, if
        # variation in temperature, pressure or box dimensions must be verified).
        #printf "Temperature" | gmx energy -f md.edr -o ./temperature.xvg
        #printf "Density" | gmx energy -f md.edr -o ./density.xvg
        #printf "Pres-XX" | gmx energy -f md.edr -o ./pressure_XX.xvg
        #printf "Pres-YY" | gmx energy -f md.edr -o ./pressure_YY.xvg
        #printf "Pres-ZZ" | gmx energy -f md.edr -o ./pressure_ZZ.xvg
        #printf "Box-XX" | gmx energy -f md.edr -o ./box_XX.xvg
        #printf "Box-YY" | gmx energy -f md.edr -o ./box_YY.xvg
        #printf "Box-ZZ" | gmx energy -f md.edr -o ./box_ZZ.xvg
        #printf "Volume" | gmx energy -f md.edr -o ./volume.xvg
        #printf "Pres-XY" | gmx energy -f md.edr -o ./pressure_XY.xvg
        #printf "Pres-XZ" | gmx energy -f md.edr -o ./pressure_XZ.xvg
        #printf "Pres-YX" | gmx energy -f md.edr -o ./pressure_YX.xvg
        #printf "Pres-YZ" | gmx energy -f md.edr -o ./pressure_YZ.xvg
        #printf "Pres-ZX" | gmx energy -f md.edr -o ./pressure_ZX.xvg
        #printf "Pres-ZY" | gmx energy -f md.edr -o ./pressure_ZY.xvg

        cd ../..

        # Copy the data
        mkdir -p runs/md_$i
        mkdir -p runs/nvt_$i
        mkdir -p runs/npt_$i
        
        mv runs/md/* runs/md_$i/
        mv runs/nvt/* runs/nvt_$i/
        mv runs/npt/* runs/npt_$i/
        
        mkdir -p runs/md
        mkdir -p runs/nvt
        mkdir -p runs/npt
 
        # End the run
        echo [Chikaterasu] Run $i finished. Yay!
done

exit 1
