#!/bin/sh

: '
*************************************************************
Chikaterasu         version dev
gmx                 2022.3

Last changes        see github
Issues              see github

Erik Walinda
Kyoto University
Graduate School of Medicine

*************************************************************
'

: '
*************************************************************
Manually setup parameters for this run
Most parameters are defined by mdp files
The starting PDB file should be placed into gromacs/coord

Debug level 
0   full production MD run [NPT]; no debug
1   topology generation [pdb2gmx]
2   solvation
3   addition of counterions
    addition of distance restraints; no debug level implemented yet
4   energy minimization
5   NVT done
6   NPT done
*************************************************************
'

# == What to simulate ? ==

protein_name="1UBQ"         # For small molecule simulations, this will be ATP.pdb etc. even though the variable name says protein.
nruns=1                     # 1 for testing; 10 for production

# == Debug level ==

debug_level=0               # Manually set debug level. Or give as argument, e.g.: ./chikaterasu 0

# == Histidine and Zn2+ stuff ==

his_manual=false            # manually specify histidine protonation state
distance_restraints=false   # for Zn2+
disulfide=false             # do we want to account for disulfide bridges when making the topology?

# === Ions ===

specify_salt_concentration=true  # Specify means in molar; not specify means to count the number of ions.
salt_concentration=0.050

pos_ions=0                  # only if manually specifying ions
neg_ions=2                  # only if manually specifying ions

magnesium=false             # Set Mg as the positive ion; if false it is Na.

# === Box shape and size ===
# Ideally, we want 3 options here:
# a) simulate only a protein:          insert_small_molecules is False
# b) simulate only small molecules     insert_small_molecules is True
# c) both protein with small molecules insert_small_molecules is True
#                                      protein_with_small_molecules is True

# Must set box_dim when working with inserting molecules
insert_small_molecules=false       # If we want to simulate small molecules like ATP or sialic acid
insert_small_molecules_number=2

protein_with_small_molecules=false
protein_added_small_molecule_name="ATP"  # File needs to be in gromacs/coord folder; i.e., together with the protein.

box_manual=false            # Specify box-size manually. 
box_empty=false             # Water-only simulation

box_dim="    7.845   6.497   6.363 "  # if necessary to specify the size. e.g. ATP, rheoMD, hydrodynamics, check PBC artifacts
cell_shape="triclinic"      # -d: triclinic, cubic, dodecahedron, octahedron

water="tip4p"               # [spce], spc, tip3p, tip4p, tip5p, ... (for tip4p2005, use tip4p and replace itp file)
water_file="tip4p.gro"      # [spc216.gro], tip4p.gro, ..., or explicitly with directory (amber03ws.ff/tip4p2005.gro)


: '
*************************************************************
If well programmed no change should be necessary from here

Setup directories for the run
*************************************************************
'

if [ -d "./runs" ]; then
    echo "[Chikaterasu] A folder named ./runs already exists. This suggests that Chikaterasu has been launched in this folder before. While this is generally no problem at all, it can lead to potential artifacts, so be warned. For example, a run that was supposed to crash due to gromacs errors in early stages can continue smoothly all the way without spotting an error just because the runs folder already contains the data for continuation. This happens for example, if an entire simulation was copied moved to start a fresh run."
    #exit 1 # untested
else
    echo "[Chikaterasu] Starting a fresh run."
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
            gmx pdb2gmx -f $protein_name.pdb -o ../top/$protein_name.pdb_processed.gro -p ../top/topol.top -water $water -chainsep interactive -ignh -rtpres -merge interactive
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

define  = -DDISRES              ; only position restrains left
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

    echo "[Chikaterasu-dev] Added distance restraints to topology file. If using Zn2+ ion, manually edit the topol.top now to add protein-Zn2+ bonds in the [bonds] section."
    read -p "[Chikaterasu-dev] Distance restraints are ON. Continue?" dummy

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

        printf "Temperature" | gmx energy -f npt.edr -o ./temperature.xvg
        printf "Density" | gmx energy -f npt.edr -o ./density.xvg
        #printf "Pres-XX" | gmx energy -f npt.edr -o ./pressure_XX.xvg
        #printf "Pres-YY" | gmx energy -f npt.edr -o ./pressure_YY.xvg
        #printf "Pres-ZZ" | gmx energy -f npt.edr -o ./pressure_ZZ.xvg
        #printf "Box-XX" | gmx energy -f npt.edr -o ./box_XX.xvg
        #printf "Box-YY" | gmx energy -f npt.edr -o ./box_YY.xvg
        #printf "Box-ZZ" | gmx energy -f npt.edr -o ./box_ZZ.xvg
        printf "Volume" | gmx energy -f npt.edr -o ./volume.xvg
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

        gmx editconf -f md.tpr -o system_before_md.pdb
        gmx editconf -f md.tpr -o md.gro

        gmx grompp -f ./md.mdp -c ../npt/npt.gro -t ../npt/npt.cpt -p ../../gromacs/top/topol.top -o md.tpr -maxwarn 1
        
        # Note (added 20220509): when using mdrun deform option for rheo. MD, a bug seems to exist somehow.
        # Thread-MPI optimization somehow does not work; however, if the next line is changed to ... it works! (no optimization though?)
        # gmx mdrun -deffnm md -ntmpi 1 -v
        gmx mdrun -deffnm md -nb gpu -v

        printf "Temperature" | gmx energy -f md.edr -o ./temperature.xvg
        printf "Density" | gmx energy -f md.edr -o ./density.xvg
        #printf "Pres-XX" | gmx energy -f md.edr -o ./pressure_XX.xvg
        #printf "Pres-YY" | gmx energy -f md.edr -o ./pressure_YY.xvg
        #printf "Pres-ZZ" | gmx energy -f md.edr -o ./pressure_ZZ.xvg
        #printf "Box-XX" | gmx energy -f md.edr -o ./box_XX.xvg
        #printf "Box-YY" | gmx energy -f md.edr -o ./box_YY.xvg
        #printf "Box-ZZ" | gmx energy -f md.edr -o ./box_ZZ.xvg
        printf "Volume" | gmx energy -f md.edr -o ./volume.xvg
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
