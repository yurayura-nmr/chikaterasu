#!/bin/sh

: '
*************************************************************
Chikaterasu_m       version dev
gmx                 tested for gmx version 2019.2
martini             tested for ff  version 2.2
                    
Last change         see github

Erik Walinda
Kyoto University
Graduate School of Medicine

*************************************************************
'

: '
*************************************************************
Manually setup parameters for this run
Most parameters are defined by mdp files in chika_mdp/
The starting PDB file should be placed into gromacs/coord

Debug level 

0   full production MD run (no debug)
1   topology generation (editconf)
2   solvation (solvate)
3   addition of counterions
4   equilibration before production
*************************************************************
'

# == Debug level ==

debug_level=1               # Manually set debug level. Or give as argument, e.g.: ./chikaterasu 0

# === Box shape and size ===
box_manual=true
box_dim="8.000   8.000   8.000"

# top/itp must be manually prepared and ready in folder
top="k48"               # name
protein_name="k48-CG"   # -CG coarse grained (martinized) pdb

: '
*************************************************************
If well programmed no change should be necessary from here
Setup directories for the run
*************************************************************
'

if [ -z "$1" ]
then
    read -p "[Chikaterasu] Command line arguments are empty. Using manually set debug level $debug_level." dummy
else
    read -p "[Chikaterasu] Command line arugments provided. Using first argument as debug level $1." dummy
    debug_level=$1
fi

mkdir -p gromacs
rm -rf gromacs/top
rm -rf gromacs/solvation
rm -rf gromacs/addions
rm -rf gromacs/emin
mkdir -p gromacs/top
mkdir -p gromacs/solvation
mkdir -p gromacs/addions
mkdir -p gromacs/emin
mkdir -p gromacs/coord
mkdir -p runs
mkdir -p runs/nvt
mkdir -p runs/npt
mkdir -p runs/md
mkdir -p custom_analysis

: '
*************************************************************
Assuming that the topology has already been set up using martinize

Example for K48 diUb
python2.7 martinize.py -f 1aar_modified.pdb -o 1aar_modified.top -x 1aar_modified-CG.pdb -dssp dssp -p backbone -merge A,B -elastic -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0

To let bonds decay
python2.7 martinize.py -f 1aar_modified.pdb -o 1aar_modified.top -x 1aar_modified-CG.pdb -dssp dssp -p backbone -merge A,B -elastic -ef 500 -el 0.5 -eu 0.9 -ea 1 -ep 1

To identify cross-bonds between subunits (we want to avoid these)
awk -v x=163 '($1 <= x) && ($2 >x)' Protein_A+Protein_B.itp > to_replace.txt
Open in vsCode side by side and comment out respective lines in Pro~Pro.itp

Additional bonds that might be too stiff (C-terminus!)
  154   162      6   0.87282 RUBBER_FC*0.688791   ; R72-G75
  154   163      6   0.76329 RUBBER_FC*0.768516   ; R72-G76
  157   163      6   0.47946 RUBBER_FC*1.020747   ; L73-G76
  
Additional bonds that might be too stiff (constraining C-terminus!)
   77   154      6   0.85887 RUBBER_FC*0.698464   ; keeps R72 artificially in place
   ...  154-157 etc.                              ; everything that links R72, L73 to the main core


Example for monoUb
python2.7 martinize.py -f 1UBQ.pdb -o single-ubq.top -x 1UBQ-CG.pdb -dssp dssp -p backbone

*************************************************************
'

: '
*************************************************************
[editconf]

Generate box around protein and perform 1 vacuum minimization
*************************************************************
'

cd gromacs
cd coord

# Create box with editconf
gmx editconf -f $protein_name.pdb -d 1.0 -bt triclinic -o $protein_name.gro

if [ "$box_manual" = true ] ; then
    sed -i '$ d' $protein_name.gro
    echo "$box_dim" >> $protein_name.gro
fi

# Use topology to prepare vacuum minimization

cd ../top/
cp ../../chika_mdp/martini.itp .
cp ../../*.top .
cp ../../*.itp .

# Vacuum minimization
gmx grompp -p $top.top -f ../../chika_mdp/minimization.mdp -c ../coord/$protein_name.gro -o ../emin/minimization-vac.tpr
cd ../emin/
gmx mdrun -deffnm minimization-vac -v

if [ "$debug_level" = 1 ] ; then
    echo "[Chikaterasu] Debug level 1 set. Exiting after initial topology generation [pdb2gmx]"
    exit 1
fi

: '
*************************************************************
Solvate the protein

Then perform minization of solvated system
*************************************************************
'

cd ../solvation
cp ../../chika_mdp/water.gro .

# Solvation
gmx solvate -cp ../emin/minimization-vac.tpr -cs water.gro -radius 0.21 -o system_solvated.gro

# Now we need to count how many waters we added
cp ../top/$top.top ../top/system_solvated.top
echo -n "\nW\t\t" >> ../top/system_solvated.top
grep -c W system_solvated.gro >> ../top/system_solvated.top

# Solution minimization
cd ../emin
gmx grompp -p ../top/system_solvated.top -c ../solvation/system_solvated.gro -f ../../chika_mdp/minimization.mdp -o minimization-sol.tpr
gmx mdrun -deffnm minimization-sol -v

if [ "$debug_level" = 2 ] ; then
    echo "[Chikaterasu] Debug level 2 set. Exiting after solvation."
    exit 1
fi

: '
*************************************************************
Add ions

...
*************************************************************
'

# IONS NOT IMPLEMENTED YET

# But seems to work mostly the same as normal gromacs

# Pseudocode
# get ions.mdp from chika_mdp
# get martini_ions.itp and somehow include it in the topology
# gmx grompp -f ./ions.mdp -c ../solvation/system_solvated.gro -p ../top/system_solvated.top -o ./ions.tpr    # -maxwarn 1
# gmx genion -s ./ions.tpr -o solv_ions.gro -p ../top/system_solvated.top -conc 0 -neutral
# and then continue with the solv ions.
# check if topology updated with ions


if [ "$debug_level" = 3 ] ; then
    echo "[Chikaterasu] Debug level 3 set. Exiting after adding counterions."
    exit 1
fi

: '
*************************************************************
Equilibration
*************************************************************
'

cd ../../runs/

mkdir -p npt
cd npt

cp ../../gromacs/top/*.top .
cp ../../gromacs/top/*.itp .
cp ../../gromacs/emin/minimization-sol.gro .

# Equilibration before production MD
gmx grompp -p system_solvated.top -c minimization-sol -f ../../chika_mdp/equilibration.mdp -o equilibration.tpr
gmx mdrun -deffnm equilibration -v 

if [ "$debug_level" = 4 ] ; then
    echo "[Chikaterasu] Debug level 4 set. Exiting after equilibration."
    exit 1
fi

: '
*************************************************************
Production MD (only 1 run)
*************************************************************
'

cd ..
mkdir -p md
cd md

cp ../npt/*.top .
cp ../npt/*.itp .
cp ../npt/equilibration.gro .
cp ../npt/equilibration.tpr .

# Run. Gromacs recognizes what to do with the GPU without specifically saying "-nb gpu"
gmx grompp -p system_solvated.top -f ../../chika_mdp/dynamic.mdp -o dynamic.tpr -c equilibration.tpr -maxwarn 1
gmx mdrun -deffnm dynamic -v 

exit 1
