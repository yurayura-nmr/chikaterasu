#!/bin/sh

: '
*************************************************************
Chikaterasu_m       version dev
gmx                 tested for versions ....
                    martini 2.2
                    
Last change         see github

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

box_manual=true
box_dim="8.000   8.000   8.000"
protein_name="k48-CG"

# top/itp must be manually prepared and ready in folder
top="k48" # name

# 1. editconf

cd gromacs
cd coord
gmx editconf -f $protein_name.pdb -d 1.0 -bt triclinic -o $protein_name.gro

if [ "$box_manual" = true ] ; then
    sed -i '$ d' $protein_name.gro
    echo "$box_dim" >> $protein_name.gro
fi

# 2. use top

cd ../top/
cp ../../chika_mdp/martini.itp .
cp ../../*.top .
cp ../../*.itp .

gmx grompp -p $top.top -f ../../chika_mdp/minimization.mdp -c ../coord/$protein_name.gro -o ../emin/minimization-vac.tpr

# Vacuum minimization
cd ../emin/
gmx mdrun -deffnm minimization-vac -v

# 3. solvation

cd ../solvation

cp ../../chika_mdp/water.gro .
gmx solvate -cp ../emin/minimization-vac.tpr -cs water.gro -radius 0.21 -o system_solvated.gro

# Now we need to count how many waters we added
cp ../top/$top.top ../top/system_solvated.top
echo -n "W\t\t" >> ../top/system_solvated.top
grep -c W system_solvated.gro >> ../top/system_solvated.top


# Solution minimization
cd ../emin
gmx grompp -p ../top/system_solvated.top -c ../solvation/system_solvated.gro -f ../../chika_mdp/minimization.mdp -o minimization-sol.tpr
gmx mdrun -deffnm minimization-sol -v

# 4. Equilibration

cd ../../runs/
mkdir -p npt
cd npt
cp ../../gromacs/top/*.top .
cp ../../gromacs/top/*.itp .
cp ../../gromacs/emin/minimization-sol.gro .

gmx grompp -p system_solvated.top -c minimization-sol -f ../../chika_mdp/equilibration.mdp -o equilibration.tpr
gmx mdrun -deffnm equilibration -v 

# 5. Production

cd ..
mkdir -p md
cd md

cp ../npt/*.top .
cp ../npt/*.itp .
cp ../npt/equilibration.gro .
cp ../npt/equilibration.tpr .

gmx grompp -p system_solvated.top -f ../../chika_mdp/dynamic.mdp -o dynamic.tpr -c equilibration.tpr -maxwarn 1
#gmx grompp -p system_solvated.top -c equlibration.tpr -f ../../chika_mdp/dynamic.mdp -o dynamic.tpr -maxwarn 1
gmx mdrun -deffnm dynamic -v 


