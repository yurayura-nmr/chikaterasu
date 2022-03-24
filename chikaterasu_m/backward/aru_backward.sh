"""
Erik Walinda
2022/3/11

Adjusted initram.sh for non-ancient gromacs versions since initram writes its own mdp files.
"""


# 1. Create gro file of CG model using trjconv (take care of PBC).

gmx trjconv -f cyclic_extended_frame1_to_backcalc.pdb -s dynamic.tpr -o cyclic_extended_frame1_to_backcalc-CG.gro 

# 2. In atomistic PDB convert GLQ back to GLY and LYQ back to LYS.

vi ..
:%s/LYQ/LYS/g
:%s/GLQ/GLY/g

# 3. Run "chikaterasu 1" to obtain atomistic topology (CHARMM27 forcefield).
# go back to chikaterasu folder (non-martini version)

cd ..
cd ..
./chikaterasu.sh 1
cp gromacs/top/topol.top chikaterasu_m/backward/

# 4. Run "chikaterasu 5" to obtain the file posre.itp (used by BACKWARD).

./chikaterasu.sh 5
cp gromacs/top/posre.itp chikaterasu_m/backward/

# 5. Run "chikaterasu 1" again to obtain topology file without any new stuff (introduced by chikaterasu 5).
cp topolo.top chikaterasu_m/backward/

#/ 6. Try to run initram:

# needs py27
# conda create --name py27 python=2.7

conda activate py27
./aru_initram.sh -f cyclic_extended_frame1_to_backcalc-CG.gro -o aa_charmm.gro -to charmm36 -p ./topol.top

# And visualize the output as usual (aa_charmm.gro).
