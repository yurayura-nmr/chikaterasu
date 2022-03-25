"""
Erik Walinda
2022/3/11

Adjusted initram.sh for non-ancient gromacs versions since initram writes its own mdp files.

Backward environment currently only set up in Mayuta (misleading path):
/home/arurun/data/md/20220324_backward/for_sorada/extended_cyclic_18A/chikaterasu-main/chikaterasu_m/backward

Required input files:

* ~~~_to_backcalc.pdb   # the CG structure to be backcalculated
* ~~~.tpr               # additional parameters [seems to be the same for cyclic/non-cyclic]
* atomistic pdb (to be converted to topology by chikaterasu)

MDP scripts might still be optimized for GROMACS 2021. They seem old. I.e., I get these warnings:

NOTE 1 [file unknown]:
  You are using constraints on all bonds, whereas the forcefield has been
  parametrized only with constraints involving hydrogen atoms. We suggest
  using constraints = h-bonds instead, this will also improve performance.

NOTE 2 [file 4-mdpr-0.0005.mdp]:
  Removing center of mass motion in the presence of position restraints
  might cause artifacts. When you are using position restraints to
  equilibrate a macro-molecule, the artifacts are usually negligible.

NOTE 3 [file 4-mdpr-0.0005.mdp]:
  You are using a plain Coulomb cut-off, which might produce artifacts.
  You might want to consider using PME electrostatics.

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

# 5. Try to run initram:

# needs py27
# conda create --name py27 python=2.7

conda activate py27
./aru_initram.sh -f cyclic_extended_frame1_to_backcalc-CG.gro -o aa_charmm.gro -to charmm36 -p ./topol.top

# And visualize the output as usual (aa_charmm.gro).
