Erik Walinda
2022/3/11

Adjusted initram.sh for non-ancient gromacs versions since initram writes its own mdp files.

1. Create gro file of CG model using trjconv (take care of PBC).
2. In atomistic PDB convert GLQ back to GLY and LYQ back to LYS.
3. Run "chikaterasu 1" to obtain atomistic topology (CHARMM27 forcefield).
4. Run "chikaterasu 5" to obtain the file posre.itp (used by BACKWARD).
5. Run "chikaterasu 1" again to obtain topology file without any new stuff (introduced by chikaterasu 5).
6. Try to run initram:

./initram.sh -f CG_posre.gro -o aa_charmm.gro -to charmm36 -p gromacs/top/topol.top

And visualize the output as usual (aa_charmm.gro).
