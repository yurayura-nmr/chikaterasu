From pymol generated file, remove TER.

conda install conda-forge::ambertools

May require python 3.11 at the moment (need to make a dedicated environment)

pdb4amber -i your_protein.pdb -o processed_your_protein.pdb

May need to rename HETATM to ATOM

