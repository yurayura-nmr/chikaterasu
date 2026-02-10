From pymol generated file, remove TER.

conda install conda-forge::ambertools

May require python 3.11 at the moment (need to make a dedicated environment)

pdb4amber -i your_protein.pdb -o processed_your_protein.pdb --no-hyd --add-missing-atoms

Rename HETATM to ATOM

Rename C-terminal NMET cap residues exactly the way they are in the demo data

