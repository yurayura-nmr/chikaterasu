Chikaterasu
===========

.. image:: logo.png
   :alt: Chikaterasu logo
   :align: right

Bash script to automate setup of multiple identical MD simulations.
MD parameters are for the amber-type forcefields such as amber99sb-ildn and amber03ws.

Requirements
------------

GROMACS is used. Any recent version should work. 
To install GROMACS, download it, extract it, and build from source with GPU support.
Check gmx --version to see if everything works (GPU/cuda properly recognized? etc.).

Usage
-----

In linux::

  git clone https://github.com/yurayura-nmr/chikaterasu.git

1. Put PDB file of your system to be simulated into ./gromacs/coord/          [e.g. 1UBQ.pdb]
2. Edit the top section of chikaterasu.sh to specify its filename             [e.g. 1UBQ]
3. Set MD parameters for NVT, NPT, and production MD in chika_mdp/~.mdp files [time, temperature, amounts of frames saved, etc.]
4. For a new system, run chikaterasu at various levels to confirm no errors in preparation/execution of MD simulation. If the system is already known to behave well, can go straight to a full run (level 0) and skip the individual levels (1-6).

Run chikaterasu at level 1::

  sh chikaterasu.sh 1 
  # Tests if pdb format to gromacs conversion works. Make sure no errors or warnings are given by GROMACS.

Run chikaterasu at level 2::
 
  sh chikaterasu.sh 2
  # Tests if protein could be solvated in the box. Make sure no errors before continuing.

Run chikaterasu at level 3::
 
  sh chikaterasu.sh 3
  # Tests if ...

Run chikaterasu at level 4::

  sh chikaterasu.sh 4
  # Tests if ...

Run chikaterasu at level 5::

  sh chikaterasu.sh 5
  # Tests if the first equilibration step (100 ps of NVT position-restrained simulation) is working without any issues.

Run chikaterasu at level 6::

  sh chikaterasu.sh 6
  # Tests if all equilibration steps including the final 100-ps NPT position-restrained simulation are working without any issues.

Finally, run chikaterasu at production level (0)::
 
  sh chikaterasu.sh 0
  # Starts production MD

Special bonds
-------------

Zinc
""""

In the case of Zn2+ ions (or similar) in the input PDB file::

  ./chikaterasu.sh 0   # this should run normally

In the chikaterasu base folder, provide an additional file containing the distance restraints for the Zn2+ ion (e.g., to coordinating His or Cys residues). As an example distance_restraints.itp file::

  #ifdef DISRES
  [ distance_restraints ]
  ; ai aj type index type’ low up1 up2 fac
    2564 3356   1   1     1     0.0 0.23 0.3 1
    2595 3356   1   2     1     0.0 0.23 0.3 1
    2766 3356   1   3     1     0.0 0.23 0.3 1
    2809 3356   1   4     1     0.0 0.23 0.3 1
  #endif


Modelling missing loops in the structure
----------------------------------------

1. Go to: https://modbase.compbio.ucsf.edu/modloop/
2. Using registered email address and license key - MODELIRANJE
3. Upload coordinate file - file.pdb

Enter loop segments (residue:chain_ID:residue:chain_ID)::

  70:A:71:A:

For that the uploaded pdb file needs to be tuned so that (in this example) ALA71 already exists. i.e., add dummy atoms like this in a text editor (positions should not matter and probably can even be 0 0 0)::

  ATOM    556  N   ALA A  71      32.763  35.831  23.090  1.00 12.71           N
  ATOM    557  CA  ALA A  71      34.145  35.472  23.481  1.00 16.06           C
  ATOM    558  C   ALA A  71      34.239  35.353  24.979  1.00 18.09           C
  ATOM    559  O   ALA A  71      33.707  36.197  25.728  1.00 19.26           O


Change log
----------

2021-10-24
""""""""""

Added just another folder for user-specific (non-automatable specific) analysis.
(not overwritten by the cleanup function)

Such as specific PCA of only atoms 1-70 of Ub2.
Or just 1 basepair of a DNA.
                    
Before that: (February)
-----------------------

Added Mg ion functionality  [tested a bit, but may still have bugs]

Added insert molecules      [tested a bit, but may still have bugs]


To do
-----

chikaterasu.sh
""""""""""""""

* Issue warning if low on disk space before starting a new run.
* ss untested and only implemented for His=false yet
* re-add dssp function: 
* gmx xpm2ps -f ss.xpm -di dssp.m2p

ana_chikaterasu.sh
""""""""""""""""""
