Chikaterasu
===========

.. image:: logo.png
   :alt: Chikaterasu logo
   :align: right

Bash script to automate setup of multiple identical MD simulations.
MD parameters are for the amber-type forcefields such as amber99sb-ildn and amber03ws.

Usage
-----

1. Put PDB file of your system to be simulated into ./gromacs/coord/          [e.g. 1UBQ.pdb]
2. Edit the top section of chikaterasu.sh to specify its filename             [e.g. 1UBQ]
3. Set MD parameters for NVT, NPT, and production MD in chika_mdp/~.mdp files [time, temperature, amounts of frames saved, etc.]
4. sh chikaterasu.sh 1 -> test if pdb format to gromacs conversion works. Make sure no errors or warnings are given by GROMACS.
5. sh chikaterasu.sh 2 -> test if protein could be solvated in the box. Make sure no errors before continuing.
6. sh chikaterasu.sh 3 ->
7. sh chikaterasu.sh 4 ->
8. sh chikaterasu.sh 5 -> Test if the first equilibration step (100 ps NVT) is working without any issues.
9. sh chikaterasu.sh 6 -> Test if all equilibration steps including 100 ps NPT are working without any issues.
10. sh chikaterasu.sh 0 -> Start production MD

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
