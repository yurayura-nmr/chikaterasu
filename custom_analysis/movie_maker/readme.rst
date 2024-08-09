Movie Maker
-----------

1. Preparations
"""""""""""""""

In anaconda etc. create an environment and install the requirements::

* nglview
* mdanalysis
* jupyter

2. Get processed trajectory
"""""""""""""""""""""""""""

Execute::

  ./ana_chikaterasu.sh

All options can be off, we just want the centered xtc file with PBC removed.

Copy these two files to a folder which we want to see::

* md_target_centered_no_PBC.xtc
* md_target.tpr
