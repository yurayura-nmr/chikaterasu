Movie Maker
-----------

1. Get processed trajectory
"""""""""""""""""""""""""""

Execute::

  ./ana_chikaterasu.sh

All options can be off, we just want the centered xtc file with PBC removed.

Copy these two files to a folder which we want to see::

* md_fit.xtc
* md_target.tpr

Convert tpr to gro::

  .\gmx.exe editconf -f md_target.tpr -o .\md_target.gro

Load both gro and xtc in Chimerax. Select and color as wished::
  
  color /?:1-76 blue
  color /?:77-152 green
  color /?:191-248 grey  
  color @ZN grey

Set snail to 1.
Press "movie record" or write (very slow to supersample)::

  movie record supersample 2

It will play the trajectory and while that is going on, capture the images and then after done encode the movie.

Using nglview
-------------

1. Preparations
"""""""""""""""

In anaconda etc. create an environment and install the requirements::

* nglview
* mdanalysis
* jupyter

