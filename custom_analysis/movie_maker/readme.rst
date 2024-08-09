Movie Maker
-----------

1. Filter down
""""""""""""""

Making movies takes time, especially in Chimera.
Best, first select only the trajectory where a movie is really necessary (i.e., do other analysis first).


1. Get processed trajectory
"""""""""""""""""""""""""""

Set dt = 1000, unless we need really detailed movies.

Execute::

  ./ana_chikaterasu.sh

All options can be off, we just want the centered xtc file with PBC removed.

Copy these two files to a folder which we want to see::

* md_fit.xtc
* md_target.tpr

Convert tpr to gro::

  .\gmx.exe editconf -f md_target.tpr -o .\md_target.gro

2. Load the data
""""""""""""""""

Load both gro and xtc in Chimerax or in VMD. Syntax of course differs - if ChimeraX, select and color as wished::
  
  color /?:1-76 blue
  color /?:77-152 green
  color /?:191-248 grey  
  color @ZN grey

Set snail to 1.
Press "movie record" or write (very slow to supersample)::

  movie record supersample 2

It will play the trajectory and while that is going on, capture the images and then after done encode the movie.

If VMD, make representations for each part of the molecule we are interested in (Selected Atoms in representations)::

  resid 1 to 76
  resid 77 to 152
  resid 191 to 248
  ... and color each by ColorID
  ... and maybe add (transparant) surface representations
  ... and maybe trajectory smoothing to 1

Using nglview
-------------

1. Preparations
"""""""""""""""

In anaconda etc. create an environment and install the requirements::

* nglview
* mdanalysis
* jupyter

