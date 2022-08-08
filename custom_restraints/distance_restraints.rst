If using distance restraints, need to modify the mdp files in the chika_mdp folder
----------------------------------------------------------------------------------

This is not automatically done.

To nvt.mdp
""""""""""

define  = -DPOSRES -DDISRES     ; position restrain the protein
disre   = simple                ; Enable Distance Restraints

To npt.mdp
""""""""""

define  = -DPOSRES -DDISRES     ; position restrain the protein
disre   = simple                ; Enable Distance Restraints

To md.mdp
"""""""""

define  = -DDISRES              ; only position restrains left
disre   = simple                ; Enable Distance Restraints
