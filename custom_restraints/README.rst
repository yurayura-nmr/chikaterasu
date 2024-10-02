Using distance restraints in GROMACS
------------------------------------

If using distance restraints, need to modify the mdp files in the chika_mdp folder as follows (this is not automatically done).

In case of Zn2+, just the Zn-protein bonds might be enough to keep the Zn ion there.
If not, additional distance restraints can be added as follows.

To nvt.mdp::

  define  = -DPOSRES -DDISRES     ; position restrain the protein (default) and additionally add distance restraints as specified by custom distance_restraints.itp file.
  disre   = simple                ; Enable Distance Restraints

To npt.mdp::

  define  = -DPOSRES -DDISRES     ; position restrain the protein
  disre   = simple                ; Enable Distance Restraints

To md.mdp::

  define  = -DDISRES              ; only distance restrains left (no need to position restraint in production MD)
  disre   = simple                ; Enable Distance Restraints
