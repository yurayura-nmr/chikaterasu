# very preliminary

# ---- Visualization
# to see motion. Center on protein and prevent jumping (for k48)
gmx editconf -f ../dynamic.tpr -o target.pdb
printf "non-Water\nq\n" | gmx make_ndx -f target.pdb -o target.ndx
gmx convert-tpr -s ../dynamic.tpr  -n target.ndx -o md_target.tpr
gmx trjconv -s ../dynamic.tpr -f ../dynamic.xtc -center -ur compact -pbc nojump -o md_target_centered_no_PBC.xtc
gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit rot+trans -o md_fit.xtc
gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit rot+trans -o test.pdb -conect
# ----

# on the fly check distance

# make_ndx as usual
gmx distance -f ../dynamic.xtc -s ../dynamic.tpr -oall distance.xvg -n index.ndx -rmpbc -tu us
