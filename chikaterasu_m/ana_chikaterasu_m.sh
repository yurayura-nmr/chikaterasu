# very preliminary

# ---- Visualization
# to see motion. Center on protein and prevent jumping (for k48)
gmx editconf -f ../dynamic.tpr -o target.pdb
printf "non-Water\nq\n" | gmx make_ndx -f target.pdb -o target.ndx
printf "Protein\nProtein\n" | gmx convert-tpr -s ../dynamic.tpr  -n target.ndx -o md_target.tpr
#printf "Protein\nq\n" | gmx convert-tpr -s ../dynamic.tpr  -n target.ndx -o md_target.tpr
printf "Protein\nProtein\n" | gmx trjconv -s ../dynamic.tpr -f ../dynamic.xtc -center -ur compact -pbc nojump -o md_target_centered_no_PBC.xtc
printf "Protein\nProtein\n" |gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit rot+trans -o md_fit.xtc
printf "Protein\nProtein\n" |gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit rot+trans -o test.pdb -conect -dt 2500
# ----

# on the fly check distance

# make_ndx as usual
gmx distance -f ../dynamic.xtc -s ../dynamic.tpr -oall distance.xvg -n index.ndx -rmpbc -tu ns
# print max. distance
awk 'NR==2{max = $2 + 0; next} {if ($2 > max) max = $2;} END {print max}' distance.xvg 

#or just
tail -n 30 distance.xvg
