# This continues a gromacs run e.g., after a power outage etc.
# Run within [chikaterasu]/runs/md folder

gmx convert-tpr -s md.tpr -extend 0 -o md_continued.tpr

# Time unit is picoseconds. So to extend by 100ns: -extend 100000
# gmx convert-tpr -s md.tpr -extend 100000 -o md_continued.tpr

gmx mdrun -s md_continued.tpr -cpi md.cpt -deffnm md -nb gp
