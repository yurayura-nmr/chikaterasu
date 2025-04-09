# This script resumes a GROMACS molecular dynamics (MD) simulation, 
# typically used after an interruption (e.g., power outage).
# Execute this script within the [chikaterasu]/runs/md directory.

# Generate a continued TPR file without extending the simulation time.
gmx convert-tpr -s md.tpr -extend 0 -o md_continued.tpr
# To extend the simulation, specify the desired duration in picoseconds.
# For example, to extend by 100 nanoseconds (100,000 ps):
# gmx convert-tpr -s md.tpr -extend 100000 -o md_continued.tpr

# Resume the MD run using checkpoint data, enabling GPU acceleration where applicable.
gmx mdrun -s md_continued.tpr -cpi md.cpt -deffnm md -nb gpu

# In some cases (sudden power outages), the md.cpt conflict with md.xtc, in that case, adjust the above to:
#gmx mdrun -s md_continued.tpr -cpi md_prev.cpt -deffnm md -nb gpu
