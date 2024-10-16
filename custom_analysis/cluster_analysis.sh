# After running the chikaterasu analysis, navigate to the results/md_1/traj/ directory.
# Run GROMACS cluster analysis on the system trajectory using the following:
# 
# - The trajectory file (-f) md_target.xtc and structure file (-s) md_target.tpr
# - GROMOS clustering method (-method gromos) with a 1 nm RMSD cutoff (-cutoff 1)
# - Output the cluster representative structures to a pdb file (-cl). You can view this with pymol.
# - Use a frame every 100 ps (-dt 100)
#
# Optionally, provide an index file (-n) to analyze a specific subgroup (e.g., C-alpha atoms).
gmx cluster -f ../md_target.xtc -s ../md_target.tpr -method gromos -cutoff 1 -cl -dt 100 
