# Step 1: Generate the target index by selecting backbone atoms and manually refining the topology.
# 1. Use GROMACS to create the index.ndx file with backbone atoms:
#    Command: gmx make_ndx -f ../md.tpr
# 2. When prompted, select option 4 for backbone atoms.
# 3. Exit the make_ndx tool by pressing 'q', which saves the index file.
# 4. Open index.ndx in a text editor (e.g., vi) and remove all atom groups except for the backbone group.
# 5. Manually edit the backbone group to exclude atoms that should not be considered in the PCA (e.g., long loops, linkers).
# 6. Execute this script from the `results/md_1/pca` directory after running Chikaterasu analysis.

echo "[ Chikaterasu ] PCA analysis will take a while ..."
echo "[ Chikaterasu ] 1. Building covariance matrix ..."

# Step 2: Build the covariance matrix from the molecular dynamics trajectory.
gmx covar -s ../md_target.tpr -f ../md_fit.xtc -o ./eigenval.xvg -av ./average.pdb -n pca_target.ndx -l covar.log -xpm covar.xpm -xpma covar_atomic.xpm -ascii covar.dat
gmx xpm2ps -f covar.xpm
gmx xpm2ps -f covar_atomic.xpm

# Step 3: Analyze the first three principal components.
echo "[ Chikaterasu ] 2. Analyze first three eigenvectors ..."

# Uncomment the following line to perform eigenvector analysis:
#gmx anaeig -s ../md_target.tpr -f ../md_fit.xtc -first 1 -last 3 -extr ./extreme_pcavec.pdb -nframes 30 -entropy -n pca_target.ndx

# Step 4: Generate 2D projections of the first three principal components.
echo "[ Chikaterasu ] 3. Plot first 3 PCA vectors as 2D maps ..."
# Uncomment the following line to generate 2D plots:
#gmx anaeig -s ../md_target.tpr -f ../md_fit.xtc -2d ./2dproj_12.xvg -first 1 -last 2 -n pca_target.ndx
