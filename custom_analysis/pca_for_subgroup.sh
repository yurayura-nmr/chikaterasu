# First, make the target by getting all backbone atoms and then excluding unwanted parts of the topology manually
# 1. Run GROMACS to get index.ndx file:
# gmx make_ndx -f ../md.tpr
# 2. Select 4 for backbone
# 3. Select q to quit and write the index file
# 4. In vi, open index.ndx and remove all groups except for backbone
# 5. In vi, remove those atoms from the backbone that should be excluded from the PCA (long loops, linkers etc.)
# 6. Run this script in the chikaterasu subdirectory results/md_1/pca after running the chikaterasu analysis script.

echo "[ Chikaterasu ] PCA analysis will take a while ..."
echo "[ Chikaterasu ] 1. Building covariance matrix ..."

gmx covar -s ../md_target.tpr -f ../md_fit.xtc -o ./eigenval.xvg -av ./average.pdb -n pca_target.ndx -l covar.log -xpm covar.xpm -xpma covar_atomic.xpm -ascii covar.dat
gmx xpm2ps -f covar.xpm
gmx xpm2ps -f covar_atomic.xpm

echo "[ Chikaterasu ] 2. Analyze first three eigenvectors ..."
#gmx anaeig -s ../md_target.tpr -f ../md_fit.xtc -first 1 -last 3 -extr ./extreme_pcavec.pdb -nframes 30 -entropy -n pca_target.ndx

echo "[ Chikaterasu ] 3. Plot first 3 PCA vectors as 2D maps ..."
#gmx anaeig -s ../md_target.tpr -f ../md_fit.xtc -2d ./2dproj_12.xvg -first 1 -last 2 -n pca_target.ndx
