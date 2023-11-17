# First, make the target by getting all backbone atoms and then excluding unwanted parts of the topology manually
# gmx make_ndx -f ../md.tpr

echo "[ Chikaterasu ] PCA analysis will take a while ..."
echo "[ Chikaterasu ] 1. Building covariance matrix ..."

gmx covar -s ../md_target.tpr -f ../md_fit.xtc -o ./eigenval.xvg -av ./average.pdb -n pca_target.ndx

echo "[ Chikaterasu ] 2. Analyze first three eigenvectors ..."
#gmx anaeig -s ../md_target.tpr -f ../md_fit.xtc -first 1 -last 3 -extr ./extreme_pcavec.pdb -nframes 30 -entropy -n pca_target.ndx

echo "[ Chikaterasu ] 3. Plot first 3 PCA vectors as 2D maps ..."
#gmx anaeig -s ../md_target.tpr -f ../md_fit.xtc -2d ./2dproj_12.xvg -first 1 -last 2 -n pca_target.ndx
