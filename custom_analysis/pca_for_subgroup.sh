gmx covar -s md_target.tpr -f md_fit.xtc -o pca/eigenval.xvg -av pca/average.pdb -n index.ndx

gmx anaeig -s md_target.tpr -f md_fit.xtc -first 1 -last 3 -extr pca/extreme_pcavec.pdb -nframes 30 -entropy -v eigenvec.trr  -n index.ndx

gmx anaeig -s md_target.tpr -f md_fit.xtc -2d pca/2dproj_12.xvg -first 1 -last 2 -n index.ndx
