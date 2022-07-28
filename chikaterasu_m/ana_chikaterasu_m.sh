#!/bin/sh

: '
*************************************************************
Chikaterasu Analysis

Erik Walinda
Kyoto University
Graduate School of Medicine

Last change: 2022-3-31

gmx version         2021
chikaterasu version martini-dev
*************************************************************
'

: '
*************************************************************
Manually setup parameters for the run
What will be analyzed?
How many runs? 
etc.
*************************************************************
'

cleanup=true
sample="noncyclic_ea1_eb1"

# == Most important parameters (most often set wrong :)) ==
dna=false           # DNA or protein?
fit=true            # rot_trans fit before analysis. AVOID for RheoMD if studying alignment. Must be true at the moment.
#pbc=               # fixed for martini

# == Where are the data? How many runs? ==
proc_folder="md_"   # PREFIX for folders, e.g. md for md_1, md_2, ...
nruns=1             # total number of runs (3 ~ 20)
anatime=0           # analyze frames starting at time t [ps]
on_the_fly=true     # analyze a currently active run?

# == What to analyze? ==
# Implement step by step
#rmsd=true
#rmsf=true
traj=true
dt=5000              # time interval for various analysis functions [ps]. E.g. traj, hbond
#vmd=false
#gyration=false      # Calculate R_gyr
#hbond=false
distance=true
#sasa=false
#pca=false
#contactmap=false    # Draw contact map (mean-smallest-distance map)
#dipole=false
#mu=2.273            # Dipole moment of water. spc: 2.273 | tip4p/2005: 2.38


# == Dimers, hexamers or multiple chains? ==
indexFileProvided=false     # give a chika.ndx file in main folder to specify what to analyze
nchains=1                   # number of chains in the chika.ndx file

: '
*************************************************************
For distance calculation.
Specify the two atoms between which the distance
is to be calculated from the topology file.
Splitting into domain of a multidomain protein
is not necessary for this.
*************************************************************
'
#distance_name="I44_I20_AC1"
atom1=93
atom2=256

: '
"""
Define sample name (easier to keep track of on-the-fly processed filenames)
"""
'
valhost=$(hostname)"_"$sample"_traj.pdb"
valdist=$(hostname)"_"$sample"_distance.xvg"

: '
*************************************************************
Setup directories for the run
*************************************************************
'

read -p "[Chikaterasu-dev] Starting analysis." dummy

if [ "$cleanup" = true ] ; then
    read -p "[Chikaterasu-dev] Chikaterasu will clean up the results folder to save disk space! Abort if necessary." dummy
    rm -rf results
fi


: '
*************************************************************
Start the analysis
Loops over all data
*************************************************************
'

for i in `seq 1 $nruns`;

do

  mkdir -p results
  mkdir -p results/$proc_folder
  mkdir -p results/$proc_folder/rmsd
  mkdir -p results/$proc_folder/rmsf
  mkdir -p results/$proc_folder/traj
  mkdir -p results/$proc_folder/dssp
  mkdir -p results/$proc_folder/energy
  mkdir -p results/$proc_folder/gyration
  mkdir -p results/$proc_folder/mindist
  mkdir -p results/$proc_folder/sasa
  mkdir -p results/$proc_folder/hbond
  mkdir -p results/$proc_folder/pca
  mkdir -p results/$proc_folder/vmd
  mkdir -p results/$proc_folder/distance
  mkdir -p results/$proc_folder/domain_angle
  mkdir -p results/$proc_folder/paxis
  mkdir -p results/$proc_folder/contact_map
  mkdir -p results/$proc_folder/dipole

  : '
  *************************************************************
  Copy trajectory
  Go to folder of interest
  *************************************************************
  '
  if [ "$on_the_fly" = true ] ; then
    read -p "[Chikaterasu-dev] Chikaterasu will try to analyze a currently on-going run! Abort if necessary." dummy
    cp ./runs/md/dynamic.xtc ./results/$proc_folder/
    cp ./runs/md/dynamic.tpr ./results/$proc_folder/
    cd ./results/$proc_folder
  fi

  if [ "$on_the_fly" = false ] ; then
    cp ./runs/$proc_folder$i/dynamic.xtc ./results/$proc_folder/
    cp ./runs/$proc_folder$i/dynamic.tpr ./results/$proc_folder/
    cd ./results/$proc_folder
  fi

  : '
  *************************************************************
  Make NDX file of system
  If multiple chains, have to split chains at this point
  See old version of nicoterasu for details

  Currently, a syntax error in make_ndx; but for non-Water
  this does not seem to matter, since group is alread there?
  *************************************************************
  '

  # Step by step
  # ---- Visualization, PBC etc ----
  # to see motion. Center on protein and prevent jumping (for k48)

  gmx editconf -f ./dynamic.tpr -o ./target.pdb

  # Select group 1 (Protein). Since I have multiple groups called protein, this will throw me an error if I select by string right now.
  printf "1\nq\n" | gmx make_ndx -f ./target.pdb -o ./target.ndx

  : '
  *************************************************************
  Make TPR file of system
  [non-water atoms]
  This is useful for keeping cofactors such as ZN ions in the
  analysis. Choosing protein only here is also OK, but it
  discards such cofactors.
  *************************************************************
  '
  # Select group 1 (Protein)
  printf "1\n1\n" | gmx convert-tpr -s ./dynamic.tpr  -n ./target.ndx -o ./md_target.tpr

  : '
  *************************************************************
  Make XTC file of system
  Create a XTC trajectory file of only the desired part of the
  system, e.g. only protein A or protein B
  This step went wrong in earlier versions (mayuterasu).
  The current implementation should work for 1 protein.
  DNA/protein complexes may require a different setup here
  Keep md_full for other analysis (e.g. protein-solvent IA)
  *************************************************************
  '

  : '
  *************************************************************
  Remove PBC
  This may be tricky for DNA and protein-prtotein complexes.
  See old version of nicoterasu for details.
  For now, only simple protein behaviour is implemented.
  
  dt untested
  *************************************************************
  '

  # Remove PBC
  printf "1\n1\n" | gmx trjconv -s ./dynamic.tpr -f ./dynamic.xtc -center -ur compact -pbc nojump -dt $dt -o ./md_target_centered_no_PBC.xtc
  #printf "1\n1\n" | gmx trjconv -s ./dynamic.tpr -f ./dynamic.xtc -center -ur compact -pbc nojump -o ./md_target_centered_no_PBC.xtc


  : '
  *************************************************************
  Rot-trans fit
  Not good for Rheo-MD, since we want to study alignment
  For normal protein studies, it should be no problem to default this.
  For DNA, sometimes tricky, see nicoterasu for details
  *************************************************************
  '
  if [ "$fit" = true ] ; then
      printf "1\n1\n" | gmx trjconv -s md_target.tpr -f ./md_target_centered_no_PBC.xtc -fit rot+trans -dt $dt -o ./md_fit.xtc
  fi

  : '
  *************************************************************
  All file preparations and manipulations done.
  Can start analysis now.
  *************************************************************
  '
  echo "[Chikaterasu] Run $i of $nruns : Starting analysis..."

  : '
  *************************************************************
  PDB frames
  Extract PDB frames for pymol analysis.
  Choose dt wisely or the file will be huge.
  Splitting the chains for pymol analysis not implemented yet.
  *************************************************************
  '
  echo "[Chikaterasu] Making PDB trajectory for viewing in Pymol etc."
  if [ "$traj" = true ] ; then
      printf "1\n1\n" | gmx trjconv -s md_target.tpr -f md_fit.xtc -fit rot+trans -o ./traj/$valhost -conect -dt $dt
  fi

  : '
  *************************************************************
  Distance analysis
  Atoms are specificed in the options section.
  *************************************************************
  : '
  if [ "$distance" = true ] ; then
      rm ./distance/chikaterasu.ndx
      echo "[ Chikaterasu Distance $atom1 to Residue $atom2 ]" > ./distance/chikaterasu.ndx
      echo "$atom1 $atom2" >> ./distance/chikaterasu.ndx

      printf "0\n" | gmx distance -f ./md_fit.xtc -s ./md_target.tpr -n ./distance/chikaterasu.ndx -oall ./distance/distance.xvg -tu ns -dt $dt
  fi

  : '
  *************************************************************
  Store folder correctly
  *************************************************************
  '
  cd ../..
  mv results/$proc_folder results/$proc_folder$i

  echo "[Chikaterasu] Run $i of $nruns finished!"

done

exit 1




# ----

# on the fly check distance

# make_ndx as usual
#gmx distance -f ../dynamic.xtc -s ../dynamic.tpr -oall $valdist -n index.ndx -rmpbc -tu ns

# print max. distance
#awk 'NR==2{max = $2 + 0; next} {if ($2 > max) max = $2;} END {print max}' $valdist

#or just
#tail -n 30 distance.xvg
