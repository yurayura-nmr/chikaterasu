#!/bin/sh

: '
*************************************************************
Chikaterasu Analysis

Erik Walinda
Kyoto University
Graduate School of Medicine

Last change: 2022-12-14
Added dipole (untested)
to add: gmx cluster -f md_protein.xtc -s md_protein.tpr -method gromos -cl cluster.pdb -g cluster.log -cutoff 0.2

gmx version         2021
chikaterasu version dev

Notes to self:

* signed principal axis may not correctly work with 3DNA for now ... (does not always find Helix axis?)
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

# == Most important parameters (most often set wrong :)) ==
dna=false            # DNA or protein?
fit=true            # rot_trans fit before analysis. AVOID for RheoMD if studying alignment
pbc=nojump          # mol or nojump

# == Where are the data? How many runs? ==
proc_folder="md_"   # PREFIX for folders, e.g. md for md_1, md_2, ...
nruns=1             # total number of runs (3 ~ 20)
anatime=0           # analyze frames starting at time t [ps]
on_the_fly=false    # analyze a currently active run?

# == What to analyze? ==
rmsd=true
rmsf=true
traj=true
dt=10               # time interval for various analysis functions [ps]. E.g. traj, hbond
vmd=false
gyration=false      # Calculate R_gyr
hbond=false
distance=false
sasa=false
pca=false
contactmap=false    # Draw contact map (mean-smallest-distance map)
dipole=false
mu=2.273            # Dipole moment of water. spc: 2.273 | tip4p/2005: 2.38


hbond_ATP=false     # Custom function for ATP research; have to specify group number "13" instead of name "ATP" bug?
domain_angle=false  # Custom function for diUb research

paxis=false         # Use gromacs to calculate principal axis. This can flip the axis since vector is unsigned.
signedPaxis=false   # Do we want to get the signed vector of the principal axis? Requires helper script [chika_paxis.py] 

dssp=false


: '
*************************************************************
Extra analysis options for dimers, multimers, or multiple chains?

If multiple similar domains, it might be difficult to get a consistent RMSD etc., 
since for each timestep the fitting to the initial frame structure will be different.

indexFileProvided: 
Are you providing a file with the atom numbers for analysis? [chika.ndx]
(if yes: provide a chika.ndx file in main folder). This specifies what atom groups to analyze
*************************************************************
'
indexFileProvided=false     
nchains=3                   # number of chains in the chika.ndx file 


: '
*************************************************************
For distance calculation.
Specify the two atoms between which the distance
is to be calculated from the topology file.
Splitting into domain of a multidomain protein
is not necessary for this.

Pair 1 is atom1-atom2 
Pair 2 is atom3-atom4
*************************************************************
'
atom1=701
atom2=1931
atom3=701
atom4=1931

: '
*************************************************************
For domain angle analysis

Not very sophisticated code at the moment.
Specify 0 / 1 for the distal and proximal subunits.

This means:

0    the first set of atoms in the PDB, e.g. atoms 1 - 1230
1    the second set of atoms in the PDB,e.g. atoms 1231 - ...

distal    GLQ is the identifier for distal Ub
              except in linear diUb, where it is GLY
proximal  LYQ is the identifier for proximal Ub
              except in linear diUb, where it is MET
*************************************************************
'
lineardiub=false    # custom for domain angle. Proximal/distal treatment differs.
distal_ID=0 
proximal_ID=1 

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
    cp ./runs/md/md.xtc ./results/$proc_folder/
    cp ./runs/md/md.tpr ./results/$proc_folder/
    cd ./results/$proc_folder  
  fi
  
  if [ "$on_the_fly" = false ] ; then
    cp ./runs/$proc_folder$i/md.xtc ./results/$proc_folder/
    cp ./runs/$proc_folder$i/md.tpr ./results/$proc_folder/
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
  gmx editconf -f md.tpr -o target.pdb 
  printf "non-Water\nq\n" | gmx make_ndx -f target.pdb -o target.ndx

  : '
  *************************************************************
  Make TPR file of system
  [non-water atoms]
  This is useful for keeping cofactors such as ZN ions in the
  analysis. Choosing protein only here is also OK, but it
  discards such cofactors.
  *************************************************************
  '
  printf "non-Water" | gmx convert-tpr -s md.tpr -n target.ndx -o md_target.tpr

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
  printf "non-Water" | gmx trjconv -s md_target.tpr -f md.xtc -n target.ndx -o md_target.xtc

  mv md.xtc md_full.xtc  # original file still contains all atoms (solvent, ions)


  : '
  *************************************************************
  Remove PBC
  This may be tricky for DNA and protein-prtotein complexes.
  See old version of nicoterasu for details.
  For now, only simple protein behaviour is implemented.
  *************************************************************
  '
  if [ "$dna" = true ] ; then
      printf "DNA\nSystem" | gmx trjconv -s md_target.tpr -f md_target.xtc -center -ur compact -pbc $pbc -o md_target_centered_no_PBC.xtc
  fi
  if [ "$dna" = false ] ; then
      printf "Protein\nSystem" | gmx trjconv -s md_target.tpr -f md_target.xtc -center -ur compact -pbc $pbc -o md_target_centered_no_PBC.xtc
  fi


  : '
  *************************************************************
  Rot-trans fit
  Not good for Rheo-MD, since we want to study alignment
  For normal protein studies, it should be no problem to default this.
  For DNA, sometimes tricky, see nicoterasu for details
  *************************************************************
  '
  if [ "$fit" = true ] ; then
      if [ "$dna" = true ] ; then
          printf "DNA\nSystem" | gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit rot+trans -o md_fit.xtc
      fi
      if [ "$dna" = false ] ; then
          printf "Backbone\nSystem" | gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit rot+trans -o md_fit.xtc
      fi
  fi
  if [ "$fit" = false ] ; then
      if [ "$dna" = true ] ; then
          printf "DNA\nSystem" | gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit trans -o md_fit.xtc
      fi
      if [ "$dna" = false ] ; then
          printf "Backbone\nSystem" | gmx trjconv -s md_target.tpr -f md_target_centered_no_PBC.xtc -fit trans -o md_fit.xtc
      fi
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
  RMSD
  C-alpha atoms
  DNA version not implemeted yet.
  *************************************************************
  '
  if [ "$rmsd" = true ] ; then
      if [ "$indexFileProvided" = false ] ; then
          printf "C-alpha\nC-alpha" | gmx rms -s md_target.tpr -f ./md_fit.xtc -o ./rmsd/rmsd.xvg -fit rot+trans
      fi
      if [ "$indexFileProvided" = true ] ; then
          for k in `seq 1 $nchains`;
          do
              gmx rms -s md_target.tpr -f ./md_fit.xtc -o ./rmsd/rmsd_$k.xvg -n ../../chika.ndx -fit rot+trans
          done
      fi
   fi

  : '
  *************************************************************
  RMSF
  C-alpha atoms
  DNA version not implemeted yet.
  *************************************************************
  '
  if [ "$rmsf" = true ] ; then
      if [ "$indexFileProvided" = false ] ; then
          printf "C-alpha\n" | gmx rmsf -s md_target.tpr -f ./md_fit.xtc -o ./rmsf/rmsf.xvg -fit -res -b $anatime
      fi
      if [ "$indexFileProvided" = true ] ; then
          for j in `seq 1 $nchains`;
          do
              gmx rmsf -s md_target.tpr -f ./md_fit.xtc -o ./rmsf/rmsf_$j.xvg -n ../../chika.ndx -fit -res -b $anatime
          done
      fi
  fi

  : '
  *************************************************************
  PDB frames
  Extract PDB frames for pymol analysis.
  Choose dt wisely or the file will be huge.
  Splitting the chains for pymol analysis not implemented yet.
  *************************************************************
  '
  if [ "$traj" = true ] ; then
      printf "System" | gmx trjconv -f md_fit.xtc -s md_target.tpr -o ./traj/traj.pdb -dt $dt # -pbc mol -pbc nojump
  fi

  : '
  *************************************************************
  VMD output
  If we want to open the trajectory in VMD, we need a gro
  file of our target (system without water) topology.
  *************************************************************
  '
  if [ "$vmd" = true ] ; then
      cp md_fit.xtc ./vmd/ 
      gmx editconf -f md_target.tpr -o ./vmd/md.gro
  fi

  : '
  *************************************************************
  Hydrogen bonding analysis
  This part of the code is not very general.
  In the future optimization is required.
  For example, could specify in the option section:
  * Protein/Water
  * Main-chain/Main-chain
  etc. 
  And then check by if statements, which one to do
  *************************************************************
  '
  if [ "$hbond" = true ] ; then
      cd hbond
      printf "MainChain+H\nMainChain+H" | gmx hbond -f ../md_target.xtc -s ../md_target.tpr -hbm -hbn hbond.ndx -num backbone_backbone.xvg -dt $dt
      
      # *************************************************************
      # Extract specific hydrogen bonds in table-form using Justin's script
      # * Uses the accessory script plot_hbmap.pl. The script is provided in the custom_analysis directory.
      # * Requires that plot_hbmap.pl is in the PATH and is executable. At present, Chikaterasu does not automatically copy that script here.
      #   (so for now, the easiest fix may be to copy plot_hbmap.pl into /usr/local/bin and make it executable.
      
      # [Caution]
      # * Syntax of plot_hbmap.pl requires that the .ndx file (e.g., analyze.ndx) # contains only the [hbonds...] section (e.g. [ hbonds_MainChain+H ]) and only the atom numbers herein.
      #   (so, there should be no other sections in this file; this would confuse the final output).
      # * To achieve that, we use awk to so that only this section is in the data is in the ndx that will be the input for the perl script.

      # read -p "[Chikaterasu-dev] Make sure that plot_hbmap.pl is in the PATH and executable." dummy
      awk '/hbond/{y=1;next}y' hbond.ndx > analyze.ndx
      plot_hbmap.pl -s ../target.pdb -map hbmap.xpm -index analyze.ndx

      # *************************************************************

      # Continue with other hydrogen bond analyses ...
      printf "Protein\nProtein" | gmx hbond -f ../md_target.xtc -s ../md_target.tpr -hbm -hbn -num protein_protein.xvg -dt $dt

      #printf "Protein\nWater" | gmx hbond -f md_full.xtc -s md.tpr -num hbond/protein_water.xvg
      if [ "$hbond_ATP" = true ] ; then
          printf "Protein\n13" | gmx hbond -f ../md_full.xtc -s ../md.tpr -num protein_ATP.xvg
          printf "Water\n13" | gmx hbond -f ../md_full.xtc -s ../md.tpr -num ATP_water.xvg
      fi
      cd ..
  fi

  : '
  *************************************************************
  Contact map
  Untested
  *************************************************************
  '
  if [ "$contactmap" = true ] ; then
      cd contact_map 

      # Make xpm file
      printf "MainChain+H" | gmx mdmat -f ../md_target.xtc -s ../md_target.tpr -dt $dt
      
      # Convert to eps for viewing
      gmx xpm2ps -f dm.xpm
                          
      cd ..
  fi

  : '
  *************************************************************
  Radius of gyration
  *************************************************************
  '
  if [ "$gyration" = true ] ; then
      if [ "$dna" = false ] ; then
          printf "Protein\n" | gmx gyrate -s md_target.tpr -f md_fit.xtc -o gyration/gyr.xvg -b $anatime
      fi
      if [ "$dna" = true ] ; then
          printf "DNA\n" | gmx gyrate -s md_target.tpr -f md_fit.xtc -o gyration/gyr.xvg -b $anatime
      fi
  fi

  : '
  *************************************************************
  DSSP
  *************************************************************
  '
  if [ "$dssp" = true ] ; then
      cd dssp
      if [ "$indexFileProvided" = false ] ; then
          gmx do_dssp -s ../md_target.tpr -f ../md_fit.xtc -b $anatime 
      fi
      if [ "$indexFileProvided" = true ] ; then
          gmx do_dssp -s ../md_target.tpr -f ../md_fit.xtc -n ../../../chika.ndx  -b $anatime 
      fi
      #exit 1
      cd ..
  fi

  : '
  *************************************************************
  Principal axes of inertia
  
  Principal axis can be plotted with the script: 
  md_plot_paxis_tip_overtime.py
  
  or MD analysis (python).
  If using GROMACS, need to remove header and time column.
  *************************************************************
  '
  if [ "$paxis" = true ] ; then
      if [ "$dna" = false ] ; then
          printf "Protein\n" | gmx principal -s md_target.tpr -f md_fit.xtc -a1 paxis/major_principal_axis.xvg -b $anatime 
      fi
      if [ "$dna" = true ] ; then
          printf "DNA\n" | gmx principal -s md_target.tpr -f md_fit.xtc -a1 paxis/major_principal_axis.xvg -b $anatime 
      fi
   fi

  : '
  *************************************************************
  SASA
  Calculate the surface accessible area for all residues of
  the protein and its standard deviation.
  The order of chains in the output file should be the same
  as in the topology, i.e. GLQ: distal, LYQ: proximal.
  *************************************************************
  '
  if [ "$sasa" = true ] ; then
      printf "Protein\n" | gmx sasa -s md_target.tpr -f md_fit.xtc -or sasa/sasa.xvg -b $anatime 
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

      printf "0\n" | gmx distance -f md_fit.xtc -s md_target.tpr -n ./distance/chikaterasu.ndx -oall ./distance/distance.xvg

      mv ./distance/distance.xvg ./distance/distance_1.xvg
      rm ./distance/chikaterasu.ndx

      echo "[ Chikaterasu Distance $atom3 to Residue $atom4 ]" > ./distance/chikaterasu.ndx
      echo "$atom3 $atom4" >> ./distance/chikaterasu.ndx

      printf "0\n" | gmx distance -f md_fit.xtc -s md_target.tpr -n ./distance/chikaterasu.ndx -oall ./distance/distance.xvg
      mv ./distance/distance.xvg ./distance/distance_2.xvg
  fi


  : '
  *************************************************************
  Dipole analysis of solvent (untested)
  *************************************************************
  : '
  if [ "$dipole" = true ] ; then
      #printf "0\n" | gmx distance -f md_fit.xtc -s md_target.tpr -n ./distance/chikaterasu.ndx -oall ./distance/distance.xvg
      gmx dipoles -P 1 -o ./dipole/dip_test -mu $mu -mumax 5.0 -f md_full.xtc -s md.tpr
  fi
  
  : '
  *************************************************************
  Principal component analysis  
  *************************************************************
  : '
  if [ "$pca" = true ] ; then
      echo "[ Chikaterasu ] PCA analysis will take a while ..."
      echo "[ Chikaterasu ] 1. Building covariance matrix ..."
      
      printf "Backbone\nBackbone" | gmx covar -s md_target.tpr -f md_fit.xtc -o pca/eigenval.xvg -av pca/average.pdb -b $anatime
      
      echo "[ Chikaterasu ] 2. Analyze first three eigenvectors ..."
      printf "Backbone\nBackbone" | gmx anaeig -s md_target.tpr -f md_fit.xtc -first 1 -last 3 -extr pca/extreme_pcavec.pdb -nframes 30 -entropy
              
      echo "[ Chikaterasu ] 3. Plot first 3 PCA vectors as 2D maps ..."
      printf "Backbone\nBackbone" | gmx anaeig -s md_target.tpr -f md_fit.xtc -2d pca/2dproj_12.xvg -first 1 -last 2
      printf "Backbone\nBackbone" | gmx anaeig -s md_target.tpr -f md_fit.xtc -2d pca/2dproj_13.xvg -first 1 -last 3
      
      #exit 1 # tested OK
  fi

  : '
  *************************************************************
  Signed principal axes analysis
  *************************************************************
  : '
  if [ "$signedPaxis" = true ] ; then
      
      # loop over trajectory
      dt="10"                 # ideally 10 for finest sampling
      steps="20000"           # ideally 20000 for 0.2 microsecond

      # copy helper script
      cp ../../chika_paxis.py ./paxis/
      cp ../../chika_paxis_DNA.py ./paxis/

      if [ "$dna" = false ] ; then
          for j in `seq 0 $steps`;
          do
            rm ./paxis/\#*
            echo "Step $j of $steps"
            currentFrame=$((j * dt))
            printf "Protein \nq\n" | gmx trjconv -s md_target.tpr -f md_fit.xtc -o ./paxis/test.pdb -b $currentFrame -e $currentFrame

            cd paxis
            python chika_paxis.py >> output.txt
            cd ..
            #exit 1
            done
      fi
      if [ "$dna" = true ] ; then
          for j in `seq 0 $steps`;
          do
            rm ./paxis/\#*
            echo "Step $j of $steps"
            currentFrame=$((j * dt))
            printf "DNA \nq\n" | gmx trjconv -s md_target.tpr -f md_fit.xtc -o ./paxis/test.pdb -b $currentFrame -e $currentFrame
            #rm *moi*
            #rm *paxis*
            #printf "DNA\n" | gmx principal -s md_target.tpr -f md_fit.xtc -a1 paxis/major_principal_axis.xvg -b $currentFrame -e $currentFrame
            cd paxis
            #python chika_paxis_DNA.py >> output.txt
            
            find_pair test.pdb bpfile.dat
            analyze bpfile.dat

            # Not in all frames will 3DNA find the helix axis (too skewed?)
            if grep -q "Helix:" test.out; then
                echo -n " " $currentFrame "  " >> output.txt
            else
                echo Could not find helix axis in frame $currentFrame
            fi
            grep "Helix:" test.out >  a.txt
            awk '{print $2 "  " $3 "  " $4}' a.txt >> output.txt
            cd ..
            #exit 1
          done
      fi
  fi
  
  : '
  *************************************************************
  Domain angle analysis
  Not very sophisticated code.
  Have to specify in the header, which GROMACS atom group is
  proximal Ub and which one is distal Ub.
  Currently, using [del 0-14] so that I have to specify only
  0 / 1 for distal proximal.
 
  untested ... continue here!

  Currently testing (2018-02-27) ... hope this goes well!
  after output.txt is written, do the rest (plotting, converg.
  outside of chikaterasu; plot in dropbox etc.)
  *************************************************************
  '
  if [ "$domain_angle" = true ] ; then
      if [ "$lineardiub" = false ] ; then
          printf "splitch 0\ndel 0-14\nq\n" | gmx make_ndx -f md_target.tpr -o ./domain_angle/test.ndx
      fi
      if [ "$lineardiub" = true ] ; then
          printf "r 1-76\nr 77-152\ndel 0-14\nq\n" | gmx make_ndx -f md_target.tpr -o ./domain_angle/test.ndx
      fi

      # static structures - for testing
      printf "$proximal_ID \nq\n" | gmx editconf -f md_target.tpr -n ./domain_angle/test.ndx -o ./domain_angle/proximal_initial.pdb
      printf "$distal_ID \nq\n" | gmx editconf -f md_target.tpr -n ./domain_angle/test.ndx -o ./domain_angle/distal_initial.pdb

      #exit 1
      : '
      Using python here to do the domain angle fitting.

      Reference   Proximal
      Mobile      Distal

      Thus, we are superimposing distal onto proximal Ub.
      Using only residues 1-70, i.e. rigid part of Ub.
      '


      # Make python script

      VAR1=$(cat <<EOF
from math import *
import transformations

from MDAnalysis import *
from MDAnalysis.analysis.align import *
from MDAnalysis.analysis.rms import rmsd

ref     = Universe('./proximal.pdb')
mobile  = Universe('./distal.pdb')

mobile0 = mobile.atoms.CA.positions - mobile.atoms.center_of_mass()
ref0    = ref.atoms.CA.positions    - ref.atoms.center_of_mass()

mobileR = mobile0[0:70]
refR    = ref0[0:70]

R, rmsd = rotation_matrix(mobileR, refR)

#print R
q = transformations.quaternion_from_matrix(R)

print rmsd, "\t", degrees(2*np.arccos(q[0]))
EOF
      )

      echo "${VAR1}" > ./domain_angle/chika_domain_angle.py
      cp ../../transformations.py ./domain_angle/
      #exit 1

      # loop over trajectory
      dt="10"                 # ideally 10 for finest sampling
      steps="10000"           # ideally 10000 for 0.1 microsecond

      for j in `seq 0 $steps`;
      do
        rm ./domain_angle/\#*
        echo "Step $j of $steps"
        currentFrame=$((j * dt))
        printf "$proximal_ID \nq\n" | gmx trjconv -s md_target.tpr -f md_fit.xtc -n ./domain_angle/test.ndx -o ./domain_angle/proximal.pdb -b $currentFrame -e $currentFrame
        printf "$distal_ID \nq\n" | gmx trjconv -s md_target.tpr -f md_fit.xtc -n ./domain_angle/test.ndx -o ./domain_angle/distal.pdb -b $currentFrame -e $currentFrame

        cd domain_angle
        #cp ../../chika_domain_angle.py .
        python chika_domain_angle.py >> output.txt
        cd ..
        #exit 1
      done
  fi
  #exit 1

  : '
  *************************************************************
  Store folder correctly
  *************************************************************
  '
  cd ../..
  mv results/$proc_folder results/$proc_folder$i

  echo "[Chikaterasu] Run $i of $nruns finished!"

done
