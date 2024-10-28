#!/bin/bash
# Script to compare two sets of clusters from different MD simulations
# Purpose: Identify if clusters from simulation A are present in simulation B, and vice versa.

# Step 1: Preparation - Generate cluster files with the script cluster_analysis.sh.

# Step 2: Split each cluster file into individual model files for RMSD comparisons.
# Usage:
# awk '/^MODEL/{i++}{print > "clusterA_model_"i".pdb"}' clustersA.pdb
# awk '/^MODEL/{i++}{print > "clusterB_model_"i".pdb"}' clustersB.pdb

# Define the number of clusters for each simulation

num_clusters_in_A=10
num_clusters_in_B=10

# Disable GROMACS backup creation for temporary files
export GMX_MAXBACKUP=-1

# Initialize output file for the RMSD matrix
echo "# RMSD Matrix for Cluster Comparison" > rmsd_matrix.xvg

# Main loop: Perform pairwise RMSD calculations between all clusters in A and B
for i in $(seq 1 $num_clusters_in_A); do
  row=""
  for j in $(seq 1 $num_clusters_in_B); do
    # Calculate RMSD between cluster i from A and cluster j from B using C-alpha atoms
    gmx rms -f clusterA_model_$i.pdb -s clusterB_model_$j.pdb -o rmsd_tmp.xvg <<EOF
3
3
EOF

    # Extract the last RMSD value from rmsd_tmp.xvg
    rmsd=$(tail -n 1 rmsd_tmp.xvg | awk '{print $2}')

    # Ensure RMSD value was captured
    if [[ -n "$rmsd" ]]; then
      row="$row $rmsd"
    else
      row="$row NaN"  # Mark as NaN if extraction failed
    fi

    # Remove temporary file to avoid backup accumulation
    rm -f rmsd_tmp.xvg
  done
  # Append row to matrix
  echo "$row" >> rmsd_matrix.xvg
done
