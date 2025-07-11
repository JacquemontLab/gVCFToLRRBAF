#!/bin/bash

##########################################################################################
#*****************************************************************************************
#
# Script name: launch_allpVCF.sh       
# Author: Mame Seynabou Diop
#
#******************************************************************************************
############################################################################################

# Description:
# This script submits the pVCF_to_plink.sh job separately
# for each chromosome (1 to 22, X, Y).
# It is useful when dealing with large pVCF files per chromosome,
# allowing parallel or sequential processing on a cluster.
#
# Usage:
#   ./launch_all.sh
#
# Note:
#   Make sure that 'pVCF_to_plink.sh' is executable and 
#   accessible in the same directory or in your PATH.
# ---------------------------------------------------------

# Loop through autosomes 1 to 22 and sex chromosomes X and Y
for chr in {1..22} X Y; do
  echo "Submitting job for chromosome $chr"
  sbatch pVCF_to_plink.sh $chr
done
