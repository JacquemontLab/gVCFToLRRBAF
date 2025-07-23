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
# This script submits the pVCF_to_plink.sh job separately by pVCF files present in outputDir
# It is useful when dealing with large pVCF files per chromosome,
# allowing parallel or sequential processing on a cluster.
#
# Usage:
#   ./launch_allpVCF.sh
#
# Note:
#   Make sure that 'pVCF_to_plink.sh' is executable and 
#   accessible in the same directory or in your PATH.
# ---------------------------------------------------------

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

inputDir=/home/$USER/projects/rrg-jacquese/flben/WGS_pipeline/test/data/subset_10_samples
outputDir=/home/$USER/projects/rrg-jacquese/flben/WGS_pipeline/test/data/plink

# Loop over all .vcf.gz files
find "$inputDir" -type f -name "*.vcf.gz" | while read -r vcf_file; do
  echo "💻 Submitting job for file: $vcf_file"
  
  # sbatch \
  # --ntasks=1 \
  # --cpus-per-task=1 \
  # --mem-per-cpu=3500MB \
  # --time=00:10:00 \
  # --account=rrg-jacquese \
  "$SCRIPT_DIR/pVCF_to_plink.sh" "$vcf_file" "$outputDir"
done











vcf_file=/home/flben/projects/rrg-jacquese/flben/WGS_pipeline/test/data/subset_10_samples_v2/wgs_12509_genome.deepvariant.chr1.vcf.gz