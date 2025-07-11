#!/bin/bash
#SBATCH --ntasks=1                  
#SBATCH --mem=48000M              # Memory requested 
#SBATCH --time=2:30:00              # Maximum run time
#SBATCH --account=rrg-jacquese      # Cluster account
#SBATCH --nodes=1                   # Number of nodes

##########################################################################################
#*****************************************************************************************
#
# Script name: pVCF_to_plink.sh      Author: Mame Seynabou Diop
#
#*****************************************************************************************
##########################################################################################

# --------------------------------------------------------------------
# Convert a pVCF to PLINK binary format (BED/BIM/FAM) for a single
# chromosome, applying a simple depth (DP) filter.
# --------------------------------------------------------------------

# Get the chromosome number passed as a command-line argument
chr=$1

# Working directory (created if it does not exist)
workDir=/home/$USER/SPARKpVCF
mkdir -p "$workDir"

# Directory containing the gzipped pVCF files
inputDir=/home/$USER/projects/rrg-jacquese/All_user_common_folder/RAW_DATA_old/SPARK_V2_V3-2021_10/iWGS_v1_gVCF/deepvariant/pvcf

# Output directory for the PLINK files
outputDir=$workDir/pvcfToPlink
mkdir -p "$outputDir"

# Load required modules
module load StdEnv/2020
module load plink/2.00a3.6

# ---------------------------------------
# Convert the pVCF to PLINK binary format
# ---------------------------------------
plink2 \
  --vcf $inputDir/wgs_12509_genome.deepvariant.chr"$chr".vcf.gz \
  --memory 47000 \
  --snps-only \
  --max-alleles 2 \
  --vcf-min-dp 1 \
  --vcf-half-call "m" \
  --make-bed \
  --out $outputDir/wgs_12509_genome.deepvariant.chr"$chr"



# Explanation of PLINK2 options:
# --vcf             = input pVCF file for the specified chromosome
# --memory          = memory allocated to PLINK in megabytes (MB)
# --snps-only       = keep only SNPs, exclude indels and other variant types
# --max-alleles 2   = keep only biallelic variants (variants with exactly 2 alleles)
# --vcf-min-dp 1    = exclude genotypes with depth (DP) less than 1
# --vcf-half-call m = treat half-called genotypes (e.g., 0/.) as missing ("m")
# --make-bed        = generate binary PLINK output files: .bed, .bim, .fam
# --out             = prefix used for naming output files
