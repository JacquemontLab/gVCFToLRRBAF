#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=75010M
#SBATCH --time=03:00:00
#SBATCH --account=rrg-jacquese
#SBATCH --nodes=1

##########################################################################################
#*****************************************************************************************
#
# Script: filter_merge_missing.sh
# Author: Mame Seynabou Diop
#
#*****************************************************************************************
##########################################################################################


# Description:
# - Applies per-chromosome filters: geno, mind, hwe, maf (MAF ≥ 1/1000), on SNP-only data
# - Merges all chromosomes using plink2 --pmerge-list
# - Computes missing data statistics (--missing) on the merged dataset. This one will help to calulate PFB file
##########################################################################################

# Load required modules
module load StdEnv/2020
module load plink/2.00a3.6


# Define directories
workDir=/home/$USER/SPARKpVCF
inputDir=$workDir/pvcfToPlink
filteredDir=$workDir/filtered_chr
mergedDir=$filteredDir/merged
statsDir=$mergedDir/stats

mkdir -p $filteredDir $mergedDir $statsDir

# Step 1: Filter each chromosome
for chr in {1..22}; do
  echo "Filtering chromosome $chr..."

  plink2 \
    --bfile $inputDir/wgs_12509_genome.deepvariant.chr${chr} \
    --memory 75000 \
    --geno 0.05 \
    --mind 0.05 \
    --hwe 0.000001 \
    --maf 0.001 \
    --make-bed \
    --out $filteredDir/wgs_12509_genome.deepvariant.chr${chr}.filtered
done


# Step 2: Create list.txt of all chromosomes except chr1 (for merge)
list_file="list.txt"
> "$list_file"

for chr in {2..22}; do
    prefix="$filteredDir/wgs_12509_genome.deepvariant.chr${chr}.filtered"
    echo -e "${prefix}.bed\t${prefix}.bim\t${prefix}.fam" >> "$list_file"
done

# Step 3: Merge chromosomes (starting from chr1)
plink2 \
  --bfile $filteredDir/wgs_12509_genome.deepvariant.chr1.filtered \
  --pmerge-list list.txt \
  --make-bed \
  --out $mergedDir/wgs_12509_genome.deepvariant.AllCHR


# Step 4: Compute missing data statistics on merged dataset
plink2 \
  --bfile $mergedDir/wgs_12509_genome.deepvariant.AllCHR \
  --memory 75000 \
  --missing \
  --out $statsDir/wgs_12509_genome.deepvariant.AllCHR.missing
