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

# ----------------------------- Argument checks ------------------------------------
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_bfile_dir> <outputDir>"
    exit 1
fi

inputDir="$1"
outputDir="$2"

# Check input dir
if [[ ! -d "$inputDir" ]]; then
    echo "Error: input directory '$inputDir' does not exist."
    exit 1
fi

# Create output dir if needed
mkdir -p "$outputDir"


# ----------------------------- Find bfiles ----------------------------------------
# Get sorted list of all .bed files
bed_files=($(find "$inputDir" -maxdepth 1 -name "*.bed" | sort))


# Use the first .bed file as the base
base_bed="${bed_files[0]}"
base_prefix="${base_bed%.bed}"

# Load plink2 if needed
if ! command -v plink2 >/dev/null 2>&1; then
    echo "'plink2' not found — loading modules..."
    module load StdEnv/2020
    module load plink/2.00a3.6

    # Check again
    if ! command -v plink2 >/dev/null 2>&1; then
        echo "Error: plink2 is still not available after loading modules." >&2
        exit 1
    fi
else
    echo "plink2 is already available."
fi


# Get number of CPUs
cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running with $ncores cores"

# Get memory safely: prefer SLURM allocated memory
if [[ -n "$SLURM_MEM_PER_NODE" ]]; then
  # Use 90% of SLURM allocated memory (as safety margin)
  plink_mem=$(( SLURM_MEM_PER_NODE * 90 / 100 ))
  echo "Detected SLURM memory: $SLURM_MEM_PER_NODE MB"
else
  # Fallback to checking system memory
  read total_mem used_mem free_mem shared_mem buff_cache available_mem <<< $(free -m | awk '/Mem:/ {print $2, $3, $4, $5, $6, $7}')
  echo "Available memory (MB): $available_mem"

  # Use 90% of available memory
  plink_mem=$(( available_mem * 90 / 100 ))
fi

echo "Setting PLINK memory to: $plink_mem MB"
echo "Setting PLINK threads to: $cpus"


# Define directories
filteredDir=$outputDir/filtered_chr
mergedDir=$filteredDir/merged
statsDir=$mergedDir/stats

mkdir -p $filteredDir $mergedDir $statsDir

# ----------------------------------------------
# Loop over all .bed files in the input directory
# ----------------------------------------------
for bed_file in "$inputDir"/*.bed; do
  base_name=$(basename "$bed_file" .bed)
  prefix="$inputDir/$base_name"

  echo "Filtering $base_name..."

  plink2 \
    --bfile "$prefix" \
    --threads "$cpus" \
    --memory "$plink_mem" \
    --geno 0.05 \
    --mind 0.05 \
    --hwe 0.000001 \
    --maf 0.001 \
    --make-bed \
    --out "$filteredDir/${base_name}.filtered"
done


# Step 2: Create list.txt of all chromosomes except first vcf (for merge)
list_file="list.txt"
> "$list_file"


skip_first=true
for bed_file in "$filteredDir"/*.bed; do
  base_name=$(basename "$bed_file" .bed)
  prefix="$filteredDir/$base_name"

  if $skip_first; then
    first_vcf="$prefix"
    skip_first=false
    continue
  fi

  echo -e "${prefix}.bed\t${prefix}.bim\t${prefix}.fam" >> "$list_file"
done

# Step 3: Merge chromosomes (starting from first vcf)
plink2 \
  --bfile "$first_vcf" \
  --threads "$cpus" \
  --memory "$plink_mem" \
  --pmerge-list list.txt \
  --make-bed \
  --out $mergedDir/merged_dataset


# Step 4: Compute missing data statistics on merged dataset
plink2 \
  --bfile $mergedDir/merged_dataset \
  --threads "$cpus" \
  --memory "$plink_mem" \
  --missing \
  --out $statsDir/merged_dataset.missing
