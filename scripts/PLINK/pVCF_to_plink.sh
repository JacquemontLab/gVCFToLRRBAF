#!/bin/bash

##########################################################################################
#*****************************************************************************************
#
# Script name: pVCF_to_plink.sh      Author: Mame Seynabou Diop
#
#*****************************************************************************************
##########################################################################################

# --------------------------------------------------------------------
# Convert a pVCF to PLINK binary format (BED/BIM/FAM) applying a simple depth (DP) filter.
# --------------------------------------------------------------------

# Explanation of PLINK2 options used:
# --vcf             = input pVCF file
# --memory          = memory allocated to PLINK in megabytes (MB)
# --snps-only       = keep only SNPs, exclude indels and other variant types
# --max-alleles 2   = keep only biallelic variants (variants with exactly 2 alleles)
# --vcf-min-dp 1    = exclude genotypes with depth (DP) less than 1
# --vcf-half-call m = treat half-called genotypes (e.g., 0/.) as missing ("m")
# --make-bed        = generate binary PLINK output files: .bed, .bim, .fam
# --out             = prefix used for naming output files


# Usage check
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf.gz> <output_dir>"
    exit 1
fi


# Get command-line arguments
inputVcf="$1"
outputDir="$2"

# Check input file
if [ ! -f "$inputVcf" ]; then
    echo "Error: Input file '$inputVcf' does not exist." >&2
    exit 1
fi

# Create output directory if needed
mkdir -p "$outputDir"

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

outPrefix=$(basename "$inputVcf" .vcf.gz)



# Get number of CPUs
cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running with $cpus cores"

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

# ---------------------------------------
# Convert the pVCF to PLINK binary format
# ---------------------------------------
plink2 \
  --vcf "$inputVcf" \
  --threads "$cpus" \
  --memory "$plink_mem" \
  --snps-only \
  --max-alleles 2 \
  --vcf-min-dp 1 \
  --vcf-half-call "m" \
  --make-bed \
  --out "$outputDir/$outPrefix"

