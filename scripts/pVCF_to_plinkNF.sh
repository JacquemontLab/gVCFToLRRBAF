#!/bin/bash

# === Usage check ===
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf.gz> <output_dir>"
    exit 1
fi

inputVcf="$1"
outputDir="$2"

# === Check input file ===
if [ ! -f "$inputVcf" ]; then
    echo "Error: Input file '$inputVcf' does not exist." >&2
    exit 1
fi

# === Load plink2 if needed ===
if ! command -v plink2 >/dev/null 2>&1; then
    echo "'plink2' not found — loading modules..."
    module load StdEnv/2023
    module load plink/2.00a5.8

    if ! command -v plink2 >/dev/null 2>&1; then
        echo "Error: plink2 is still not available after loading modules." >&2
        exit 1
    fi
else
    echo "plink2 is already available."
fi

# === Get prefix ===
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


echo "💻 Running with $cpus cores"
echo "🧠 PLINK memory: $plink_mem MB"

# === Run PLINK ===
plink2 \
  --vcf "$inputVcf" \
  --threads "$cpus" \
  --memory "$plink_mem" \
  --chr 1-22 \
  --snps-only \
  --max-alleles 2 \
  --vcf-min-dp 1 \
  --vcf-half-call "m" \
  --make-bed \
  --out "$outputDir/$outPrefix"  # <--- NE PAS mettre $outputDir ici

# === (Optionnel) copier dans outputDir pour archivage ===
# mv "${outPrefix}".{bed,bim,fam,log} "$outputDir/"
