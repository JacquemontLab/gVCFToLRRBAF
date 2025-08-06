#!/bin/bash

# Environment setup to avoid locale prel in batch environments
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
    
#======================================================================
# Batch parallel CNV caller with cpu managment. Relies on cnv_caling.sh 
#======================================================================
usage() {
  echo "Usage: $0 --batch_list FILE --output_dir PATH"
  echo ""
  echo "Required options:"
  echo "  --batch_list   FILE    TSV file with two columns per line:"
  echo "                         1) Sample ID"
  echo "                         2) Path to the gVCF file"
  echo "  --output_dir   PATH    Path to the output directory"
  echo "  --help                 Show this help message and exit"
  exit 1
}
# Get the absolute path of the script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

#Parse Options
while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch_list)
      batch_list="$2"
      shift 2
      ;;
    --output_dir)
      output_dir="$2"
      shift 2
      ;;
    --help)
      usage
      ;;
    *)
      echo "❌ Unknown option: $1"
      usage
      ;;
  esac
done


# Cancel if required arguments are missing
if [[ -z "$batch_list" || -z "$output_dir" ]]; then
  echo "❌ Error: --batch_list and --output_dir are required."
  usage
  exit 1
fi


# Load bcftools if needed
if ! command -v bcftools >/dev/null 2>&1; then
    echo "'bcftools' not found — loading modules..."
    module load bcftools

    # Check again
    if ! command -v bcftools >/dev/null 2>&1; then
        echo "Error: bcftools is still not available after loading modules." >&2
        exit 1
    fi
else
    echo "bcftools is already available."
fi

# Get the absolute path of the script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

mkdir -p "$output_dir"

# Determine number of CPUs
cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running with $cpus cores"

# Export necessary vars for GNU Parallel
export SCRIPT_DIR
export output_dir

# Run in parallel
cat "$batch_list" | parallel -j "$cpus" --colsep '\t' --eta --line-buffer '
  sample={1}
  gvcf={2}
  echo "🔄 Processing $sample"
  perl "$SCRIPT_DIR/fromgVCFToSignalIntensity.pl" "$gvcf" \
    --sample_id "$sample" \
    --output_dir "$output_dir" \
    --genome_version "GRCh38"
'