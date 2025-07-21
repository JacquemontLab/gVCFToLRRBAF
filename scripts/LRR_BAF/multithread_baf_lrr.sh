#!/bin/bash
    
#======================================================================
# Batch parallel CNV caller with cpu managment. Relies on cnv_caling.sh 
#======================================================================
usage() {
  echo "Usage: $0 --batch_list FILE"
  echo ""
  echo "Required options:"
  echo "  --batch_list   FILE    Text file with one file path per line"
  echo "  --output_dir   PATH    Path to the output directory"
  echo "  --help                 Show this help message and exit"
  exit 1
}

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

# Get the absolute path of the script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

mkdir -p "$output_dir"

# Determine number of CPUs
cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running with $cpus cores"

# Export necessary vars for GNU Parallel
export SCRIPT_DIR output_dir

# Run in parallel
cat "$batch_list" | parallel -j "$cpus" --colsep '\t' --eta --line-buffer '
  sample={1}
  gvcf={2}
  echo "🔄 Processing $sample"
  perl "$SCRIPT_DIR/fromgVCFToSignalIntensity.pl" "$gvcf" \
    --sample_id "$output_dir/$sample" \
    --genome_version "GRCh38"
'