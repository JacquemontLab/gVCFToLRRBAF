#!/usr/bin/env bash
set -euo pipefail

if (( $# < 3 )); then
  echo "Usage: $0 <iid_list.txt> <BAF_DIR> <n_batches> [<outdir=.>]" >&2
  exit 1
fi

iid_list=$1          # ex: keep_unrel_topN.iid.txt
baf_dir=$2           # ex: /.../LRR_BAF
n_batches=$3         # ex: 100
outdir=${4:-.}       # ex: .

mkdir -p "$outdir/batches"

# 1) Nettoyage: retire un éventuel header 'IID' & garde la 1ère colonne
grep -Ev '^[[:space:]]*IID[[:space:]]*$' "$iid_list" \
  | awk 'BEGIN{FS="[[:space:]]+"} {print $1}' \
  > "$outdir/iid_clean.txt"

# 2) IDs réellement présents dans BAF_DIR (supporte *_filtered.baf_lrr.tsv)
#find "$baf_dir" -maxdepth 1 -type f -name '*.baf_lrr.tsv' -printf '%f\n' \
find -L "$baf_dir" -maxdepth 1 -type f -name '*.baf_lrr.tsv' -printf '%f\n' \
  | sed -E 's/(_filtered)?\.baf_lrr\.tsv$//' \
  | sort -u > "$outdir/have_ids.txt"

# 3) Intersection: garder uniquement les IIDs qui existent vraiment dans BAF_DIR
awk 'NR==FNR{h[$1]=1; next} ($1 in h){print $1}' \
  "$outdir/have_ids.txt" "$outdir/iid_clean.txt" \
  > "$outdir/iid_present.txt"

N=$(wc -l < "$outdir/iid_present.txt" || echo 0)
if (( N == 0 )); then
  echo "No IID present in BAF dir: $baf_dir" >&2
  exit 1
fi

# 4) Split équilibré en n_batches fichiers (iid_batch_000, 001, ...)
split -n l/"$n_batches" -d -a 3 \
  "$outdir/iid_present.txt" "$outdir/batches/iid_batch_"

echo "[make_iid_batches] IIDs totaux: $(wc -l < "$outdir/iid_clean.txt") | présents: $N"
echo "[make_iid_batches] Batches écrits dans: $outdir/batches"
