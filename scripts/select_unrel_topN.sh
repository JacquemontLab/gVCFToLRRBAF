#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <kin0_file> <smiss_file> <topN> <out_iid_file>" >&2
  exit 2
fi

kin0="$1"          # stats/merged_dataset_unrelated.kin0
smiss="$2"         # stats/merged_dataset.missing.smiss  (#FID IID MISSING_CT OBS_CT F_MISS)
topN="$3"          # e.g. 1000
out="$4"           # e.g. keep_unrel_topN.iid.txt

KIN_THRESH="${REL_KIN_CUTOFF:-0.0442}"   # >= 3e degré
tmpdir="$(mktemp -d)"; trap 'rm -rf "$tmpdir"' EXIT

# IID -> F_MISS
awk 'NR>1{print $2, $5}' "$smiss" > "$tmpdir/iid_fmiss.tsv"   # IID F_MISS

# IIDs à retirer (dans chaque paire KINSHIP>=seuil, on retire le pire call rate)
awk -v kin="$KIN_THRESH" '
  NR==FNR { fmiss[$1]=$2; next }
  $1!~/^#/ && $8>=kin {
    i1=$2; i2=$4;
    if ( (i1 in fmiss) && (i2 in fmiss) ) {
      if      (fmiss[i1]>fmiss[i2]) print i1;
      else if (fmiss[i2]>fmiss[i1]) print i2;
      else                          print (i1<i2? i2 : i1);
    }
  }
' "$tmpdir/iid_fmiss.tsv" "$kin0" | sort -u > "$tmpdir/to_remove.txt"

# Garde non-apparentés, trie par meilleur call rate (F_MISS croissant), prend topN → IID only
awk 'NR==FNR{rm[$1]=1; next} NR>1 && !($2 in rm){ print $2, $5 }' \
    "$tmpdir/to_remove.txt" "$smiss" > "$tmpdir/candidates.tsv"

# 4) Trier par meilleur call rate (F_MISS croissant), tie-break par IID, puis prendre topN
LC_ALL=C sort -k2,2n -k1,1 "$tmpdir/candidates.tsv" > "$tmpdir/candidates.sorted.tsv"
head -"$topN" "$tmpdir/candidates.sorted.tsv" | awk '{print $1}' > "$out"

# 5) Petit log
echo "[KIN]  threshold  = $KIN_THRESH"
echo "[KEEP] selected N = $(wc -l < "$out")"
echo "[OUT]  -> $out"
