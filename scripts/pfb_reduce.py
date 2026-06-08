#!/usr/bin/env python3
import sys, os, glob
import polars as pl

def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: pfb_reduce.py <out.tsv> <partial1.parquet> [partial2.parquet ...] | <dir>\n")
        sys.exit(1)

    out_tsv = sys.argv[1]
    inputs  = sys.argv[2:]

    # Autoriser d'appeler avec un dossier contenant des .parquet
    if len(inputs) == 1 and os.path.isdir(inputs[0]):
        inputs = sorted(glob.glob(os.path.join(inputs[0], "*.parquet")))

    if not inputs:
        sys.stderr.write("No partial parquet provided.\n")
        sys.exit(1)

    # Lazy scan → agrégation → collecte SANS streaming=True
    lf = pl.scan_parquet(inputs).select("Name", "Chr", "Position", "sum_baf", "n")

    df = (
        lf.group_by(["Name", "Chr", "Position"])
          .agg(
              pl.col("sum_baf").sum().alias("SUM"),
              pl.col("n").sum().alias("N"),
          )
          .with_columns((pl.col("SUM") / pl.col("N")).alias("PFB"))
          .select("Name", "Chr", "Position", pl.col("PFB").round(3))
          .sort("Name")
          .collect()                      # <-- plus de streaming=True
    )

    df.write_csv(out_tsv, separator="\t")
    print(f"[pfb_reduce] Wrote {out_tsv}  (rows={df.height})")

if __name__ == "__main__":
    main()
