#!/usr/bin/env python3
"""
pfb_partial.py  —  calcule sum_baf + n par (Name, Chr, Position) pour un batch d'IIDs.
Traite les fichiers un à un pour éviter les problèmes mémoire et le streaming Polars.
"""
import os, sys, glob
import polars as pl


def read_iids(path):
    ids = []
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            if not line.strip():
                continue
            tok = line.split()[0]
            if ln == 1 and tok.upper() == "IID":
                continue
            ids.append(tok)
    return ids


def detect_baf_col(tsv, iid):
    with open(tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    for c in (f"{iid}.B Allele Freq", "B Allele Freq", "BAF"):
        if c in hdr:
            return c
    return None


def read_and_agg_one(tsv, baf_col):
    """Lit un fichier BAF et retourne les stats agrégées par position (sum_baf, n)."""
    return (
        pl.read_csv(
            tsv,
            separator="\t",
            has_header=True,
            schema_overrides={"Chr": pl.Utf8, "Position": pl.Int64, baf_col: pl.Float64},
            ignore_errors=True,
        )
        .select(["Name", "Chr", "Position", pl.col(baf_col).alias("BAF")])
        .drop_nulls(["BAF"])
        .group_by(["Name", "Chr", "Position"])
        .agg(
            pl.col("BAF").sum().alias("sum_baf"),
            pl.len().alias("n"),
        )
    )


def main():
    if len(sys.argv) < 4:
        sys.stderr.write(
            "Usage: pfb_partial.py <iid_batch.txt> <baf_dir> <out_partial.parquet>\n"
        )
        sys.exit(1)

    iid_batch, baf_dir, out_pq = sys.argv[1:4]
    iids = read_iids(iid_batch)

    missing = 0
    nobaf = 0
    acc_df = None

    for i, iid in enumerate(iids):
        # Cherche le fichier BAF (nom exact ou suffixe _filtered)
        candidates = [
            os.path.join(baf_dir, f"{iid}.baf_lrr.tsv"),
            os.path.join(baf_dir, f"{iid}_filtered.baf_lrr.tsv"),
        ]
        tsv = next((p for p in candidates if os.path.exists(p)), None)
        if tsv is None:
            hits = glob.glob(os.path.join(baf_dir, f"{iid}*.baf_lrr.tsv"))
            tsv = hits[0] if hits else None
        if tsv is None:
            missing += 1
            continue

        baf_col = detect_baf_col(tsv, iid)
        if baf_col is None:
            nobaf += 1
            continue

        partial = read_and_agg_one(tsv, baf_col)

        if acc_df is None:
            acc_df = partial
        else:
            # Fusionne avec l'accumulateur : re-agrège sum_baf et n
            acc_df = (
                pl.concat([acc_df, partial], how="vertical")
                .group_by(["Name", "Chr", "Position"])
                .agg(
                    pl.col("sum_baf").sum(),
                    pl.col("n").sum(),
                )
            )

        if (i + 1) % 10 == 0:
            sys.stderr.write(f"[pfb_partial] {i+1}/{len(iids)} fichiers traités\n")

    if acc_df is None:
        raise SystemExit("No usable files in this batch.")

    sys.stderr.write(
        f"[pfb_partial] terminé : {len(iids)-missing-nobaf} ok, "
        f"{missing} manquants, {nobaf} sans col BAF\n"
    )

    tmp = out_pq + ".tmp"
    acc_df.write_parquet(tmp, compression="zstd")
    os.replace(tmp, out_pq)


if __name__ == "__main__":
    main()
