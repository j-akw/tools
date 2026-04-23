#!/usr/bin/env python3
"""
Convert PLINK --freq --within output (.frq.strat) to TreeMix input format.
"""
import argparse
import gzip
import sys
import pandas as pd


def smart_open(path, mode):
    if path == "-":
        return sys.stdout if "w" in mode else sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("infile", help="PLINK .frq.strat file (.gz OK)")
    ap.add_argument("-o", "--out", default="-",
                    help="Output file; use '-' for stdout (default). "
                         ".gz extension triggers gzip output.")
    ap.add_argument("--keep-missing", action="store_true",
                    help="Emit '0,0' for populations missing a SNP instead of "
                         "dropping the SNP entirely.")
    ap.add_argument("--sort-pops", action="store_true",
                    help="Alphabetize populations in the header (default: "
                         "order of first appearance in input).")
    args = ap.parse_args()

    df = pd.read_csv(args.infile, sep=r"\s+")

    required = {"SNP", "CLST", "MAC", "NCHROBS"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: input missing required columns: {sorted(missing)}")

    pops = sorted(df["CLST"].unique()) if args.sort_pops \
        else list(pd.unique(df["CLST"]))
    snp_order = pd.unique(df["SNP"])  # preserve input order

    sys.stderr.write(f"Populations ({len(pops)}): {' '.join(map(str, pops))}\n")
    sys.stderr.write(f"SNPs in input: {len(snp_order)}\n")

    # Build "A2,A1" count strings in one vectorized pass.
    a2 = df["NCHROBS"] - df["MAC"]
    df["counts"] = a2.astype(str) + "," + df["MAC"].astype(str)

    wide = (df.pivot(index="SNP", columns="CLST", values="counts")
              .reindex(index=snp_order, columns=pops))

    if args.keep_missing:
        wide = wide.fillna("0,0")
    else:
        before = len(wide)
        wide = wide.dropna()
        if before - len(wide):
            sys.stderr.write(
                f"Dropped {before - len(wide)} SNPs missing in >=1 population\n"
            )

    sys.stderr.write(f"SNPs written: {len(wide)}\n")

    with smart_open(args.out, "wt") as out:
        out.write(" ".join(map(str, pops)) + "\n")
        # itertuples is much faster than iterrows for this volume
        for row in wide.itertuples(index=False, name=None):
            out.write(" ".join(row) + "\n")


if __name__ == "__main__":
    main()
