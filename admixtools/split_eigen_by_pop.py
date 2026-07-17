#!/usr/bin/env python3
"""
split_eigen_by_pop.py — build one EIGENSTRAT dataset per population, each
containing a specific subset of individuals from that pop plus a fixed set
of reference pops.

Takes a flat list of sample IDs, looks each one up in the .ind file to discover its pop, groups
them by pop, and writes one dataset per group.

For every pop P that has one or more listed individuals, creates:
    <outdir>/<P>/<P>.geno
    <outdir>/<P>/<P>.snp
    <outdir>/<P>/<P>.ind

The .ind of each output contains the listed individuals from P first, then
all individuals from the reference pops.

Usage:
    python split_by_individuals.py \\
        -g data.geno -s data.snp -i data.ind \\
        -l selected_individuals.txt \\
        -r ref_pops.txt \\
        -o selected_datasets/

Requires numpy.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

try:
    import numpy as np
except ImportError:
    sys.exit("error: numpy is required. Install with: pip install numpy")


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_id_list(path: Path) -> List[str]:
    ids: List[str] = []
    seen: Set[str] = set()
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split()[0]
            if tok not in seen:
                ids.append(tok)
                seen.add(tok)
    return ids


def parse_pop_list(file_arg: Path | None, inline_arg: str | None) -> List[str]:
    pops: List[str] = []
    seen: Set[str] = set()
    if file_arg:
        with file_arg.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                tok = line.split()[0]
                if tok not in seen:
                    pops.append(tok)
                    seen.add(tok)
    if inline_arg:
        for tok in inline_arg.split(","):
            tok = tok.strip()
            if tok and tok not in seen:
                pops.append(tok)
                seen.add(tok)
    return pops


def parse_ind(path: Path) -> List[Tuple[str, str, str]]:
    rows: List[Tuple[str, str, str]] = []
    with path.open() as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 3:
                rows.append((parts[0], parts[1], parts[2]))
    return rows


def write_ind_fixed(fh, sid: str, sex: str, pop: str) -> None:
    fh.write(f"{sid:>20} {sex:1} {pop:>10}\n")


# ---------------------------------------------------------------------------
# Geno memmap helpers
# ---------------------------------------------------------------------------

def inspect_geno(path: Path, n_samples_expected: int) -> Tuple[int, int]:
    """
    Peek at the .geno file to determine (n_snps, row_bytes).
    row_bytes includes the newline. Validates that the file is a clean
    rectangle and that the row width matches the .ind file.
    """
    file_size = path.stat().st_size
    with path.open("rb") as fh:
        first = fh.readline()
    if not first:
        sys.exit(f"error: {path} is empty")
    row_bytes = len(first)
    # Detect CRLF — we don't support it here (would need a conversion pass).
    if first.endswith(b"\r\n"):
        sys.exit(f"error: {path} has CRLF line endings; convert with "
                 f"`dos2unix {path}` first")
    if not first.endswith(b"\n"):
        sys.exit(f"error: {path} first line is not newline-terminated")

    data_chars = row_bytes - 1  # strip newline
    if data_chars != n_samples_expected:
        sys.exit(
            f"error: .geno row width is {data_chars} chars but .ind has "
            f"{n_samples_expected} samples"
        )
    if file_size % row_bytes != 0:
        sys.exit(
            f"error: {path} size ({file_size}) is not a multiple of row size "
            f"({row_bytes}); file is not a clean rectangle"
        )
    n_snps = file_size // row_bytes
    return n_snps, row_bytes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-g", "--geno", type=Path, required=True, help="input .geno")
    p.add_argument("-s", "--snp",  type=Path, required=True, help="input .snp")
    p.add_argument("-i", "--ind",  type=Path, required=True, help="input .ind")
    p.add_argument("-o", "--outdir", type=Path, required=True,
                   help="output directory (one subdir per pop with listed samples)")

    p.add_argument("-l", "--list", type=Path, required=True,
                   help="file of sample IDs to include (one per line)")

    p.add_argument("-r", "--refs-file", type=Path,
                   help="file of reference pop names, one per line")
    p.add_argument("--refs", type=str,
                   help="comma-separated ref pops (alternative to -r)")

    p.add_argument("--snp-mode", choices=["copy", "symlink", "hardlink"],
                   default="symlink",
                   help="how to materialise the per-target .snp file (default: symlink)")
    p.add_argument("--dry-run", action="store_true",
                   help="report plan without writing anything")

    args = p.parse_args()

    # --- Load inputs ---
    wanted_ids = parse_id_list(args.list)
    if not wanted_ids:
        sys.exit(f"error: no sample IDs found in {args.list}")
    print(f"loaded {len(wanted_ids)} sample IDs from {args.list}", file=sys.stderr)

    refs = parse_pop_list(args.refs_file, args.refs)
    if not refs:
        print("warning: no reference pops given; outputs will contain listed individuals only",
              file=sys.stderr)
    ref_set = set(refs)

    ind_rows = parse_ind(args.ind)
    n_samples_total = len(ind_rows)

    id_to_col: Dict[str, int] = {}
    id_to_pop: Dict[str, str] = {}
    for col, (sid, _sex, pop) in enumerate(ind_rows):
        if sid not in id_to_col:
            id_to_col[sid] = col
            id_to_pop[sid] = pop
    print(f"loaded {n_samples_total} samples from {args.ind}", file=sys.stderr)

    # --- Group wanted IDs by pop ---
    missing: List[str] = []
    by_pop: Dict[str, List[int]] = defaultdict(list)
    for sid in wanted_ids:
        col = id_to_col.get(sid)
        if col is None:
            missing.append(sid)
            continue
        pop = id_to_pop[sid]
        if pop in ref_set:
            print(f"warning: listed individual {sid!r} belongs to reference pop "
                  f"{pop!r}; including anyway", file=sys.stderr)
        by_pop[pop].append(col)

    if missing:
        print(f"warning: {len(missing)} listed IDs not found in .ind "
              f"(first 10: {missing[:10]})", file=sys.stderr)

    if not by_pop:
        sys.exit("error: none of the listed IDs matched samples in the .ind file")

    ref_cols: List[int] = [
        i for i, (_, _, pop) in enumerate(ind_rows) if pop in ref_set
    ]
    print(f"reference pops resolve to {len(ref_cols)} samples", file=sys.stderr)

    # Plan per-pop output: listed individuals first, then refs (deduped)
    plans: Dict[str, List[int]] = {}
    for pop in sorted(by_pop):
        listed = by_pop[pop]
        seen = set(listed)
        merged = listed + [c for c in ref_cols if c not in seen]
        plans[pop] = merged
        print(f"  {pop}: {len(listed)} listed + "
              f"{len(merged) - len(listed)} ref = {len(merged)} total",
              file=sys.stderr)

    if args.dry_run:
        print("(dry run — no files written)", file=sys.stderr)
        return 0

    # --- Inspect and memmap the input .geno ---
    n_snps, row_bytes = inspect_geno(args.geno, n_samples_total)
    print(f"geno: {n_snps} SNPs × {n_samples_total} samples "
          f"({args.geno.stat().st_size / 1e9:.2f} GB)", file=sys.stderr)

    # Shape: (n_snps, row_bytes) where the last column is '\n'
    arr = np.memmap(args.geno, dtype=np.uint8, mode="r",
                    shape=(n_snps, row_bytes))

    # --- Write per-pop outputs ---
    args.outdir.mkdir(parents=True, exist_ok=True)
    newline = np.uint8(ord("\n"))

    for pop, cols in plans.items():
        sub = args.outdir / pop
        sub.mkdir(parents=True, exist_ok=True)
        out_root = sub / pop

        # .ind
        with (out_root.with_suffix(".ind")).open("w") as ind_out:
            for c in cols:
                sid, sex, p = ind_rows[c]
                write_ind_fixed(ind_out, sid, sex, p)

        # .snp (shared across all outputs — copy/link from source)
        snp_dst = out_root.with_suffix(".snp")
        if snp_dst.exists() or snp_dst.is_symlink():
            snp_dst.unlink()
        src_abs = args.snp.resolve()
        if args.snp_mode == "copy":
            shutil.copy2(src_abs, snp_dst)
        elif args.snp_mode == "hardlink":
            os.link(src_abs, snp_dst)
        else:
            os.symlink(src_abs, snp_dst)

        # .geno — vectorized column extraction
        col_idx = np.asarray(cols, dtype=np.int64)
        n_out_cols = len(cols)
        # Build output array with trailing newline column
        out = np.empty((n_snps, n_out_cols + 1), dtype=np.uint8)
        # Fancy-index the memmap. This materializes one pop's worth of data
        # in RAM: n_snps * n_out_cols bytes. For 5M SNPs × 100 samples = 500 MB.
        out[:, :n_out_cols] = arr[:, col_idx]
        out[:, n_out_cols] = newline
        out.tofile(str(out_root.with_suffix(".geno")))
        print(f"  wrote {pop}: {n_snps} SNPs × {n_out_cols} samples "
              f"({out.nbytes / 1e6:.1f} MB)", file=sys.stderr)

    print(f"done. wrote {len(plans)} pop datasets", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
