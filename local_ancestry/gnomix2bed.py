#!/usr/bin/env python3
"""
Convert GNOMIX .fb + .msp output into per-individual BED files of
high-probability ancestry regions (posterior > threshold, default 0.99).

For each individual, contiguous windows where a given (haplotype, ancestry)
exceeds the threshold are merged into a single BED interval. The score
column holds the minimum probability across the merged interval.

Output BED columns:
    chrom   start(0-based)   end   name(ANC_hap{1,2})   score(min posterior)
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_fb(fb_path):
    """
    Returns:
        ancestries: list of population codes, e.g. ['EUT','AFT','ASI'] (index = ancestry id)
        sample_hap_anc_cols: dict[(sample, hap, anc_idx)] -> column index
        samples: ordered list of unique sample IDs
        rows: list of [chrom, ppos, gpos, marker, *float_probs] (one per window)
    """
    with open(fb_path) as f:
        lines = f.readlines()

    # Line 0: "#reference_panel_population:\tEUT\tAFT\tASI"
    pop_line = lines[0].rstrip('\n').split('\t')
    ancestries = pop_line[1:]
    anc_idx = {a: i for i, a in enumerate(ancestries)}

    # Line 1: header. First 4 cols are chrom/ppos/gpos/marker, then SAMPLE:::hapN:::ANC
    header = lines[1].rstrip('\n').split('\t')
    sample_order = []
    seen = set()
    sample_hap_anc_cols = {}
    for col_i, col in enumerate(header[4:], start=4):
        parts = col.split(':::')
        if len(parts) != 3:
            continue
        sample, hap, anc = parts
        if sample not in seen:
            seen.add(sample)
            sample_order.append(sample)
        hap_num = 1 if hap == 'hap1' else 2
        sample_hap_anc_cols[(sample, hap_num, anc_idx[anc])] = col_i

    # Data rows
    rows = []
    for ln in lines[2:]:
        if not ln.strip():
            continue
        parts = ln.rstrip('\n').split('\t')
        rows.append(parts)

    return ancestries, sample_hap_anc_cols, sample_order, rows


def parse_msp_windows(msp_path):
    """
    Returns list of (chrom, spos, epos) tuples in window order.
    """
    windows = []
    with open(msp_path) as f:
        lines = f.readlines()
    # Line 0: #Subpopulation order
    # Line 1: header
    for ln in lines[2:]:
        if not ln.strip():
            continue
        parts = ln.rstrip('\n').split('\t')
        chrom = parts[0]
        spos = int(parts[1])
        epos = int(parts[2])
        windows.append((chrom, spos, epos))
    return windows


def extract_regions(fb_rows, windows, sample_hap_anc_cols, ancestries,
                    samples, threshold):
    """
    For each (sample, hap, anc), walk windows and emit merged BED intervals
    where posterior > threshold.

    Returns: dict[sample] -> list of (chrom, start_0based, end, name, min_prob)
    """
    assert len(fb_rows) == len(windows), \
        f".fb has {len(fb_rows)} rows but .msp has {len(windows)} windows"

    per_sample = defaultdict(list)

    for sample in samples:
        for hap in (1, 2):
            for a_idx, a_name in enumerate(ancestries):
                col = sample_hap_anc_cols[(sample, hap, a_idx)]
                run_start = None     # spos of first window in current run
                run_end = None       # epos of last window in current run
                run_chrom = None
                run_min_prob = None

                for win_i, (chrom, spos, epos) in enumerate(windows):
                    prob = float(fb_rows[win_i][col])
                    if prob > threshold:
                        if run_start is None:
                            run_start = spos
                            run_end = epos
                            run_chrom = chrom
                            run_min_prob = prob
                        else:
                            # contiguous if same chromosome (always true within one .fb file)
                            run_end = epos
                            if prob < run_min_prob:
                                run_min_prob = prob
                    else:
                        if run_start is not None:
                            per_sample[sample].append((
                                run_chrom,
                                run_start - 1,   # BED is 0-based
                                run_end,
                                f"{a_name}_hap{hap}",
                                run_min_prob,
                            ))
                            run_start = run_end = run_chrom = run_min_prob = None

                # Flush trailing run
                if run_start is not None:
                    per_sample[sample].append((
                        run_chrom,
                        run_start - 1,
                        run_end,
                        f"{a_name}_hap{hap}",
                        run_min_prob,
                    ))

    return per_sample


def write_beds(per_sample, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for sample, intervals in per_sample.items():
        # Sort by chrom, then start, then name (stable, predictable output)
        intervals.sort(key=lambda x: (x[0], x[1], x[3]))
        out_path = os.path.join(out_dir, f"{sample}.bed")
        with open(out_path, 'w') as out:
            for chrom, start, end, name, score in intervals:
                out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score:.6f}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fb', required=True, help='GNOMIX .fb file')
    ap.add_argument('--msp', required=True, help='GNOMIX .msp file')
    ap.add_argument('--out-dir', required=True, help='Output directory for per-sample BED files')
    ap.add_argument('--threshold', type=float, default=0.99,
                    help='Posterior probability threshold (strictly greater than). Default: 0.99')
    args = ap.parse_args()

    ancestries, cols, samples, fb_rows = parse_fb(args.fb)
    windows = parse_msp_windows(args.msp)

    sys.stderr.write(
        f"Parsed {len(fb_rows)} windows, {len(samples)} individuals, "
        f"ancestries: {ancestries}\n"
    )

    per_sample = extract_regions(fb_rows, windows, cols, ancestries,
                                 samples, args.threshold)

    write_beds(per_sample, args.out_dir)

    total = sum(len(v) for v in per_sample.values())
    sys.stderr.write(
        f"Wrote {len(per_sample)} BED files to {args.out_dir} "
        f"({total} intervals total at threshold > {args.threshold})\n"
    )


if __name__ == '__main__':
    main()
