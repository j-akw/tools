#!/usr/bin/env python3
"""
Mask a VCF so that specified individuals retain genotypes ONLY at sites where
both haplotypes were called as a target ancestry (default: AFT) with high
confidence. Everywhere else, the genotype is set to './.' for those individuals.

Inputs:
  - A VCF (plain or .gz)
  - Per-individual BED files (output of gnomix_to_bed.py), with the NAME
    column formatted as '<ancestry>_hap{1,2}'
  - One or more sample IDs to mask, each with the path to their BED file

Logic (strict):
  For each target sample, compute per-chromosome intersection of
  ANCESTRY_hap1 intervals AND ANCESTRY_hap2 intervals. A VCF site is kept
  for that sample iff its position falls inside the intersection. Otherwise
  the sample's GT is replaced with './.' (other FORMAT fields untouched).

Other samples in the VCF are passed through unchanged.
"""

import argparse
import gzip
import shutil
import subprocess
import sys
from bisect import bisect_right
from collections import defaultdict


# ---------- BED parsing ----------

def load_bed_ancestry_intervals(bed_path, ancestry):
    """
    Returns dict[chrom] -> {1: [(start, end), ...], 2: [(start, end), ...]}
    Coordinates kept in BED convention (0-based, half-open).
    Only rows whose NAME == f'{ancestry}_hap1' or f'{ancestry}_hap2' are loaded.
    """
    want = {f"{ancestry}_hap1": 1, f"{ancestry}_hap2": 2}
    by_chrom = defaultdict(lambda: {1: [], 2: []})
    with open(bed_path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith('#'):
                continue
            parts = ln.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            name = parts[3]
            if name not in want:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            by_chrom[chrom][want[name]].append((start, end))
    # Sort + merge overlapping/adjacent intervals within each hap
    for chrom in by_chrom:
        for hap in (1, 2):
            by_chrom[chrom][hap] = merge_intervals(sorted(by_chrom[chrom][hap]))
    return by_chrom


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def intersect_intervals(a, b):
    """Both lists sorted, non-overlapping. Returns list of overlap intervals."""
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def build_keep_regions(bed_intervals_by_chrom):
    """
    Given dict[chrom]->{1:[...], 2:[...]}, return dict[chrom]->(starts, ends)
    where starts/ends are sorted parallel lists of the hap1 ∩ hap2 intersection,
    ready for binary-search lookup.
    """
    keep = {}
    for chrom, haps in bed_intervals_by_chrom.items():
        inter = intersect_intervals(haps[1], haps[2])
        starts = [s for s, _ in inter]
        ends = [e for _, e in inter]
        keep[chrom] = (starts, ends)
    return keep


def position_in_keep(keep_for_chrom, pos_1based):
    """
    Returns True iff the VCF position falls in any kept interval.
    Intervals are 0-based half-open [start, end); VCF pos is 1-based,
    so we check if (pos-1) is within some [start, end).
    """
    if keep_for_chrom is None:
        return False
    starts, ends = keep_for_chrom
    if not starts:
        return False
    p0 = pos_1based - 1
    # Largest start <= p0
    idx = bisect_right(starts, p0) - 1
    if idx < 0:
        return False
    return p0 < ends[idx]


# ---------- VCF I/O ----------

def open_vcf_in(path):
    """Read-side: plain gzip is fine for input (handles both .gz and BGZF)."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


class BgzipWriter:
    """Write-side: pipes to `bgzip -c` so tabix can index the result."""
    def __init__(self, path):
        if shutil.which('bgzip') is None:
            sys.exit(
                "ERROR: --out ends in .gz but 'bgzip' is not on PATH. "
                "Install htslib, or write to a plain .vcf and bgzip it yourself."
            )
        self._fout = open(path, 'wb')
        self._proc = subprocess.Popen(
            ['bgzip', '-c'], stdin=subprocess.PIPE, stdout=self._fout
        )

    def write(self, s):
        self._proc.stdin.write(s.encode())

    def close(self):
        self._proc.stdin.close()
        rc = self._proc.wait()
        self._fout.close()
        if rc != 0:
            sys.exit(f"ERROR: bgzip exited with code {rc}")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


def open_vcf_out(path):
    if path.endswith('.gz'):
        return BgzipWriter(path)
    return open(path, 'wt')


def process_vcf(vcf_in, vcf_out, sample_keep_regions):
    """
    sample_keep_regions: dict[sample_id] -> dict[chrom] -> (starts, ends)
    Streams the VCF and rewrites GT for the listed samples.
    """
    sample_col_for = {}    # sample_id -> column index in the data line
    n_records = 0
    n_masked_calls = defaultdict(int)   # per sample, count of sites masked
    n_kept_calls = defaultdict(int)

    with open_vcf_in(vcf_in) as fin, open_vcf_out(vcf_out) as fout:
        for line in fin:
            if line.startswith('##'):
                fout.write(line)
                continue
            if line.startswith('#CHROM'):
                fout.write(line)
                header = line.rstrip('\n').split('\t')
                # FORMAT is column 8 (0-indexed); samples start at column 9
                for i in range(9, len(header)):
                    if header[i] in sample_keep_regions:
                        sample_col_for[header[i]] = i
                missing = set(sample_keep_regions) - set(sample_col_for)
                if missing:
                    sys.exit(f"ERROR: samples not found in VCF: {sorted(missing)}")
                sys.stderr.write(
                    f"Masking samples: {list(sample_col_for.keys())} "
                    f"(VCF columns {list(sample_col_for.values())})\n"
                )
                continue
            if not line.strip():
                fout.write(line)
                continue

            # Data line
            fields = line.rstrip('\n').split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            fmt_keys = fields[8].split(':')
            try:
                gt_idx = fmt_keys.index('GT')
            except ValueError:
                # No GT field — write as-is
                fout.write(line)
                continue

            for sample_id, col in sample_col_for.items():
                keep_for_chrom = sample_keep_regions[sample_id].get(chrom)
                if position_in_keep(keep_for_chrom, pos):
                    n_kept_calls[sample_id] += 1
                    continue
                # Mask GT to ./.
                sample_fields = fields[col].split(':')
                if gt_idx < len(sample_fields):
                    sample_fields[gt_idx] = './.'
                    fields[col] = ':'.join(sample_fields)
                    n_masked_calls[sample_id] += 1

            fout.write('\t'.join(fields) + '\n')
            n_records += 1

    sys.stderr.write(f"Processed {n_records} VCF records.\n")
    for s in sample_col_for:
        sys.stderr.write(
            f"  {s}: kept {n_kept_calls[s]}, masked {n_masked_calls[s]}\n"
        )


# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--vcf', required=True, help='Input VCF (.vcf or .vcf.gz)')
    ap.add_argument('--out', required=True, help='Output VCF path (use .gz suffix for gzipped output)')
    ap.add_argument('--ancestry', default='AFT',
                    help='Ancestry label to retain (must match the NAME prefix in the BEDs). Default: AFT')
    ap.add_argument(
        '--sample',
        action='append',
        required=True,
        metavar='SAMPLE_ID=BED_PATH',
        help='Sample ID and path to its BED file, e.g. --sample BAR4050=beds/BAR4050.bed. '
             'Repeat the flag for multiple samples.',
    )
    args = ap.parse_args()

    sample_keep_regions = {}
    for spec in args.sample:
        if '=' not in spec:
            sys.exit(f"--sample expects SAMPLE_ID=BED_PATH, got: {spec}")
        sid, bed_path = spec.split('=', 1)
        bed_intervals = load_bed_ancestry_intervals(bed_path, args.ancestry)
        keep = build_keep_regions(bed_intervals)
        n_chr = len(keep)
        n_intervals = sum(len(v[0]) for v in keep.values())
        sys.stderr.write(
            f"Sample {sid}: {n_intervals} {args.ancestry}-hap1 ∩ {args.ancestry}-hap2 "
            f"intervals across {n_chr} chromosome(s)\n"
        )
        sample_keep_regions[sid] = keep

    process_vcf(args.vcf, args.out, sample_keep_regions)


if __name__ == '__main__':
    main()
