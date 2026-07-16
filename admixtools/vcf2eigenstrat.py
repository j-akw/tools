#!/usr/bin/env python3
"""
vcf2eigenstrat.py — convert a VCF into EIGENSTRAT (.snp/.ind/.geno) format.

Fast text-parsing version. No pysam dependency — reads VCF directly via the
standard library (plain or .gz). For region filtering, pre-slice with:
    bcftools view -r chr1:1-1000000 in.vcf.gz -Oz -o slice.vcf.gz

Features:
  * Python 3 stdlib only.
  * Precomputed GT -> EIGENSTRAT lookup table (fast path).
  * Correct indel vs multiallelic counters.
  * --split-multiallelic to keep multiallelic SNPs as multiple biallelic rows.
  * Optional PLINK genetic map (--genetic-map) for the cM column.
  * Optional synthetic all-ref sample (-r NAME).
  * Fixed-width .ind and .snp output matching convertf convention.

EIGENSTRAT geno encoding (count of REFERENCE alleles):
    2 = hom-ref, 1 = het, 0 = hom-alt, 9 = missing
"""

from __future__ import annotations

import argparse
import gzip
import io
import logging
import sys
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Tuple

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

log = logging.getLogger("vcf2eigenstrat")


class _FieldCountMismatch(Exception):
    def __init__(self, n: int) -> None:
        self.n = n


# ---------------------------------------------------------------------------
# GT -> EIGENSTRAT lookup (fast path, biallelic only)
# ---------------------------------------------------------------------------

GT_FAST: Dict[str, str] = {
    "0/0": "2", "0/1": "1", "1/0": "1", "1/1": "0",
    "./.": "9", "./0": "9", "0/.": "9", "./1": "9", "1/.": "9",
    "0|0": "2", "0|1": "1", "1|0": "1", "1|1": "0",
    ".|.": "9", ".|0": "9", "0|.": "9", ".|1": "9", "1|.": "9",
    "0": "2", "1": "0", ".": "9",
}


def encode_gt_slow(gt: str, alt_idx: int = 1) -> str:
    """Fallback for GTs the fast table doesn't cover (multiallelic, higher ploidy)."""
    parts = gt.replace("|", "/").split("/")
    ref_count = 0
    alt_str = str(alt_idx)
    for a in parts:
        if a == "0":
            ref_count += 1
        elif a == alt_str:
            pass
        else:
            return "9"
    ploidy = len(parts)
    if ploidy == 1:
        return "2" if ref_count == 1 else "0"
    if ploidy == 2:
        return ("0", "1", "2")[ref_count]
    if ref_count == ploidy:
        return "2"
    if ref_count == 0:
        return "0"
    return "1"


def geno_row_numpy(gt_block: str, n_samples: int) -> Optional[str]:
    """
    Vectorized biallelic diploid GT encoder. Returns None if the block isn't
    homogeneous 3-char-per-sample diploid (caller should fall back).

    Expects `gt_block` to be just the GT portions joined with tabs, e.g.
    "0/0\\t0/1\\t1/1\\t./.". Handles both / and | separators.
    """
    # Fast structural check: n samples of 3 chars each + (n-1) tabs
    if len(gt_block) != 4 * n_samples - 1:
        return None

    arr = np.frombuffer(gt_block.encode("ascii"), dtype=np.uint8)
    # first allele at offsets 0, 4, 8, ...; second allele at 2, 6, 10, ...
    a1 = arr[0::4]
    a2 = arr[2::4]
    if a1.shape[0] != n_samples or a2.shape[0] != n_samples:
        return None

    # ASCII: '0'=48, '1'=49, '.'=46
    is_ref_1 = a1 == 48
    is_alt_1 = a1 == 49
    is_ref_2 = a2 == 48
    is_alt_2 = a2 == 49
    valid = (is_ref_1 | is_alt_1) & (is_ref_2 | is_alt_2)

    # Count REF alleles -> ASCII '0'/'1'/'2'
    count = is_ref_1.view(np.uint8) + is_ref_2.view(np.uint8) + 48
    # Missing -> '9' (57)
    np.putmask(count, ~valid, 57)
    return count.tobytes().decode("ascii")


# ---------------------------------------------------------------------------
# Sample metadata
# ---------------------------------------------------------------------------

@dataclass
class SampleMeta:
    sex: str = "U"
    pop: str = "POP"


def load_ind_map(path: Path) -> Dict[str, SampleMeta]:
    meta: Dict[str, SampleMeta] = {}
    with path.open() as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.split()
            if len(parts) < 3:
                if parts:
                    log.warning("ind map line %d malformed, skipping: %r", lineno, line)
                continue
            meta[parts[0]] = SampleMeta(sex=parts[1], pop=parts[2])
    return meta


def resolve_samples(
    samples: Tuple[str, ...],
    ind_map: Optional[Dict[str, SampleMeta]],
    ind_as_pop: bool,
) -> List[SampleMeta]:
    out: List[SampleMeta] = []
    missing = 0
    for s in samples:
        if ind_map is not None:
            if s in ind_map:
                out.append(ind_map[s])
            else:
                missing += 1
                out.append(SampleMeta())
        elif ind_as_pop:
            out.append(SampleMeta(pop=s))
        else:
            out.append(SampleMeta())
    if missing:
        log.warning("%d VCF samples not found in ind map; defaulted to U/POP", missing)
    return out


# ---------------------------------------------------------------------------
# Genetic map
# ---------------------------------------------------------------------------

class GeneticMap:
    """PLINK-format .map: `chrom rsid cM bp`. Linear interpolation in bp."""

    def __init__(self, path: Path) -> None:
        self._by_chrom: Dict[str, List[Tuple[int, float]]] = {}
        with path.open() as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 4:
                    continue
                try:
                    self._by_chrom.setdefault(parts[0], []).append(
                        (int(parts[3]), float(parts[2]))
                    )
                except ValueError:
                    continue
        for v in self._by_chrom.values():
            v.sort()

    def cm(self, chrom: str, bp: int) -> float:
        entries = self._by_chrom.get(chrom)
        if not entries:
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            entries = self._by_chrom.get(alt)
        if not entries:
            return 0.0
        if bp <= entries[0][0]:
            return entries[0][1]
        if bp >= entries[-1][0]:
            return entries[-1][1]
        lo, hi = 0, len(entries) - 1
        while lo < hi - 1:
            mid = (lo + hi) // 2
            if entries[mid][0] <= bp:
                lo = mid
            else:
                hi = mid
        bp0, cm0 = entries[lo]
        bp1, cm1 = entries[hi]
        return cm0 + (cm1 - cm0) * (bp - bp0) / (bp1 - bp0)


# ---------------------------------------------------------------------------
# VCF reader
# ---------------------------------------------------------------------------

def open_vcf(path: Path) -> TextIO:
    if str(path).endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")  # type: ignore
    return path.open("r", encoding="utf-8")


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def convert(
    vcf_path: Path,
    out_root: Path,
    ind_map_path: Optional[Path],
    ind_as_pop: bool,
    ref_sample: Optional[str],
    genetic_map_path: Optional[Path],
    split_multiallelic: bool,
) -> None:
    ind_map = load_ind_map(ind_map_path) if ind_map_path else None
    gmap = GeneticMap(genetic_map_path) if genetic_map_path else None

    gt_fast_get = GT_FAST.get  # local alias for hot loop
    encode_slow = encode_gt_slow

    n_written = 0
    n_indel = 0
    n_multi = 0
    n_noalt = 0
    progress_every = 100_000

    with ExitStack() as stack:
        vcf = stack.enter_context(open_vcf(vcf_path))
        snp_fh = stack.enter_context(open(f"{out_root}.snp", "w", buffering=1 << 20))
        ind_fh = stack.enter_context(open(f"{out_root}.ind", "w"))
        geno_fh = stack.enter_context(open(f"{out_root}.geno", "w", buffering=1 << 20))

        # --- header ---
        samples: Tuple[str, ...] = ()
        for line in vcf:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = tuple(line.rstrip("\n").split("\t")[9:])
                break
        else:
            sys.exit("error: no #CHROM header line found in VCF")

        metas = resolve_samples(samples, ind_map, ind_as_pop)

        def write_ind(sid: str, sex: str, pop: str) -> None:
            ind_fh.write(f"{sid:>20} {sex:1} {pop:>10}\n")

        if ref_sample:
            write_ind(ref_sample, "U", "REF")
        for sid, m in zip(samples, metas):
            write_ind(sid, m.sex, m.pop)
        ind_fh.flush()

        ref_prefix = "2" if ref_sample else ""
        n_samples = len(samples)

        # --- body ---
        for line in vcf:
            if line[0] == "#":
                continue

            # Slice the 9 fixed fields by hand; avoids a full split.
            tab1 = line.find("\t")
            tab2 = line.find("\t", tab1 + 1)
            tab3 = line.find("\t", tab2 + 1)
            tab4 = line.find("\t", tab3 + 1)
            tab5 = line.find("\t", tab4 + 1)
            tab6 = line.find("\t", tab5 + 1)
            tab7 = line.find("\t", tab6 + 1)
            tab8 = line.find("\t", tab7 + 1)
            tab9 = line.find("\t", tab8 + 1)
            if tab9 == -1:
                continue

            chrom = line[:tab1]
            pos_s = line[tab1 + 1 : tab2]
            rsid_in = line[tab2 + 1 : tab3]
            ref = line[tab3 + 1 : tab4]
            alt = line[tab4 + 1 : tab5]
            fmt = line[tab8 + 1 : tab9]

            if len(ref) != 1:
                n_indel += 1
                continue

            if "," in alt:
                if not split_multiallelic:
                    n_multi += 1
                    continue
                alts: Tuple[str, ...] = tuple(alt.split(","))
            else:
                if not alt or alt == ".":
                    n_noalt += 1
                    continue
                alts = (alt,)

            sample_block = line[tab9 + 1 :].rstrip("\r\n")

            # Lazy split: we only pay for the split if we actually need
            # per-sample fields (fmt != "GT" or numpy fast path fails).
            sample_fields: Optional[List[str]] = None

            def get_fields() -> List[str]:
                nonlocal sample_fields
                if sample_fields is None:
                    sample_fields = sample_block.split("\t")
                    if len(sample_fields) != n_samples:
                        raise _FieldCountMismatch(len(sample_fields))
                return sample_fields

            # GT is always the first FORMAT field per VCF spec.
            if fmt == "GT":
                gt_block_for_numpy: Optional[str] = sample_block
                gts: Optional[List[str]] = None  # defer split until needed
            else:
                try:
                    fields = get_fields()
                except _FieldCountMismatch as e:
                    log.warning("line has %d sample columns, expected %d; skipping",
                                e.n, n_samples)
                    continue
                gts = [f[: f.find(":")] if ":" in f else f for f in fields]
                gt_block_for_numpy = None

            for alt_idx, this_alt in enumerate(alts, start=1):
                if len(this_alt) != 1:
                    n_indel += 1
                    continue

                geno_row: Optional[str] = None

                # Fast path: numpy vectorized (biallelic diploid only)
                if HAS_NUMPY and len(alts) == 1:
                    if gt_block_for_numpy is None and gts is not None:
                        gt_block_for_numpy = "\t".join(gts)
                    if gt_block_for_numpy is not None:
                        geno_row = geno_row_numpy(gt_block_for_numpy, n_samples)

                # Fallback: per-sample dict / slow encoder
                if geno_row is None:
                    if gts is None:
                        # fmt was "GT" but numpy fast path bailed (mixed ploidy etc.)
                        try:
                            gts = get_fields()
                        except _FieldCountMismatch as e:
                            log.warning("line has %d sample columns, expected %d; skipping",
                                        e.n, n_samples)
                            break
                    row_chars: List[str] = []
                    append = row_chars.append
                    if len(alts) == 1:
                        for gt in gts:
                            c = gt_fast_get(gt)
                            if c is None:
                                c = encode_slow(gt, 1)
                            append(c)
                    else:
                        for gt in gts:
                            append(encode_slow(gt, alt_idx))
                    geno_row = "".join(row_chars)

                if ref_prefix:
                    geno_fh.write(ref_prefix)
                geno_fh.write(geno_row)
                geno_fh.write("\n")

                if rsid_in and rsid_in != ".":
                    rsid = rsid_in if len(alts) == 1 else f"{rsid_in}_{this_alt}"
                else:
                    rsid = f"{chrom}_{pos_s}_{ref}_{this_alt}"

                pos_i = int(pos_s)
                cm = gmap.cm(chrom, pos_i) if gmap else 0.0
                snp_fh.write(
                    f"{rsid:>20}{chrom:>6}{cm:>16.6f}{pos_i:>16} {ref} {this_alt}\n"
                )
                n_written += 1

                if n_written % progress_every == 0:
                    log.info("wrote %d sites", n_written)

    log.info("done. wrote %d sites", n_written)
    log.info("skipped: %d indel, %d multiallelic, %d no-alt", n_indel, n_multi, n_noalt)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert VCF to EIGENSTRAT format (fast).")
    p.add_argument("-v", "--vcf", required=True, type=Path, help="input VCF (.vcf or .vcf.gz)")
    p.add_argument("-o", "--out", required=True, type=Path, help="output root (writes .snp/.ind/.geno)")
    p.add_argument("-i", "--ind-map", type=Path, help=".ind-style file: id sex pop")
    p.add_argument("--ind-as-pop", action="store_true", help="use sample ID as population")
    p.add_argument("-r", "--ref", help="add a synthetic all-ref sample with this name")
    p.add_argument("-g", "--genetic-map", type=Path, help="PLINK-format genetic map for cM column")
    p.add_argument("--split-multiallelic", action="store_true",
                   help="split multiallelic SNPs into multiple biallelic rows instead of dropping")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s %(message)s")
    if args.ind_map and args.ind_as_pop:
        log.warning("--ind-map takes precedence over --ind-as-pop")
    convert(
        vcf_path=args.vcf,
        out_root=args.out,
        ind_map_path=args.ind_map,
        ind_as_pop=args.ind_as_pop,
        ref_sample=args.ref,
        genetic_map_path=args.genetic_map,
        split_multiallelic=args.split_multiallelic,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
