#!/usr/bin/env python3
"""
VCF -> Eigenstrat (.geno/.snp/.ind) using cyvcf2 + NumPy
pip install cyvcf2 numpy
"""

import argparse
import sys
from cyvcf2 import VCF
import numpy as np

def get_args():
    p = argparse.ArgumentParser(description="VCF to Eigenstrat via cyvcf2.")
    p.add_argument("--vcf", required=True, help="Input VCF/VCF.gz")
    p.add_argument("--out", required=True, help="Output prefix")
    p.add_argument("--keep", help="File with sample IDs to keep (one per line)")
    p.add_argument("--chunk-size", type=int, default=1000000, help="Variants per flush")
    return p.parse_args()

def read_keep(path):
    if not path:
        return None
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]

def main():
    a = get_args()
    keep = read_keep(a.keep)

    # Subset at the reader level to avoid Python-side indexing
    vcf = VCF(a.vcf, lazy=True, gts012=True, samples=keep if keep else None)

    samples = vcf.samples
    n = len(samples)
    if keep and n == 0:
        raise ValueError("None of the --keep samples are in the VCF")

    # .ind
    with open(f"{a.out}.ind", "w") as f:
        f.writelines(f"{s}\tU\t{s}\n" for s in samples)

    # Open outputs with large buffers
    geno = open(f"{a.out}.geno", "wb", buffering=8 * 1024 * 1024)  # bytes
    snp  = open(f"{a.out}.snp",  "w",  buffering=8 * 1024 * 1024)  # text

    # Map cyvcf2 0/1/2/3 to ASCII bytes for Eigenstrat 2/1/0/9
    # index: 0,1,2,3 -> values: '2','1','0','9'
    encode_ascii = np.frombuffer(b"2109", dtype=np.uint8)

    # Buffers
    geno_buf = bytearray()
    snp_lines = []
    flush_every = max(1, a.chunk_size)

    # Small helpers
    def norm_chrom(c):
        cu = c.upper()
        if cu.startswith("CHR"):
            cu = cu[3:]
        return "23" if cu in ("X", "23") else cu

    wrote = 0
    try:
        for v in vcf:
            # skip multiallelics and indels fast
            if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
                continue

            gt = v.gt_types  # np.ndarray of uint8 (0,1,2,3)
            # encode to ASCII digits without per-element Python work
            encoded = encode_ascii[gt].tobytes()
            geno_buf.extend(encoded)
            geno_buf.append(0x0A)  # '\n'

            sid = v.ID if v.ID else f"{v.CHROM}:{v.POS}"
            snp_lines.append(f"{sid}\t{norm_chrom(v.CHROM)}\t0.0\t{v.POS}")

            wrote += 1
            if wrote % flush_every == 0:
                geno.write(geno_buf)
                geno_buf.clear()
                snp.write("\n".join(snp_lines) + "\n")
                snp_lines.clear()
                sys.stdout.write(f"\rProcessed {wrote} variants")
                sys.stdout.flush()

        if geno_buf:
            geno.write(geno_buf)
        if snp_lines:
            snp.write("\n".join(snp_lines) + "\n")
    finally:
        geno.close()
        snp.close()
        vcf.close()
        sys.stdout.write(f"\rProcessed {wrote} variants\n")

if __name__ == "__main__":
    main()
