# admixtools

A small collection of Python utilities for preparing and manipulating
[EIGENSTRAT](https://reich.hms.harvard.edu/software/InputFileFormats)
datasets (the `.geno` / `.snp` / `.ind` triple used by ADMIXTOOLS,
smartpca, and related population-genetics software).

All scripts are Python 3 and self-contained command-line tools. Two of
them (`split_eigen_by_pop.py` and, optionally, `vcf2eigenstrat.py`) use
`numpy` for speed; the rest are pure standard library.

## The EIGENSTRAT format in one paragraph

A dataset is three files sharing a prefix:

- **`.geno`** — a matrix with one row per SNP and one character per
  individual. Each character is the count of *reference* alleles:
  `2` = hom-ref, `1` = het, `0` = hom-alt, `9` = missing. Rows correspond
  to `.snp` lines; columns correspond to `.ind` lines.
- **`.snp`** — one row per SNP: id, chromosome, genetic position (cM),
  physical position (bp), reference allele, alternate allele.
- **`.ind`** — one row per individual: id, sex, population label.

Because the `.geno` matrix is positional, any tool that adds, removes, or
reorders individuals must edit `.geno` columns and `.ind` rows in lockstep.

---

## Tools

### `vcf2eigenstrat.py` — convert a VCF into EIGENSTRAT

Reads a `.vcf` or `.vcf.gz` and writes `<out>.geno/.snp/.ind`. Pure stdlib
parsing (no `pysam`); uses `numpy` for a vectorized fast path when
available. Keeps only biallelic SNPs by default; indels and multiallelic
sites are dropped (or, with `--split-multiallelic`, expanded into biallelic
rows).

```bash
python vcf2eigenstrat.py -v input.vcf.gz -o mydata \
    -i sample_metadata.ind \
    -g genetic_map.map \
    --split-multiallelic
```

Key options:

| Option | Purpose |
| --- | --- |
| `-v/--vcf` | input VCF (`.vcf` or `.vcf.gz`) |
| `-o/--out` | output prefix (writes `.geno/.snp/.ind`) |
| `-i/--ind-map` | `id sex pop` file to fill the `.ind` columns |
| `--ind-as-pop` | use each sample's ID as its population label |
| `-r/--ref` | prepend a synthetic all-reference individual |
| `-g/--genetic-map` | PLINK-format map for the cM column (interpolated) |
| `--split-multiallelic` | keep multiallelic sites as multiple biallelic rows |

To restrict to a region, pre-slice the VCF with `bcftools view -r ...`
first.

### `split_eigen_by_pop.py` — one dataset per population

Given a flat list of sample IDs, looks up each ID's population in the
`.ind` file, groups the IDs by population, and writes one EIGENSTRAT
dataset per population. Each per-pop dataset contains that population's
listed individuals **first**, followed by a fixed set of reference
populations — the layout many f-statistic / qpAdm workflows expect.

Uses a `numpy` memmap over the `.geno` file and vectorized column
extraction, so it stays fast on large (multi-GB) matrices.

```bash
python split_eigen_by_pop.py \
    -g data.geno -s data.snp -i data.ind \
    -l selected_individuals.txt \
    -r ref_pops.txt \
    -o selected_datasets/
```

Output goes to `selected_datasets/<POP>/<POP>.{geno,snp,ind}`. Since the
`.snp` file is identical across outputs, `--snp-mode` controls whether it
is `symlink`ed (default), `hardlink`ed, or `copy`ied. `--dry-run` prints
the plan without writing.

### `eigenutils.py` — check / extract / remove individuals

A pure-Python multi-tool for individual-level edits on an existing
EIGENSTRAT dataset. One of three mutually exclusive modes:

- `-C/--Check <other.ind>` — report individuals shared between the input
  `.ind` and another `.ind` (duplicate detection across datasets).
- `-E/--Extract` — write a new dataset containing only the named
  individuals.
- `-R/--Remove` — write a new dataset with the named individuals removed.

Individuals are named either inline (`-S ID`, repeatable) or via a list
file (`-L list.txt`).

```bash
# Which individuals appear in both datasets?
python eigenutils.py -g a.geno -s a.snp -i a.ind -C b.ind

# Pull three individuals into a new dataset
python eigenutils.py -g data.geno -s data.snp -i data.ind \
    -E -S IND1 -S IND2 -S IND3 -o subset

# Drop everyone listed in a file
python eigenutils.py -g data.geno -s data.snp -i data.ind \
    -R -L drop_these.txt -o pruned
```

### `update_ind.py` — relabel populations in a `.ind`

Rewrites the population column of a `.ind` file from a two-column
`id  new_pop` mapping. Individuals not in the map pass through unchanged.
Only the `.ind` is touched — the `.geno` and `.snp` are unaffected, since
relabelling does not move any columns.

```bash
# Write a new .ind
python update_ind.py -i input.ind -m pop_map.txt -o relabelled.ind

# Or edit in place (a .bak backup is written first)
python update_ind.py -i input.ind -m pop_map.txt --in-place
```

---

## Typical pipeline

```
VCF ──vcf2eigenstrat.py──▶ .geno/.snp/.ind
                               │
              update_ind.py ◀──┤  (fix / assign population labels)
                               │
        eigenutils.py -C ◀─────┤  (sanity-check for overlaps with other panels)
                               │
   split_eigen_by_pop.py  ◀────┘  (carve per-population analysis datasets)
   / eigenutils.py -E/-R          (or a single subset / prune)
```

---
