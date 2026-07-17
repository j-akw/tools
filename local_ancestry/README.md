# local_ancestry

A small collection of command-line tools for working with
[GNOMIX](https://github.com/AI-sandbox/gnomix) local-ancestry inference
output: turning per-window posterior probabilities into confident ancestry
intervals, masking a VCF down to those intervals, and plotting a per-sample
karyogram.

Two scripts are Python 3 (standard library only — no `numpy`, no `pysam`);
the plotting script is R and needs `ggplot2` and `dplyr`. Writing a
bgzipped VCF additionally requires `bgzip` (htslib) on `PATH`.

## The GNOMIX output files in one paragraph

A GNOMIX run produces two files per chromosome that these tools read:

- **`.msp`** — one row per window: `chrom`, `spos`, `epos`, `sgpos`,
  `egpos`, `n snps`, then two columns per individual (`SAMPLE.0` and
  `SAMPLE.1`, one per haplotype) holding the *hard call* — an integer
  ancestry code. Line 1 is a `#Subpopulation order/codes` comment mapping
  labels to codes; line 2 is the column header.
- **`.fb`** — the same windows, but with the *posterior probabilities*
  behind those calls: four leading columns (`chrom`, physical position,
  genetic position, marker index) followed by one column per
  `SAMPLE:::hap{1,2}:::ANCESTRY` combination. Line 1 lists the reference
  panel populations; line 2 is the column header.

The `.fb` and `.msp` for a given chromosome have the same windows in the
same order, which is what lets `gnomix2bed.py` attach genomic coordinates
to `.fb` rows.

---

## Tools

### `gnomix2bed.py` — confident ancestry calls → per-sample BED

Reads a matched `.fb` / `.msp` pair and, for every
(individual, haplotype, ancestry) combination, walks the windows and emits
runs where the posterior exceeds a threshold (default `> 0.99`). Contiguous
passing windows are merged into a single interval; the BED score column
holds the *minimum* posterior across the merged run, so a coarse interval
still reports its weakest window.

```bash
python gnomix2bed.py \
    --fb chr1.fb --msp chr1.msp \
    --out-dir beds/ \
    --threshold 0.99
```

Writes `beds/<SAMPLE>.bed`, one file per individual, with columns:

| Column | Meaning |
| --- | --- |
| `chrom` | chromosome |
| `start` | 0-based interval start (BED convention) |
| `end` | interval end |
| `name` | `<ANCESTRY>_hap{1,2}`, e.g. `AFT_hap1` |
| `score` | minimum posterior across the merged interval |

Each invocation handles one chromosome's file pair. To cover a genome, run
it per chromosome into separate directories and concatenate the matching
per-sample BEDs.

### `maskvcf.py` — keep genotypes only in confident ancestry tracts

Masks a VCF so that named individuals retain their genotypes **only** where
*both* haplotypes were confidently called as a target ancestry; everywhere
else their `GT` becomes `./.`. Other FORMAT fields, other samples, and all
header lines pass through untouched.

The rule is deliberately strict: for each sample the tool intersects that
ancestry's `_hap1` intervals with its `_hap2` intervals, per chromosome, and
keeps a site only if it falls inside the intersection. A site where one
haplotype is confidently `AFT` and the other is anything else is masked.

```bash
python maskvcf.py \
    --vcf input.vcf.gz \
    --out masked.vcf.gz \
    --ancestry AFT \
    --sample BAR4050=beds/BAR4050.bed \
    --sample BAR4051=beds/BAR4051.bed
```

`--sample` is repeatable and takes `SAMPLE_ID=BED_PATH`; the BED is the
output of `gnomix2bed.py`, and `--ancestry` must match the `NAME` prefix in
it (default `AFT`). An `--out` ending in `.gz` is piped through `bgzip`, so
the result can be `tabix`-indexed. The tool exits with an error if any named
sample is missing from the VCF, and reports per-sample kept/masked counts on
stderr.

### `plot_karyogram.r` — per-sample ancestry karyogram

Reads the `.msp` hard calls for one individual across all chromosomes and
draws a two-bar-per-chromosome karyogram (one bar per haplotype), coloured
by ancestry, as a vector PDF.

```bash
# All chromosomes for one sample
Rscript plot_karyogram.r -d msp_files/ -s BRH03

# Custom ancestry colours
Rscript plot_karyogram.r -d msp_files/ -s BRH03 -0 '#2166AC' -1 '#D6604D'

# Highlight one ancestry, fade the rest
Rscript plot_karyogram.r -d msp_files/ -s BRH03 -hl 1 --fade_alpha 0.10
```

Key options:

| Option | Purpose |
| --- | --- |
| `-d/--input_dir` | directory of `.msp` files (required) |
| `-s/--sample` | sample ID, without the `.0`/`.1` haplotype suffix (required) |
| `-o/--output` | output PDF (default `<sample>_karyogram.pdf`) |
| `-p/--pattern` | glob for `.msp` files (default `*.msp`) |
| `-0` / `-1` / `-2` | hex colour for ancestry code 0 / 1 / 2 |
| `-hl/--highlight` | show one ancestry code at full opacity, fade the others |
| `--fade_alpha` | opacity for faded ancestries (default `0.12`) |
| `--figwidth` / `--figheight` | figure size in mm (default 260 × 200) |
| `--title` / `--no-title` | set or hide the plot title |
| `--no-legend` / `--no-chr` | hide the legend / chromosome labels |

Ancestry labels come from the `#Subpopulation` line of the first `.msp`
read; codes without a label fall back to `Ancestry <n>`. If the sample ID is
not found, the script lists partial matches (or the first 30 sample IDs) to
help catch typos. Run `Rscript plot_karyogram.r --help` for the full list.

---

## Typical pipeline

```
GNOMIX .msp + .fb
     │
     ├──plot_karyogram.r──▶ <sample>_karyogram.pdf   (eyeball the calls)
     │
     └──gnomix2bed.py──▶ beds/<SAMPLE>.bed           (posterior > 0.99 tracts)
                              │
                              └──maskvcf.py──▶ masked.vcf.gz
                                 (GT kept only where hap1 ∩ hap2 == ancestry)
```

---
