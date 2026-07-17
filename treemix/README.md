# treemix

Utilities for running and plotting [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home)
analyses: one Python script to build TreeMix input from PLINK allele-frequency
output, and an R plotting library plus driver script that turn a directory of
TreeMix runs into publication-ready figures.

The R side is a rewrite of Pickrell & Pritchard's `plotting_funcs.R`. The tree
layout maths is unchanged; what is new is the styling, the replicate-aware run
discovery, and the vector-export helpers.

---

## `plink2treemix.py` — PLINK frequencies to TreeMix input

Converts the `.frq.strat` file produced by `plink --freq --within pops.clust`
into the whitespace-delimited allele-count table TreeMix expects: one header
line of population names, then one line per SNP with `refCount,altCount` per
population.

```bash
plink --bfile mydata --freq --within pops.clust --out mydata
python plink2treemix.py mydata.frq.strat.gz -o treemix_input.gz
treemix -i treemix_input.gz -m 1 -o out
```

The input needs the `SNP`, `CLST`, `MAC`, and `NCHROBS` columns; counts are
written as `NCHROBS - MAC` (reference) and `MAC` (minor). Reads and writes
plain or gzipped files, and `-` writes to stdout. Population and SNP counts go
to stderr so they don't contaminate a piped stream.

| Option | Purpose |
| --- | --- |
| `-o/--out` | output file; `.gz` triggers gzip, `-` is stdout (default) |
| `--keep-missing` | emit `0,0` for a population missing a SNP instead of dropping the SNP |
| `--sort-pops` | alphabetise the header; default is order of first appearance |

By default any SNP absent from at least one population is dropped, which is
usually what you want — TreeMix treats `0,0` as no data rather than as an
error, but a matrix full of them will skew the covariance estimates.

Requires Python 3 and `pandas`.

## `treemix_plotting_funcs.R` — the plotting library

Sourced by the driver; also usable on its own. Every plotter takes a *stem* —
the filename prefix before `.vertices.gz`, `.edges.gz`, `.cov.gz`, and so on —
and errors early with a clear message if any required file is missing.

```r
source("treemix_plotting_funcs.R")

style <- treemix_style()                       # Zissou1 palettes by default
generate_pop_order("out.m1.rep01", "pop_order.txt")

save_plot("tree.pdf", width = 9, height = 6.5, style = style, {
  plot_tree("out.m1.rep01", style = style, title = "TreeMix — m = 1")
})
```

- `plot_tree(stem, ...)` — the tree, with migration arrows coloured by weight,
  a drift-parameter axis, a 10 s.e. scale bar, and a migration-weight legend.
- `plot_resid(stem, pop_order, ...)` — the residual covariance matrix
  (data − model) as a lower-triangular heatmap scaled in SE units.
- `find_treemix_runs(dir, prefix, ...)` — scans a directory for `.llik` files
  and returns a data frame of `m`, `rep`, `stem`, and `llik`.
- `best_replicate(runs_df, m)` — the highest-likelihood replicate for a given
  `m`, or `NULL` if no complete file set exists for it.
- `generate_pop_order(stem, outfile)` — writes the population order from a
  `.cov.gz` header, which `plot_resid` needs.
- `save_plot(file, expr, ...)` — renders `expr` to `cairo_pdf`, `png`, or
  `svg`, picked from the file extension.
- `treemix_style(...)` / `treemix_par(style, mar)` — colours, fonts, and base-R
  graphics parameters. Pass any colour vector to override the defaults, e.g.
  `treemix_style(mig_palette = viridisLite::viridis(250))`.

Depends on `RColorBrewer` and `wesanderson`; both are installed on first use if
absent.

## `plot_treemix.R` — the driver

Runs the whole figure set for a directory of TreeMix output. Everything you
need to change lives in the `CONFIG` block at the top of the file; nothing
below it should need editing for a new dataset.

```bash
Rscript plot_treemix.R
```

It writes four figures into `plots_dir`:

- `TreeMix_Trees_Facet` — a facet panel of the best-likelihood tree per `m`
- `TreeMix_Residuals_Facet` — the matching residual panel
- `TreeMix_Tree_m<best_m>_best` — the main-figure tree for your chosen `m`
- `TreeMix_Residuals_m<best_m>_best` — residuals for that same replicate

Key `CONFIG` fields:

| Field | Purpose |
| --- | --- |
| `runs_dir` / `plots_dir` | where the TreeMix output lives, where figures go |
| `prefix`, `m_prefix`, `has_replicates`, `rep_pattern` | the filename convention (see below) |
| `m_range` | which `m` values appear in the facet panels |
| `best_m` | the `m` used for the main figures |
| `output_ext`, `png_dpi` | `"pdf"` (vector), `"png"`, or `"svg"` |
| `mig_palette`, `resid_palette`, `font_family` | aesthetics; `NULL` keeps the defaults |

Filenames are assumed to follow `<prefix>.<m_prefix><m>[<rep_suffix>].<ext>`,
so `prefix = "dml_data"`, `m_prefix = ""`, `has_replicates = TRUE` matches
`dml_data.0.rep01.cov.gz`, while `prefix = "out"`, `m_prefix = "m"` matches
`out.m0.rep01.cov.gz`. Set `has_replicates = FALSE` for unreplicated runs like
`run.0.cov.gz`.

Replicates are handled for you: for each `m` the driver reads every `.llik`,
picks the highest, and confirms the full set of plot files exists before using
it. An `m` with no usable replicate becomes a greyed-out "Missing" panel rather
than aborting the run, and a plot that throws becomes a red "Error" panel — so
one bad run never costs you the whole figure. Choose `best_m` yourself, from an
[OptM](https://cran.r-project.org/package=OptM) Deltam peak or a variance
criterion.

---