# =============================================================================
# Prepare constituent selection statistics for CSS
# -----------------------------------------------------------------------------
# Reads Fst (VCFtools), dSAF (VCFtools --freq), and XP-EHH (selscan) into tidy
# per-SNP tables and merges them on (chromosome, position) for input to css.R.
#
# Design notes vs. an ad-hoc script:
#   - paths and key column names are arguments, not hard-coded,
#   - chromosome is kept as CHARACTER throughout (never coerced to numeric), so
#     "X", scaffolds, and "chr1"-style names survive -> species-agnostic,
#   - tables are merged on the (chromosome, position) pair, not a pasted
#     "chr:pos" string, so a type mismatch in one source can't silently drop
#     every SNP,
#   - columns are selected by name, so joins reordering columns can't scramble
#     which statistic is which,
#   - readers are pure (return a data frame); checkpoint to disk yourself
#     between steps if you want crash safety (write.table / fwrite).
# =============================================================================


library(data.table)
library(dplyr)
library(tidyr)




# --- Fst ----------------------------------------------------------------------
# Generate with VCFtools, e.g.:
#   vcftools --vcf wgs_data.vcf \
#            --weir-fst-pop target_pop.txt --weir-fst-pop ref_pop.txt \
#            --out target_ref_Fst
# -> target_ref_Fst.weir.fst  (CHROM  POS  WEIR_AND_COCKERHAM_FST)
read_fst <- function(path) {
  fread(path, na.strings = c("NA", "nan", "-nan")) |>
    as_tibble() |>
    rename_with(~ c("chromosome", "position", "fst"), .cols = 1:3) |>
    transmute(chromosome = as.character(chromosome),
              position   = as.integer(position),
              fst        = as.numeric(fst))
}




# --- dSAF (directional change in selected-allele frequency) -------------------
# Generate per-population allele frequencies with VCFtools --freq on BIALLELIC
# sites (add --min-alleles 2 --max-alleles 2), once per population:
#   vcftools --vcf wgs_data.vcf --keep target_pop.txt --freq --out target
#   vcftools --vcf wgs_data.vcf --keep ref_pop.txt    --freq --out ref
# dSAF = (freq of the target's major allele in target)
#        - (freq of that SAME allele in the reference).
# NOTE: pass TWO DIFFERENT files here (target vs reference) -- the common
# copy-paste bug is pointing both reads at the same .frq.
read_deltaSAF <- function(target_frq, reference_frq, standardize = TRUE) {
  
  parse_frq <- function(p) {
    fread(p, header = FALSE, skip = 1,
          col.names = c("chromosome", "position",
                        "n_alleles", "n_chr", "af1", "af2")) |>
      as_tibble() |>
      separate(af1, c("allele1", "freq1"), sep = ":", convert = TRUE) |>
      separate(af2, c("allele2", "freq2"), sep = ":", convert = TRUE) |>
      mutate(chromosome = as.character(chromosome),
             position   = as.integer(position))
  }
  
  tgt <- parse_frq(target_frq) |>
    mutate(major_allele   = if_else(freq1 >= freq2, allele1, allele2),
           major_freq_tgt = pmax(freq1, freq2)) |>
    select(chromosome, position, major_allele, major_freq_tgt)
  
  ref <- parse_frq(reference_frq) |>
    select(chromosome, position, allele1, freq1, allele2, freq2)
  
  out <- tgt |>
    left_join(ref, by = c("chromosome", "position")) |>
    mutate(major_freq_ref = case_when(major_allele == allele1 ~ freq1,
                                      major_allele == allele2 ~ freq2,
                                      TRUE ~ NA_real_),
           deltaSAF = major_freq_tgt - major_freq_ref) |>
    select(chromosome, position, deltaSAF)
  
  # Standardisation matches Randhawa et al.; it's optional for CSS itself, whose
  # rank transform is invariant to any monotone rescaling. Kept for parity and
  # for any non-CSS use of the column.
  if (standardize) out$deltaSAFz <- as.numeric(scale(out$deltaSAF))
  out
}




# --- XP-EHH -------------------------------------------------------------------
# Per-chromosome with selscan (target as --vcf, reference as --vcf-ref, so
# positive = selection toward the target), then concatenate:
#   selscan --xpehh --vcf target_chr1.vcf.gz --vcf-ref ref_chr1.vcf.gz \
#           --map chr1.map --out chr1_target-ref
#
#   awk 'FNR==1 && NR!=1{next}{print}' *_target-ref.xpehh.out \
#     | sort -k1,1 -k2,2n > all_target-ref.xpehh.out
# `id_sep` splits the locus id into chromosome/position -- change it if your ids
# aren't "chrom-pos", and beware separators that also occur in scaffold names.


read_xpehh <- function(path, id_sep = "-", xpehh_col = "xpehh",
                       standardize = TRUE) {
  x <- fread(path) |> as_tibble()
  names(x)[1] <- "chrpos"
  out <- x |>
    separate(chrpos, c("chromosome", "position"), sep = id_sep, convert = FALSE) |>
    transmute(chromosome = as.character(chromosome),
              position   = as.integer(position),
              xpehh      = as.numeric(.data[[xpehh_col]]))
  if (standardize) out$XPEHHz <- as.numeric(scale(out$xpehh))
  out
}




# --- merge --------------------------------------------------------------------
# Joins any number of the tables above on (chromosome, position). `how="inner"`
# keeps SNPs typed in every input (equivalent to full_join + drop_na, but
# explicit about intent); `how="full"` keeps the union.
merge_css_inputs <- function(..., how = c("inner", "full")) {
  how   <- match.arg(how)
  join  <- if (how == "inner") inner_join else full_join
  parts <- list(...)
  
  merged <- Reduce(function(a, b) join(a, b, by = c("chromosome", "position")),
                   parts)
  if (how == "full") merged <- drop_na(merged)
  
  merged |>
    arrange(chromosome, position) |>          # cosmetic; CSS re-sorts internally
    mutate(SNPid = paste0(chromosome, ":", position)) |>
    relocate(SNPid, .after = position)
}




# --- example ------------------------------------------------------------------
if (sys.nframe() == 0L) {
  fst   <- read_fst("target_ref_Fst.weir.fst")
  saf   <- read_deltaSAF("target.frq", "ref.frq")      # <- two DIFFERENT files
  xpehh <- read_xpehh("all_target-ref.xpehh.out")
  
  df <- merge_css_inputs(fst, saf, xpehh)              # inner join by default
  write.table(df, "merged_statistics.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Hand off to css.R. All three are oriented toward the target population, so
  # every direction is "high" (flip XPEHHz/deltaSAFz to "low" for a reverse
  # scan, or use "abs" for a folded two-sided scan -- keep the two signed stats
  # consistent with each other).
  #
  # source("css.R")
  # res <- run_css(df,
  #                stat_cols = c("fst", "deltaSAFz", "XPEHHz"),
  #                chrom_col = "chromosome", pos_col = "position",
  #                direction = c(fst = "high", deltaSAFz = "high", XPEHHz = "high"))
  # res$regions
}
