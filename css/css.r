# =============================================================================
# Composite Selection Signals (CSS)
# -----------------------------------------------------------------------------
# Rank-based composite selection signal of Randhawa et al. (2014) BMC Genetics
# 15:34. Given >= 2 (typically 3+) constituent selection statistics computed at
# a set of genomic positions, CSS:
#
#   1. ranks each statistic genome-wide,
#   2. maps fractional ranks R'/(n+1) to standard-normal scores via qnorm,
#   3. averages those scores across statistics at each position,
#   4. converts the mean score to a one-sided upper-tail p-value under the null
#      mean-Z ~ N(0, 1/m), and reports CSS = -log10(p).
#
# Positions are then smoothed within chromosome over a physical sliding window,
# and significant regions are called from the top tail of the smoothed signal.
#
# Species-agnostic: chromosome identifiers are arbitrary labels; nothing assumes
# a particular count, naming, or ordering. Smoothing and region calling never
# cross chromosome boundaries.
#
# ---- Orientation --------------------------------------------------------------
# The rank->normal-score step is one-sided: larger statistic -> larger CSS
# contribution. Each statistic must therefore be oriented so that "more evidence
# of selection" = "larger value". Use the `direction` argument:


#   "high"  value is already oriented (e.g. Fst, XP-EHH toward target)   -> x
#   "low"   smaller means more selection (e.g. a per-SNP p-value)         -> -x
#   "abs"   two-sided, magnitude matters (unpolarised dDAF/dSAF, iHS)     -> |x|


# Because ranks are invariant to monotone rescaling, you do NOT need to
# standardise inputs to N(0,1) beforehand; only orientation matters.
#
# ---- Usage --------------------------------------------------------------------


#   source("css.R")
#   res <- run_css(df, stat_cols = c("fst", "xpehh", "dSAF"),
#                  chrom_col = "chr", pos_col = "pos")
#   head(res$snps)      # per-SNP: z_*, css_meanZ, css_p, css, *_smooth
#   res$regions         # called regions (chr, start, end, peak, ...)
#   manhattan_css(res$snps, threshold = res$threshold)


# =============================================================================




# --- orient a single statistic so that larger = more selection ----------------
.orient_stat <- function(x, direction) {
  switch(direction,
         high = x,
         low  = -x,
         abs  = abs(x),
         stop("direction must be one of 'high', 'low', 'abs'"))
}




# --- per-SNP composite score --------------------------------------------------
# Returns df with added columns:
#   z_<stat>     normal score for each constituent statistic
#   css_n_tests  number of statistics contributing at that SNP (m_j)
#   css_meanZ    mean normal score across contributing statistics
#   css_p        one-sided p-value under mean-Z ~ N(0, 1/m_j)
#   css          -log10(css_p)
css_scores <- function(df,
                       stat_cols,
                       direction   = "high",
                       ties_method = "average",
                       min_tests   = length(stat_cols)) {
  
  stopifnot(is.data.frame(df))
  if (length(stat_cols) < 2L)
    stop("CSS needs at least two constituent statistics.")
  miss <- setdiff(stat_cols, names(df))
  if (length(miss))
    stop("Columns not found in df: ", paste(miss, collapse = ", "))
  
  m <- length(stat_cols)
  
  # Match `direction` to each statistic (scalar recycles; else must be named).
  if (length(direction) == 1L) {
    direction <- setNames(rep(direction, m), stat_cols)
  } else {
    if (is.null(names(direction)))
      stop("When 'direction' has length > 1 it must be named by stat_cols.")
    direction <- direction[stat_cols]
    if (anyNA(names(direction)))
      stop("'direction' must name every column in stat_cols.")
  }
  
  # Per-statistic normal scores from fractional ranks. NAs are ranked out
  # (na.last = "keep") and the fractional-rank denominator uses that
  # statistic's own count of non-missing SNPs.
  z <- matrix(NA_real_, nrow = nrow(df), ncol = m,
              dimnames = list(NULL, stat_cols))
  for (s in stat_cols) {
    x  <- .orient_stat(df[[s]], direction[[s]])
    ok <- !is.na(x)
    ni <- sum(ok)
    if (ni < 2L) next
    r  <- rank(x[ok], ties.method = ties_method)   # 1..ni; larger = more selection
    z[ok, s] <- qnorm(r / (ni + 1))                # fractional rank -> N(0,1)
  }
  
  # Mean score across available statistics at each SNP.
  mj    <- rowSums(!is.na(z))
  meanZ <- ifelse(mj > 0, rowSums(z, na.rm = TRUE) / mj, NA_real_)
  
  # Upper-tail p under mean-Z ~ N(0, 1/mj). Work in the log domain so the deep
  # tail gives finite -log10(p) instead of collapsing to Inf at p == 0.
  logp <- pnorm(sqrt(mj) * meanZ, lower.tail = FALSE, log.p = TRUE)
  css  <- -logp / log(10)
  
  drop <- mj < min_tests                            # too few statistics -> undefined
  meanZ[drop] <- NA_real_
  css[drop]   <- NA_real_
  
  out <- df
  out[paste0("z_", stat_cols)] <- as.data.frame(z)
  out$css_n_tests <- mj
  out$css_meanZ   <- meanZ
  out$css_p       <- ifelse(drop, NA_real_, exp(logp))
  out$css         <- css
  out
}




# --- within-chromosome sliding-window smoothing -------------------------------
# Averages `value_col` over all SNPs whose position lies within +/- window/2 of
# each SNP, on the same chromosome. Windows with fewer than `min_snps`
# contributing SNPs are set to NA. O(n log n) via sorted cumulative sums.
smooth_css <- function(df,
                       value_col = "css",
                       chrom_col = "chr",
                       pos_col   = "pos",
                       window    = 1e6,
                       min_snps  = 5) {
  
  for (cc in c(value_col, chrom_col, pos_col))
    if (!cc %in% names(df)) stop("Column not found: ", cc)
  
  half <- window / 2
  n    <- nrow(df)
  sm   <- rep(NA_real_,    n)
  cnt  <- rep(NA_integer_, n)
  
  for (ch in unique(df[[chrom_col]])) {
    idx <- which(df[[chrom_col]] == ch)
    o   <- order(df[[pos_col]][idx])       # sort SNPs on this chromosome
    ps  <- df[[pos_col]][idx][o]
    vs  <- df[[value_col]][idx][o]
    
    # NA-safe running sums: NA values contribute 0 to the sum and 0 to the count.
    csv <- c(0, cumsum(ifelse(is.na(vs), 0, vs)))
    csn <- c(0, cumsum(as.numeric(!is.na(vs))))
    
    hi <- findInterval(ps + half,        ps)          # last SNP with pos <= hi bound
    lo <- findInterval(ps - half - 1e-9, ps)          # last SNP strictly < lo bound
    
    wsum <- csv[hi + 1] - csv[lo + 1]
    wcnt <- csn[hi + 1] - csn[lo + 1]
    res  <- ifelse(wcnt >= min_snps, wsum / wcnt, NA_real_)
    
    sm[idx[o]]  <- res
    cnt[idx[o]] <- wcnt
  }
  
  out <- df
  out[[paste0(value_col, "_smooth")]] <- sm
  out$window_n_snps <- cnt
  out
}




# --- region calling from the smoothed signal ----------------------------------
# SNPs above the top-`top_frac` genome-wide threshold are "hits"; hits within
# `merge_gap` of each other on a chromosome form a cluster; clusters with at
# least `min_cluster_snps` hits become regions, extended by `flank` each side.
css_regions <- function(df,
                        value_col        = "css_smooth",
                        chrom_col        = "chr",
                        pos_col          = "pos",
                        top_frac         = 0.001,
                        min_cluster_snps = 3,
                        merge_gap        = 1e6,
                        flank            = 0.5e6,
                        threshold        = NULL) {
  
  v <- df[[value_col]]
  if (is.null(threshold))
    threshold <- as.numeric(quantile(v, 1 - top_frac, na.rm = TRUE))
  
  hit  <- !is.na(v) & v >= threshold
  hits <- df[hit, c(chrom_col, pos_col, value_col)]
  names(hits) <- c("chr", "pos", "val")
  
  empty <- structure(
    data.frame(chr = character(), start = numeric(), end = numeric(),
               n_snps = integer(), peak_pos = numeric(), peak_value = numeric(),
               stringsAsFactors = FALSE),
    threshold = threshold)
  if (!nrow(hits)) return(empty)
  
  regions <- do.call(rbind, lapply(split(hits, hits$chr), function(h) {
    h <- h[order(h$pos), ]
    h$cluster <- cumsum(c(TRUE, diff(h$pos) > merge_gap))   # single-linkage by gap
    do.call(rbind, lapply(split(h, h$cluster), function(g) {
      if (nrow(g) < min_cluster_snps) return(NULL)
      pk <- which.max(g$val)
      data.frame(chr        = g$chr[1],
                 start      = max(1, min(g$pos) - flank),
                 end        = max(g$pos) + flank,
                 n_snps     = nrow(g),
                 peak_pos   = g$pos[pk],
                 peak_value = g$val[pk],
                 stringsAsFactors = FALSE)
    }))
  }))
  
  if (is.null(regions)) regions <- empty else rownames(regions) <- NULL
  attr(regions, "threshold") <- threshold
  regions
}




# --- end-to-end wrapper -------------------------------------------------------
run_css <- function(df,
                    stat_cols,
                    chrom_col        = "chr",
                    pos_col          = "pos",
                    direction        = "high",
                    ties_method      = "average",
                    min_tests        = length(stat_cols),
                    window           = 1e6,
                    min_snps         = 5,
                    top_frac         = 0.001,
                    min_cluster_snps = 3,
                    merge_gap        = 1e6,
                    flank            = 0.5e6) {
  
  s <- css_scores(df, stat_cols, direction, ties_method, min_tests)
  s <- smooth_css(s, "css", chrom_col, pos_col, window, min_snps)
  r <- css_regions(s, "css_smooth", chrom_col, pos_col,
                   top_frac, min_cluster_snps, merge_gap, flank)
  list(snps = s, regions = r, threshold = attr(r, "threshold"))
}




# --- optional FDR (q-values) --------------------------------------------------
# The paper calibrates the empirical CSS p-values with ConReg-R (Li et al. 2011)
# BEFORE fdrtool; that calibration step is not on CRAN, so this helper runs
# fdrtool directly on css_p. Feed it calibrated p-values to match the paper.
css_fdr <- function(css_p) {
  if (!requireNamespace("fdrtool", quietly = TRUE))
    stop("Install 'fdrtool' to estimate q-values.")
  p <- css_p[!is.na(css_p)]
  q <- rep(NA_real_, length(css_p))
  q[!is.na(css_p)] <- fdrtool::fdrtool(p, statistic = "pvalue",
                                       plot = FALSE, verbose = FALSE)$qval
  q
}




# --- base-R Manhattan plot of the (smoothed) CSS ------------------------------
manhattan_css <- function(df,
                          value_col = "css_smooth",
                          chrom_col = "chr",
                          pos_col   = "pos",
                          threshold = NULL,
                          col       = c("grey30", "grey60"),
                          ...) {
  
  d    <- df[!is.na(df[[value_col]]), ]
  chrs <- unique(d[[chrom_col]])
  num  <- suppressWarnings(as.numeric(as.character(chrs)))     # numeric order if possible
  chrs <- if (!anyNA(num)) chrs[order(num)] else chrs
  d[[chrom_col]] <- factor(d[[chrom_col]], levels = chrs)
  d    <- d[order(d[[chrom_col]], d[[pos_col]]), ]
  
  offset <- 0
  xpos   <- numeric(nrow(d))
  cols   <- character(nrow(d))
  ticks  <- numeric(length(chrs))
  for (i in seq_along(chrs)) {
    sel <- d[[chrom_col]] == chrs[i]
    pp  <- d[[pos_col]][sel]
    xpos[sel] <- pp + offset
    cols[sel] <- col[(i - 1L) %% length(col) + 1L]
    ticks[i]  <- offset + (min(pp) + max(pp)) / 2
    offset    <- offset + max(pp) + 1
  }
  
  plot(xpos, d[[value_col]], pch = 20, cex = 0.5, col = cols,
       xaxt = "n", xlab = "Chromosome", ylab = "CSS  (-log10 p)", ...)
  axis(1, at = ticks, labels = chrs, las = 2, cex.axis = 0.7)
  if (!is.null(threshold)) abline(h = threshold, col = "red", lty = 2)
  invisible(NULL)
}
