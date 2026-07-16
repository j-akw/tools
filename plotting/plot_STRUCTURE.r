#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# plot_admixture.R
#
# Plot the results of an ADMIXTURE run across a set of K values as a single
# stacked-barplot figure (one facet row per K), using ggplot2 and a muted
# palette assembled from the wesanderson colour schemes.
#
# ADMIXTURE writes one ancestry-fraction matrix per K, named
#   {prefix}.{K}.Q
# where each row is an individual (in the same order as the PLINK .fam / .bed)
# and each column is that individual's membership fraction for one cluster.
#
# Usage
#   Rscript plot_admixture.R --prefix PATH --fam PATH.fam --K 2,3,4,5 [options]
#
# Required
#   --prefix PATH     Prefix of the .Q files. e.g. "run/mydata" reads
#                     run/mydata.2.Q, run/mydata.3.Q, ...
#   --fam    PATH     The PLINK .fam file used to produce the run (row order
#                     of the .Q files matches this file).
#   --K      LIST     K values to plot. Comma-separated ("2,3,4,5") and/or
#                     ranges ("2-6"). Combinations are fine ("2,4-6").
#
# Options
#   --out    PATH     Output image. Extension picks the device
#                     (.png / .pdf / .tiff / .svg). Default: admixture.png
#   --group-by WHICH  Column of the .fam that holds the group/population label
#                     when no --meta is given. One of: fid (col 1, default),
#                     iid (col 2), none.
#   --meta   PATH     A metadata CSV used to label individuals (overrides
#                     --group-by). Rows are matched to samples by joining the
#                     .fam IID (col 2) to the --meta-id column, and the
#                     --meta-group column supplies the label. Samples with no
#                     match are labelled 'Unassigned'.
#   --meta-id COL     Metadata column matching the .fam IID.
#                     Default: "Sample ID".
#   --meta-group COL  Metadata column to use as the group label.
#                     Default: "Broad Population Group (full)".
#   --group-order SPEC  Order of the groups along the x-axis. Either a
#                     comma-separated list ("European taurine,Indicine,...") or a
#                     path to a text file with one group per line. Groups not
#                     listed are appended (alphabetically) after the listed ones.
#   --label-style S   How population labels are drawn beside the plot.
#                     leader (default): names evenly spaced below the plot and
#                     joined to each group by a thin leader line, so labels never
#                     overlap even when a group has only a couple of individuals.
#                     leader-split: same, but labels alternate above and below
#                     the plot (every other group on top), doubling the room
#                     between neighbours for very crowded runs.
#                     strip: a label directly under each group's own panel.
#                     none: no labels.
#   --label-size N    Text size for the population labels. Default: 2.6.
#   --sort   WHICH    Ordering of individuals within each group.
#                     membership (default) sorts by dominant-cluster fraction
#                     for a clean gradient; fam keeps the .fam file order.
#   --order-k K       Which K's Q matrix defines the membership sort.
#                     Default: the largest K requested.
#   --width  IN       Figure width in inches. Default: auto from sample count.
#   --height IN       Figure height in inches. Default: auto from number of K.
#   --dpi    N        Raster resolution. Default: 300.
#   --palettes LIST   Comma-separated wesanderson palette names to draw colours
#                     from (muted-first by default).
#   --install         Install any missing CRAN packages before plotting.
#
# Example
#   Rscript plot_admixture.R \
#     --prefix results/pop --fam data/pop.fam --K 2-6 --out pop_admixture.png
# ------------------------------------------------------------------------------

# ---- tiny argument parser ----------------------------------------------------

parse_args <- function(args) {
  opt <- list(
    prefix   = NULL,
    fam      = NULL,
    K        = NULL,
    out      = "admixture.png",
    group_by = "fid",
    meta       = NULL,
    meta_id    = "Sample ID",
    meta_group = "Broad Population Group (full)",
    group_order = NULL,
    label_style = "leader",
    label_size  = 2.6,
    sort     = "membership",
    order_k  = NULL,
    width    = NA_real_,
    height   = NA_real_,
    dpi      = 300,
    palettes = "Moonrise2,Cavalcanti1,Royal1,GrandBudapest2,Rushmore1,Zissou1,Darjeeling2,FantasticFox1,Moonrise3,IsleofDogs1",
    install  = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    take <- function() { i <<- i + 1L; args[[i]] }
    switch(a,
      "--prefix"   = { opt$prefix   <- take() },
      "--fam"      = { opt$fam      <- take() },
      "--K"        = { opt$K        <- take() },
      "-K"         = { opt$K        <- take() },
      "--out"      = { opt$out      <- take() },
      "--group-by" = { opt$group_by <- tolower(take()) },
      "--meta"        = { opt$meta        <- take() },
      "--meta-id"     = { opt$meta_id     <- take() },
      "--meta-group"  = { opt$meta_group  <- take() },
      "--group-order" = { opt$group_order <- take() },
      "--label-style" = { opt$label_style <- tolower(take()) },
      "--label-size"  = { opt$label_size  <- as.numeric(take()) },
      "--sort"     = { opt$sort     <- tolower(take()) },
      "--order-k"  = { opt$order_k  <- as.integer(take()) },
      "--width"    = { opt$width    <- as.numeric(take()) },
      "--height"   = { opt$height   <- as.numeric(take()) },
      "--dpi"      = { opt$dpi      <- as.numeric(take()) },
      "--palettes" = { opt$palettes <- take() },
      "--install"  = { opt$install  <- TRUE },
      "--help"     = { cat_usage(); quit(status = 0) },
      "-h"         = { cat_usage(); quit(status = 0) },
      stop("Unknown argument: ", a)
    )
    i <- i + 1L
  }
  opt
}

cat_usage <- function() {
  cat("Usage: Rscript plot_admixture.R --prefix PATH --fam PATH.fam --K 2,3,4,5 [options]\n",
      "See the header of this script for the full list of options.\n", sep = "")
}

# ---- expand a K spec like "2,4-6" into c(2,4,5,6) ----------------------------

parse_k <- function(spec) {
  if (is.null(spec)) stop("--K is required (e.g. --K 2,3,4,5 or --K 2-6)")
  parts <- strsplit(spec, ",", fixed = TRUE)[[1]]
  ks <- integer(0)
  for (p in trimws(parts)) {
    if (grepl("-", p, fixed = TRUE)) {
      rng <- as.integer(strsplit(p, "-", fixed = TRUE)[[1]])
      ks <- c(ks, seq(rng[1], rng[2]))
    } else {
      ks <- c(ks, as.integer(p))
    }
  }
  sort(unique(ks))
}

# ---- packages ----------------------------------------------------------------

load_pkgs <- function(install, need) {
  missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    if (install) {
      message("Installing missing packages: ", paste(missing, collapse = ", "))
      install.packages(missing, repos = "https://cloud.r-project.org")
    } else {
      stop("Missing required packages: ", paste(missing, collapse = ", "),
           "\nInstall them, or re-run with --install.")
    }
  }
  suppressPackageStartupMessages(for (pkg in need) library(pkg, character.only = TRUE))
}

# ---- a muted, extensible colour ramp built from wesanderson ------------------

admix_colours <- function(n, palette_names) {
  base <- unlist(lapply(palette_names, function(p) as.character(wesanderson::wes_palette(p))))
  base <- unique(base)
  if (n <= length(base)) base[seq_len(n)]
  else grDevices::colorRampPalette(base)(n)
}

# ---- read one ADMIXTURE .Q file ----------------------------------------------

read_q <- function(prefix, k) {
  f <- sprintf("%s.%d.Q", prefix, k)
  if (!file.exists(f)) stop("Q file not found: ", f)
  as.matrix(utils::read.table(f, header = FALSE, colClasses = "numeric"))
}

# ---- main --------------------------------------------------------------------

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opt$prefix)) stop("--prefix is required")
  if (is.null(opt$fam))    stop("--fam is required")
  ks <- parse_k(opt$K)

  load_pkgs(opt$install, c("ggplot2", "wesanderson"))

  # --- fam file: FID IID PID MID SEX PHENO ---
  fam <- utils::read.table(opt$fam, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(fam) < 2) stop("The .fam file needs at least two columns (FID, IID).")
  fid <- as.character(fam[[1]])
  iid <- as.character(fam[[2]])
  n   <- nrow(fam)
  id  <- paste(fid, iid, sep = "‖")            # guaranteed-unique row key

  # group label per individual: from a metadata sheet, or a .fam column
  if (!is.null(opt$meta)) {
    meta <- utils::read.csv(opt$meta, check.names = FALSE,
                            stringsAsFactors = FALSE, colClasses = "character")
    if (!opt$meta_id %in% names(meta))
      stop("metadata has no column '", opt$meta_id, "'. Set --meta-id.")
    if (!opt$meta_group %in% names(meta))
      stop("metadata has no column '", opt$meta_group, "'. Set --meta-group.")
    key <- trimws(as.character(meta[[opt$meta_id]]))
    val <- trimws(as.character(meta[[opt$meta_group]]))
    keep <- nzchar(key) & !duplicated(key)          # first row wins per ID
    lut  <- val[keep]; names(lut) <- key[keep]
    group <- unname(lut[trimws(iid)])
    miss  <- is.na(group) | !nzchar(group)
    if (any(miss)) {
      message(sum(miss), " of ", n,
              " samples not found in metadata -> labelled 'Unassigned'")
      group[miss] <- "Unassigned"
    }
    message("matched ", sum(!miss), "/", n, " samples to metadata groups")
    meta_used <- TRUE
  } else {
    group <- switch(opt$group_by,
      "fid"  = fid,
      "iid"  = iid,
      "none" = rep("", n),
      stop("--group-by must be one of: fid, iid, none")
    )
    meta_used <- FALSE
  }

  # --- read every requested Q matrix, checking row counts ---
  qmats <- setNames(lapply(ks, function(k) read_q(opt$prefix, k)), ks)
  for (k in names(qmats)) {
    if (nrow(qmats[[k]]) != n) {
      stop(sprintf("Row mismatch: %s.%s.Q has %d rows but the .fam file has %d.",
                   opt$prefix, k, nrow(qmats[[k]]), n))
    }
  }

  # --- group order: explicit (--group-order) else appearance order ---
  present <- unique(group)
  if (!is.null(opt$group_order)) {
    spec <- if (file.exists(opt$group_order))          # a file (one per line)...
      readLines(opt$group_order, warn = FALSE)
    else strsplit(opt$group_order, ",")[[1]]           # ...or a comma list
    spec <- trimws(spec); spec <- spec[nzchar(spec)]
    unknown <- setdiff(spec, present)
    if (length(unknown))
      message("--group-order names not present in data (ignored): ",
              paste(unknown, collapse = ", "))
    # specified order first, then any remaining groups (e.g. 'Unassigned')
    group_levels <- c(spec[spec %in% present], sort(setdiff(present, spec)))
  } else {
    group_levels <- present
  }

  # --- global individual ordering, then a numeric x position per individual ---
  if (opt$sort == "membership") {
    order_k <- if (is.null(opt$order_k)) max(ks) else opt$order_k
    if (!order_k %in% ks) stop("--order-k must be one of the requested K values.")
    Qref   <- qmats[[as.character(order_k)]]
    dom    <- apply(Qref, 1, which.max)             # dominant cluster per individual
    domval <- apply(Qref, 1, max)                   # its fraction
    ord    <- order(match(group, group_levels), dom, -domval)
  } else if (opt$sort == "fam") {
    ord <- order(match(group, group_levels), seq_len(n))
  } else {
    stop("--sort must be one of: membership, fam")
  }
  xpos <- integer(n); xpos[ord] <- seq_len(n)        # 1..n along the axis

  # --- long data frame across all K (numeric x) ---
  long <- do.call(rbind, lapply(ks, function(k) {
    Q <- qmats[[as.character(k)]]
    data.frame(
      x       = rep(xpos, times = ncol(Q)),
      K       = factor(sprintf("K = %d", k), levels = sprintf("K = %d", ks)),
      cluster = factor(rep(seq_len(ncol(Q)), each = n)),
      prop    = as.vector(Q),
      stringsAsFactors = FALSE
    )
  }))

  # --- each group is contiguous after sorting: get its x-extent and midpoint ---
  gx     <- lapply(group_levels, function(g) range(xpos[group == g]))
  gend   <- vapply(gx, `[`, numeric(1), 2)
  gmid   <- vapply(gx, function(r) mean(r), numeric(1))
  bounds <- utils::head(gend, -1) + 0.5              # separators between groups

  # --- colours: enough for the largest K ---
  pal_names <- trimws(strsplit(opt$palettes, ",", fixed = TRUE)[[1]])
  cols <- admix_colours(max(ks), pal_names)

  label_style <- if (!meta_used && opt$group_by == "none") "none" else opt$label_style
  ksL <- sprintf("K = %d", ks)

  base_theme <- theme_minimal(base_size = 11) +
    theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank(),
      panel.spacing.y  = grid::unit(4, "pt"),
      strip.background = element_blank(),
      strip.text.y     = element_text(size = 10, colour = "grey30", angle = 0),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(6, 12, 6, 12),
      legend.position  = "none"
    )

  # white separators between groups, restricted to the given (bar) panels
  sep_layer <- function(panel_lvls) {
    if (!length(bounds)) return(NULL)
    geom_segment(
      data = expand.grid(x = bounds, panel = factor(ksL, panel_lvls)),
      aes(x = x, xend = x, y = 0, yend = 1),
      inherit.aes = FALSE, colour = "white", linewidth = 0.4)
  }

  # leader lines + labels for a set of groups, tagged to their own facet row
  # (`panel`). Anchors are spread evenly across the full width; the groups keep
  # their x-order, so the leaders fan out without ever crossing. Living in one
  # ggplot means the x scale is shared with the bars automatically.
  track_layers <- function(groups, mids, side, panel_lvls) {
    Gk <- length(groups)
    anchor <- seq(0.5, n + 0.5, length.out = Gk + 2)[2:(Gk + 1)]
    if (side == "bottom") {
      seg <- data.frame(x = c(mids, mids, anchor), xend = c(mids, anchor, anchor),
                        y = c(rep(1, Gk), rep(.82, Gk), rep(.5, Gk)),
                        yend = c(rep(.82, Gk), rep(.5, Gk), rep(.46, Gk)))
      txt <- data.frame(x = anchor, y = .44, pop = groups); hj <- 1
    } else {                                           # top: mirror vertically
      seg <- data.frame(x = c(mids, mids, anchor), xend = c(mids, anchor, anchor),
                        y = c(rep(0, Gk), rep(.18, Gk), rep(.5, Gk)),
                        yend = c(rep(.18, Gk), rep(.5, Gk), rep(.54, Gk)))
      txt <- data.frame(x = anchor, y = .56, pop = groups); hj <- 0
    }
    seg$panel <- factor(side, panel_lvls); txt$panel <- factor(side, panel_lvls)
    list(
      geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
                   inherit.aes = FALSE, colour = "grey55", linewidth = 0.3),
      geom_text(data = txt, aes(x = x, y = y, label = pop), inherit.aes = FALSE,
                angle = 90, hjust = hj, vjust = 0.5,
                size = opt$label_size, colour = "grey25")
    )
  }

  if (label_style %in% c("leader", "leader-split")) {
    G      <- length(group_levels)
    split  <- label_style == "leader-split" && G >= 2
    panels <- if (split) c("top", ksL, "bottom") else c(ksL, "bottom")
    long$panel <- factor(long$K, levels = panels)

    if (split) {
      bi <- seq(1, G, by = 2); ti <- seq(2, G, by = 2)   # bottom / top groups
      trk <- c(track_layers(group_levels[bi], gmid[bi], "bottom", panels),
               track_layers(group_levels[ti], gmid[ti], "top",    panels))
    } else {
      trk <- track_layers(group_levels, gmid, "bottom", panels)
    }

    strip_lab <- stats::setNames(panels, panels)         # blank the track strips
    strip_lab[c("top", "bottom")] <- ""

    p <- ggplot(long, aes(x = x, y = prop, fill = cluster)) +
      geom_col(width = 1) + sep_layer(panels) + trk +
      scale_fill_manual(values = cols, guide = "none") +
      scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = NULL, title = NULL) +
      facet_grid(rows = vars(panel), labeller = as_labeller(strip_lab)) +
      coord_cartesian(clip = "off") + base_theme
    n_rows <- length(ks) + (if (split) 2 else 1)
    height <- if (is.na(opt$height)) 0.85 * n_rows + 0.6 else opt$height

  } else if (label_style == "strip") {
    # a label directly under each group's own panel (overlaps when groups small)
    long$group <- factor(group[match(long$x, xpos)], levels = group_levels)
    p <- ggplot(long, aes(x = factor(x), y = prop, fill = cluster)) +
      geom_col(width = 1) +
      scale_fill_manual(values = cols, guide = "none") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = NULL, y = NULL, title = NULL) +
      facet_grid(K ~ group, scales = "free_x", space = "free_x", switch = "x") +
      base_theme +
      theme(panel.spacing.x = grid::unit(1.5, "pt"),
            strip.placement = "outside",
            strip.text.x    = element_text(size = 9, colour = "grey30",
                                           angle = 90, hjust = 1, vjust = 0.5))
    height <- if (is.na(opt$height)) 0.9 * length(ks) + 1.4 else opt$height

  } else {                                            # no labels
    long$panel <- factor(long$K, levels = ksL)
    p <- ggplot(long, aes(x = x, y = prop, fill = cluster)) +
      geom_col(width = 1) + sep_layer(ksL) +
      scale_fill_manual(values = cols, guide = "none") +
      scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = NULL, title = NULL) +
      facet_grid(rows = vars(panel)) + base_theme
    height <- if (is.na(opt$height)) 0.9 * length(ks) + 0.6 else opt$height
  }

  width <- if (is.na(opt$width)) max(6, min(40, n * 0.035)) else opt$width
  ggsave(opt$out, plot = p, width = width, height = height, dpi = opt$dpi,
         limitsize = FALSE)
  message(sprintf("Wrote %s  (%.1f x %.1f in, %d individuals, %d groups, K = %s)",
                  opt$out, width, height, n, length(group_levels),
                  paste(ks, collapse = ",")))
}

main()
