# ==============================================================================
# TREEMIX PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------
# A publication-ready replacement for Pickrell & Pritchard's plotting_funcs.R.
#
# Public interface (all take a `stem` = the filename prefix BEFORE .vertices.gz
# etc., so that <stem>.vertices.gz, <stem>.edges.gz etc. exist):
#
#   Plotters:
#     plot_tree(stem, ...)
#     plot_resid(stem, pop_order, ...)
#
#   Run discovery (for replicate-aware plotting):
#     read_treemix_llik(llik_file)     -> numeric log-likelihood
#     find_treemix_runs(dir, ...)      -> data.frame of runs with m, rep, stem, llik
#     best_replicate(runs_df, m)       -> single-row data.frame for best rep
#
#   Utilities:
#     generate_pop_order(stem, outfile)
#     save_plot(file, expr, width, height, dpi, style)
#     treemix_style(mig_palette, resid_palette, family, ...)
#     treemix_par(style, mar)
# ==============================================================================

.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
.need("RColorBrewer")
.need("wesanderson")

# ==============================================================================
# STYLE
# ------------------------------------------------------------------------------
# Default palettes: wesanderson Zissou1 (diverging, blue → yellow → red).
# Override by passing any vector of hex colours:
#   treemix_style(mig_palette = viridisLite::viridis(250),
#                 resid_palette = rev(RColorBrewer::brewer.pal(11, "RdBu")))
# ==============================================================================

treemix_style <- function(mig_palette   = NULL,
                          resid_palette = NULL,
                          family        = "sans",
                          tree_col      = "#2b2b2b",
                          tree_lwd      = 1.6,
                          axis_col      = "#444444",
                          frame_col     = "#888888") {
  if (is.null(mig_palette))
    mig_palette   <- wesanderson::wes_palette("Zissou1", 250, type = "continuous")
  if (is.null(resid_palette))
    resid_palette <- wesanderson::wes_palette("Zissou1", 21,  type = "continuous")

  list(family      = family,
       tree_col    = tree_col,
       tree_lwd    = tree_lwd,
       axis_col    = axis_col,
       frame_col   = frame_col,
       mig_palette = mig_palette,
       resid_pal   = resid_palette)
}

treemix_par <- function(style = treemix_style(), mar = c(4, 2, 2, 2)) {
  par(family   = style$family, mar = mar,
      mgp      = c(2.2, 0.6, 0), tcl = -0.3,
      col.axis = style$axis_col, col.lab = style$axis_col,
      cex.axis = 0.85, cex.lab = 0.95)
}

# ==============================================================================
# COORDINATES  (Pickrell & Pritchard — UNCHANGED)
# ==============================================================================

set_y_coords = function(d){
  i = which(d[,3]=="ROOT")
  y = d[i,8]/ (d[i,8]+d[i,10])
  d[i,]$y = 1-y; d[i,]$ymin = 0; d[i,]$ymax = 1
  c1 = d[i,7]; c2 = d[i,9]
  ni = which(d[,1]==c1)
  ny = d[ni,8]/ (d[ni,8]+d[ni,10])
  d[ni,]$ymin = 1-y; d[ni,]$ymax = 1; d[ni,]$y = 1 - ny*(y)
  ni = which(d[,1]==c2)
  ny = d[ni,8]/ (d[ni,8]+d[ni,10])
  d[ni,]$ymin = 0; d[ni,]$ymax = 1-y; d[ni,]$y = (1-y) - ny*(1-y)
  for (j in 1:nrow(d)) d = set_y_coord(d, j)
  return(d)
}

set_y_coord = function(d, i){
  index = d[i,1]; parent = d[i,6]
  if (!is.na(d[i,]$y)) return(d)
  tmp = d[d[,1] == parent,]
  if (is.na(tmp[1,]$y)){
    d = set_y_coord(d, which(d[,1]==parent))
    tmp = d[d[,1]== parent,]
  }
  py = tmp[1,]$y; pymin = tmp[1,]$ymin; pymax = tmp[1,]$ymax
  f = d[i,8]/( d[i,8]+d[i,10])
  if (tmp[1,7] == index){
    d[i,]$ymin = py; d[i,]$ymax = pymax
    d[i,]$y = pymax - f*(pymax-py)
    if (d[i,5]== "TIP") d[i,]$y = (py+pymax)/2
  } else {
    d[i,]$ymin = pymin; d[i,]$ymax = py
    d[i,]$y = py - f*(py-pymin)
    if (d[i,5]== "TIP") d[i,]$y = (pymin+py)/2
  }
  return(d)
}

set_x_coords = function(d, e){
  i = which(d[,3]=="ROOT"); index = d[i,1]
  d[i,]$x = 0
  c1 = d[i,7]; c2 = d[i,9]
  ni = which(d[,1]==c1)
  tmpx = e[e[,1]==index & e[,2] == c1,3]
  if (length(tmpx) == 0){
    tmp = e[e[,1] == index,]
    tmpc1 = tmp[1,2]
    if (d[d[,1]==tmpc1,4] != "MIG") tmpc1 = tmp[2,2]
    tmpx = get_dist_to_nmig(d, e, index, tmpc1)
  }
  if (tmpx < 0) tmpx = 0
  d[ni,]$x = tmpx
  ni = which(d[,1]==c2)
  tmpx = e[e[,1]==index & e[,2] == c2,3]
  if (length(tmpx) == 0){
    tmp = e[e[,1] == index,]
    tmpc2 = tmp[2,2]
    if (d[d[,1]==tmpc2,4] != "MIG") tmpc2 = tmp[1,2]
    tmpx = get_dist_to_nmig(d, e, index, tmpc2)
  }
  if (tmpx < 0) tmpx = 0
  d[ni,]$x = tmpx
  for (j in 1:nrow(d)) d = set_x_coord(d, e, j)
  return(d)
}

set_x_coord = function(d, e, i){
  index = d[i,1]; parent = d[i,6]
  if (!is.na(d[i,]$x)) return(d)
  tmp = d[d[,1] == parent,]
  if (is.na(tmp[1,]$x)){
    d = set_x_coord(d, e, which(d[,1]==parent))
    tmp = d[d[,1]== parent,]
  }
  tmpx = e[e[,1]==parent & e[,2] == index,3]
  if (length(tmpx) == 0){
    tmp2 = e[e[,1] == parent,]
    tmpc2 = tmp2[2,2]
    if (d[d[,1]==tmpc2,4] != "MIG") tmpc2 = tmp2[1,2]
    tmpx = get_dist_to_nmig(d, e, parent, tmpc2)
  }
  if (tmpx < 0) tmpx = 0
  d[i,]$x = tmp[1,]$x + tmpx
  return(d)
}

set_mig_coords = function(d, e){
  for (j in 1:nrow(d)){
    if (d[j,4] == "MIG"){
      p = d[d[,1] == d[j,6],]; c = d[d[,1] == d[j,7],]
      tmpe = e[e[,1] == d[j,1],]
      mf = tmpe[1,6]; if (is.nan(mf)) mf = 0
      d[j,]$y = p[1,]$y + (c[1,]$y - p[1,]$y)*mf
      d[j,]$x = p[1,]$x + (c[1,]$x - p[1,]$x)*mf
    }
  }
  return(d)
}

get_dist_to_nmig = function(d, e, n1, n2){
  toreturn = e[e[,1] == n1 & e[,2] == n2,3]
  while (d[d[,1] ==n2,4] == "MIG"){
    tmp = e[e[,1] == n2 & e[,5] == "NOT_MIG",]
    toreturn = toreturn + tmp[1,3]
    n2 = tmp[1,2]
  }
  return(toreturn)
}

flip_node = function(d, n){
  i = which(d[,1] == n)
  t1 = d[i,7]; t2 = d[i,8]
  d[i,7] = d[i,9]; d[i,8] = d[i,10]
  d[i,9] = t1;    d[i,10] = t2
  return(d)
}

# ==============================================================================
# TREE PLOTTER
# ==============================================================================

plot_tree_internal = function(d, e, o = NA, cex = 1, disp = 0.005, plus = 0.005,
                              arrow = 0.05, ybar = 0.01, scale = TRUE, mbar = TRUE,
                              mse = 0.01, plotmig = TRUE, plotnames = TRUE,
                              xmin = 0, lwd = 1, font = 1,
                              style = treemix_style(), title = NULL){

  plot(d$x, d$y, axes = FALSE, ylab = "", xlab = "Drift parameter",
       xlim = c(xmin, max(d$x) + plus), pch = "")
  axis(1, col = style$frame_col, col.ticks = style$axis_col,
       lwd = 0.8, lwd.ticks = 0.8)

  has_mig <- any(e[, 5] == "MIG")
  mw <- if (has_mig) max(e[e[, 5] == "MIG", 4]) else 0
  mcols <- style$mig_palette

  for (i in 1:nrow(e)){
    col = style$tree_col
    if (e[i,5] == "MIG"){
      w = floor(e[i,4] * 200) + 50
      if (mw > 0.5) w = floor(e[i,4] * 100) + 50
      w = max(1, min(w, length(mcols)))
      col = mcols[w]
      if (is.na(col)) col = "#3a3a3a"
    }
    v1 = d[d[,1] == e[i,1], ]; v2 = d[d[,1] == e[i,2], ]
    if (e[i,5] == "MIG"){
      if (plotmig)
        arrows(v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y,
               col = col, length = arrow, lwd = lwd + 0.3, angle = 22)
    } else {
      lines(c(v1[1,]$x, v2[1,]$x), c(v1[1,]$y, v2[1,]$y),
            col = col, lwd = lwd)
    }
  }

  tmp = d[d[,5] == "TIP", ]
  if (!is.na(o[1])){
    for (i in 1:nrow(tmp)){
      tcol = o[o[,1] == tmp[i,2], 2]
      if (!length(tcol)) tcol = style$tree_col
      if (plotnames)
        text(tmp[i,]$x + disp, tmp[i,]$y, labels = tmp[i,2],
             adj = 0, cex = cex, col = tcol, font = font)
    }
  } else if (plotnames) {
    text(tmp$x + disp, tmp$y, labels = tmp[,2],
         adj = 0, cex = cex, font = font, col = style$tree_col)
  }

  if (scale){
    segments(0, ybar, mse * 10, ybar, col = style$axis_col, lwd = 1.2)
    segments(0, ybar - 0.012, 0, ybar + 0.012, col = style$axis_col, lwd = 1.2)
    segments(mse * 10, ybar - 0.012, mse * 10, ybar + 0.012,
             col = style$axis_col, lwd = 1.2)
    text(0, ybar - 0.04, "10 s.e.", adj = 0, cex = 0.75, col = style$axis_col)
  }

  if (mbar && has_mig){
    mcols_bar = style$mig_palette[50:length(style$mig_palette)]
    ymi = ybar + 0.15; yma = ybar + 0.35
    l = 0.2; w = l / length(mcols_bar)
    xma = max(d$x / 20)
    rect(rep(0, length(mcols_bar)),
         ymi + (0:(length(mcols_bar) - 1)) * w,
         rep(xma, length(mcols_bar)),
         ymi + (1:length(mcols_bar)) * w,
         col = mcols_bar, border = NA)
    rect(0, ymi, xma, yma, border = style$frame_col, lwd = 0.5)
    text(xma + disp, ymi, "0", adj = 0, cex = 0.7, col = style$axis_col)
    text(xma + disp, yma, if (mw > 0.5) "1" else "0.5",
         adj = 0, cex = 0.7, col = style$axis_col)
    text(0, yma + 0.05, "Migration", adj = 0, cex = 0.65,
         col = style$axis_col, font = 2)
    text(0, yma + 0.025, "weight", adj = 0, cex = 0.65,
         col = style$axis_col, font = 2)
  }

  if (!is.null(title))
    mtext(title, side = 3, line = 0.4, adj = 0,
          cex = 0.95, font = 2, col = style$tree_col)
}

plot_tree = function(stem, o = NA, cex = 1, disp = 0.003, plus = 0.01,
                     flip = vector(), arrow = 0.05, scale = TRUE, ybar = 0.1,
                     mbar = TRUE, plotmig = TRUE, plotnames = TRUE,
                     xmin = 0, lwd = 1, font = 1,
                     style = treemix_style(), title = NULL){
  if (!is.character(stem) || length(stem) != 1 || !nzchar(stem))
    stop("plot_tree: 'stem' must be a single non-empty string.")
  for (ext in c(".vertices.gz", ".edges.gz", ".covse.gz"))
    if (!file.exists(paste0(stem, ext)))
      stop("plot_tree: missing file ", stem, ext)

  d = read.table(gzfile(paste0(stem, ".vertices.gz")),
                 as.is = TRUE, comment.char = "", quote = "")
  e = read.table(gzfile(paste0(stem, ".edges.gz")),
                 as.is = TRUE, comment.char = "", quote = "")
  if (!is.na(o[1]))
    o = read.table(o, as.is = TRUE, comment.char = "", quote = "")
  e[,3] = e[,3] * e[,4]; e[,3] = e[,3] * e[,4]
  se = read.table(gzfile(paste0(stem, ".covse.gz")),
                  as.is = TRUE, comment.char = "", quote = "")
  m = mean(apply(se, 1, mean))
  for (i in seq_along(flip)) d = flip_node(d, flip[i])
  d$x = NA_real_; d$y = NA_real_; d$ymin = NA_real_; d$ymax = NA_real_
  d = set_y_coords(d); d = set_x_coords(d, e); d = set_mig_coords(d, e)
  plot_tree_internal(d, e, o = o, cex = cex, xmin = xmin, disp = disp,
                     plus = plus, arrow = arrow, ybar = ybar, mbar = mbar,
                     mse = m, scale = scale, plotmig = plotmig,
                     plotnames = plotnames, lwd = lwd, font = font,
                     style = style, title = title)
  invisible(list(d = d, e = e))
}

# ==============================================================================
# RESIDUAL PLOTTER
# ==============================================================================

plot_resid = function(stem, pop_order, min = -0.009, max = 0.009,
                     cex = 1, usemax = TRUE,
                     style = treemix_style(), title = NULL){
  if (!is.character(stem) || length(stem) != 1 || !nzchar(stem))
    stop("plot_resid: 'stem' must be a single non-empty string.")
  for (ext in c(".cov.gz", ".modelcov.gz", ".covse.gz"))
    if (!file.exists(paste0(stem, ext)))
      stop("plot_resid: missing file ", stem, ext)

  c = read.table(gzfile(paste0(stem, ".cov.gz")),
                 as.is = TRUE, head = TRUE, quote = "", comment.char = "",
                 check.names = FALSE)
  m = read.table(gzfile(paste0(stem, ".modelcov.gz")),
                 as.is = TRUE, head = TRUE, quote = "", comment.char = "",
                 check.names = FALSE)
  names(c) = rownames(c); names(m) = rownames(m)
  o = read.table(pop_order, as.is = TRUE, comment.char = "", quote = "")
  se = read.table(gzfile(paste0(stem, ".covse.gz")),
                  as.is = TRUE, comment.char = "", quote = "")
  mse = mean(apply(se, 1, mean))
  c = c[order(names(c)), order(names(c))]
  m = m[order(names(m)), order(names(m))]
  tmp = c - m
  toplot = data.frame(matrix(nrow = nrow(tmp), ncol = ncol(tmp)))
  for (i in 1:nrow(o)){
    for (j in 1:nrow(o)){
      if (!(o[i,1] %in% names(tmp)))
        warning("Population not in data: ", o[i,1])
      toplot[i, j] = tmp[which(names(tmp) == o[i,1]),
                         which(names(tmp) == o[j,1])]
    }
  }
  if (usemax){
    m1 = max(abs(toplot), na.rm = TRUE)
    max = m1 * 1.02; min = -(m1 * 1.02)
  }
  names(toplot) = o[,1]
  plot_resid_internal(toplot, max = max, min = min, mse = mse,
                      cex = cex, style = style, title = title)
  invisible(toplot)
}

plot_resid_internal = function(d, max = 0.009, min = -0.009, mse = NA,
                               cex = 0.5, style = treemix_style(), title = NULL){
  npop = nrow(d); w = 1 / npop
  pal = colorRampPalette(style$resid_pal)(200)
  val2col = function(v){
    if (is.na(v)) return("#f5f5f5")
    x = (v - min) / (max - min); x = pmin(pmax(x, 0), 1)
    pal[1 + floor(x * (length(pal) - 1))]
  }
  plot("NA", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  for (i in 1:npop){
    for (j in 1:i){
      v = d[i, j]; col = val2col(v)
      xmin_i = j/npop - 1/npop; xmax_i = j/npop
      ymin_i = 1 - (i/npop);    ymax_i = 1 - (i/npop) + 1/npop
      rect(xmin_i, ymin_i, xmax_i, ymax_i,
           col = col, border = "#e8e8e8", lwd = 0.4)
    }
    mtext(names(d)[i], side = 2, at = 1 - i/npop + 0.5/npop,
          las = 1, cex = cex, col = style$tree_col, line = 0.2)
    mtext(names(d)[i], side = 1, at = i/npop - 0.5/npop,
          las = 3, cex = cex, col = style$tree_col, line = 0.2)
  }
  if (!is.na(mse)){
    k = length(pal); ymi = 0.5; yma = 0.9
    w = (yma - ymi) / k
    xmi_l = 0.78; xma_l = 0.82
    lmi = round(min / mse, 1); lma = round(max / mse, 1)
    rect(rep(xmi_l, k), ymi + (0:(k-1)) * w,
         rep(xma_l, k), ymi + (1:k) * w,
         col = pal, border = NA)
    rect(xmi_l, ymi, xma_l, yma, border = style$frame_col, lwd = 0.5)
    text(xma_l + 0.01, ymi, paste(lmi, "SE"), adj = 0,
         cex = 0.75, col = style$axis_col)
    text(xma_l + 0.01, yma, paste(lma, "SE"), adj = 0,
         cex = 0.75, col = style$axis_col)
    text(xmi_l, yma + 0.04, "Residual", adj = 0,
         cex = 0.75, col = style$axis_col, font = 2)
  }
  if (!is.null(title))
    mtext(title, side = 3, line = 0.4, adj = 0,
          cex = 0.95, font = 2, col = style$tree_col)
  invisible(d)
}

# ==============================================================================
# RUN DISCOVERY  (for replicate-aware plotting)
# ------------------------------------------------------------------------------
# TreeMix's .llik file is human-readable text. Extract the final log-likelihood
# (which is typically the largest number in the file, appearing on the
# "Exiting ln(likelihood)..." line).
# ==============================================================================

read_treemix_llik <- function(llik_file) {
  if (!file.exists(llik_file)) return(NA_real_)
  txt <- paste(readLines(llik_file, warn = FALSE), collapse = " ")
  nums <- regmatches(txt,
                     gregexpr("-?[0-9]+\\.[0-9]+(?:[eE][+-]?[0-9]+)?",
                              txt, perl = TRUE))[[1]]
  if (!length(nums)) return(NA_real_)
  as.numeric(tail(nums, 1))
}

#' Scan a directory for TreeMix output sets.
#'
#' @param dir             Directory to scan.
#' @param prefix          Shared filename prefix before the m value.
#' @param m_prefix        Optional string inserted between prefix and m in
#'                        filenames. Common values:
#'                          ""    for files like  prefix.0.*   (no m prefix)
#'                          "m"   for files like  prefix.m0.*
#'                          ".m"  for files like  prefix.m0.*  when there is
#'                                an implicit dot (usually same as "m" here).
#' @param has_replicates  TRUE if filenames include `.rep<NN>`, FALSE if not.
#' @param rep_pattern     Regex for the replicate suffix. Default ".rep[0-9]+"
#'                        matches .rep01, .rep02, ..., .rep100. Override for
#'                        e.g. ".boot[0-9]+" or ".run[0-9]+".
#'
#' @return A data.frame with columns: m (int), rep (int or NA), stem (char),
#'         llik (numeric). One row per .llik file found.
find_treemix_runs <- function(dir, prefix,
                              m_prefix = "",
                              has_replicates = TRUE,
                              rep_pattern = "\\.rep[0-9]+") {
  if (has_replicates) {
    pat <- sprintf("^%s\\.%s([0-9]+)(%s)\\.llik$",
                   prefix, m_prefix, rep_pattern)
  } else {
    pat <- sprintf("^%s\\.%s([0-9]+)\\.llik$", prefix, m_prefix)
  }
  hits <- list.files(dir, pattern = pat, full.names = FALSE)
  if (!length(hits))
    return(data.frame(m = integer(), rep = integer(),
                      stem = character(), llik = numeric(),
                      stringsAsFactors = FALSE))

  ma <- regmatches(hits, regexec(pat, hits))
  out <- lapply(seq_along(hits), function(i) {
    caps <- ma[[i]]
    m_val   <- as.integer(caps[2])
    rep_val <- if (has_replicates) {
      # extract trailing integer from the rep suffix capture
      as.integer(regmatches(caps[3],
                            regexpr("[0-9]+", caps[3])))
    } else NA_integer_
    stem <- file.path(dir, sub("\\.llik$", "", hits[i]))
    data.frame(m    = m_val,
               rep  = rep_val,
               stem = stem,
               llik = read_treemix_llik(paste0(stem, ".llik")),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, out)
  df[order(df$m, df$rep), ]
}

#' Return the best-likelihood replicate for a given m, or NULL if none.
best_replicate <- function(runs_df, m) {
  sub <- runs_df[runs_df$m == m & !is.na(runs_df$llik), , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  best <- sub[which.max(sub$llik), , drop = FALSE]
  # Also verify the plot files exist (not just the .llik)
  required <- c(".vertices.gz", ".edges.gz", ".covse.gz",
                ".cov.gz", ".modelcov.gz")
  if (!all(file.exists(paste0(best$stem, required)))) return(NULL)
  best$n_reps <- nrow(sub)
  best
}

# ==============================================================================
# POP-ORDER HELPER
# ==============================================================================

generate_pop_order = function(stem, outfile = "pop_order.txt"){
  cov_file = paste0(stem, ".cov.gz")
  if (!file.exists(cov_file)) stop("Could not find file: ", cov_file)
  d = read.table(gzfile(cov_file), nrows = 1, header = TRUE, check.names = FALSE)
  pops = names(d)
  write.table(pops, file = outfile, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message("Wrote ", outfile, " with ", length(pops), " populations.")
  invisible(pops)
}

# ==============================================================================
# CAIRO EXPORT HELPER
# ==============================================================================

save_plot = function(file, expr, width = 8, height = 6, dpi = 600,
                     style = treemix_style()){
  ext = tolower(tools::file_ext(file))
  switch(ext,
    pdf = grDevices::cairo_pdf(file, width = width, height = height,
                               family = style$family),
    png = grDevices::png(file, width = width, height = height, units = "in",
                         res = dpi, type = "cairo-png", family = style$family),
    svg = grDevices::svg(file, width = width, height = height,
                         family = style$family),
    stop("Unknown extension: ", ext)
  )
  on.exit(grDevices::dev.off())
  treemix_par(style)
  force(expr)
  invisible(file)
}
