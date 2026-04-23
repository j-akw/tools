# ==============================================================================
# TreeMix plotting driver  —  configurable, replicate-aware
# ------------------------------------------------------------------------------
# Produces publication-ready PDFs for:
#   * A faceted panel showing the best-likelihood tree for each m
#   * A faceted panel showing the residual covariance for each m
#   * A main-figure tree for your chosen m (best replicate by log-likelihood)
#   * A main-figure residual matrix for that same replicate
#
# All user-configurable options live in the CONFIG block below. Nothing else
# in this script needs to be edited for a new dataset.
# ==============================================================================

# -------------------- 1. Load the function library ----------------------------
# Edit this path if the functions file lives somewhere else.
FUNCS_FILE <- "treemix_plotting_funcs.R"
source(FUNCS_FILE)

# ============================== CONFIG ========================================
# Everything below is meant to be edited per-project. The defaults describe a
# replicated run with files named   <prefix>.<m>.rep<NN>.<ext>   in runs_dir.
# ==============================================================================

CONFIG <- list(

  # ---- paths -----------------------------------------------------------------
  runs_dir    = "optm_runs",           # directory containing TreeMix outputs
  plots_dir   = "PLOTS",               # output directory (created if missing)
  pop_order   = "pop_order.txt",       # written by this script, used by resid plot

  # ---- filename convention ---------------------------------------------------
  # Files are assumed to be named:   <prefix>.<m_prefix><m>[<rep_suffix>].<ext>
  # Examples:
  #    prefix = "dml_data", m_prefix = "",  has_replicates = TRUE
  #        -> matches  dml_data.0.rep01.cov.gz, dml_data.1.rep02.cov.gz, ...
  #    prefix = "out",      m_prefix = "m", has_replicates = TRUE
  #        -> matches  out.m0.rep01.cov.gz, out.m1.rep02.cov.gz, ...
  #    prefix = "run",      m_prefix = "",  has_replicates = FALSE
  #        -> matches  run.0.cov.gz, run.1.cov.gz, ...
  prefix          = "dml_data",
  m_prefix        = "",
  has_replicates  = TRUE,
  rep_pattern     = "\\.rep[0-9]+",    # regex for replicate suffix

  # ---- analysis --------------------------------------------------------------
  m_range   = 0:5,        # which migration-edge values to include in facets
  best_m    = 1,          # the "chosen m" for the main figures
                          # (set based on OptM Deltam peak or variance criterion)

  # ---- figure dimensions (inches) --------------------------------------------
  facet_nrow      = 2,
  facet_ncol      = 3,
  facet_width     = 16,
  facet_tree_h    = 9,
  facet_resid_h   = 10,
  main_tree_wh    = c(9, 6.5),
  main_resid_wh   = c(7, 6.5),

  # ---- aesthetics ------------------------------------------------------------
  # Pass NULL to use the built-in Zissou1 defaults, or override with any colour
  # vector. Examples:
  #   mig_palette = viridisLite::viridis(250)
  #   resid_palette = rev(RColorBrewer::brewer.pal(11, "RdBu"))
  mig_palette     = NULL,
  resid_palette   = NULL,
  font_family     = "sans",          # try "Helvetica" or "Arial" if installed

  # ---- text size multipliers -------------------------------------------------
  facet_tree_cex   = 0.8,
  facet_resid_cex  = 0.55,
  main_tree_cex    = 1.0,
  main_resid_cex   = 0.9,

  # ---- output format ---------------------------------------------------------
  # PDF is vector (cairo_pdf); PNG is raster at the given DPI; SVG is editable.
  output_ext   = "pdf",                # "pdf", "png", or "svg"
  png_dpi      = 600
)

# ============================== END CONFIG ====================================

# -------------------- 2. Derived setup ---------------------------------------

dir.create(CONFIG$plots_dir, showWarnings = FALSE, recursive = TRUE)

style <- treemix_style(mig_palette   = CONFIG$mig_palette,
                       resid_palette = CONFIG$resid_palette,
                       family        = CONFIG$font_family)

# Helper: build an output file path with the chosen extension.
out <- function(name) file.path(CONFIG$plots_dir,
                                paste0(name, ".", CONFIG$output_ext))

# -------------------- 3. Discover runs ---------------------------------------

cat("Scanning ", CONFIG$runs_dir, " ...\n", sep = "")
runs <- find_treemix_runs(CONFIG$runs_dir,
                          prefix         = CONFIG$prefix,
                          m_prefix       = CONFIG$m_prefix,
                          has_replicates = CONFIG$has_replicates,
                          rep_pattern    = CONFIG$rep_pattern)

if (!nrow(runs)) {
  stop("No TreeMix runs found in ", CONFIG$runs_dir,
       " matching prefix='", CONFIG$prefix,
       "', m_prefix='", CONFIG$m_prefix,
       "', has_replicates=", CONFIG$has_replicates,
       ". Check CONFIG.")
}

# Summary table
cat("\nFound ", nrow(runs), " runs.\n", sep = "")
cat("Summary (best replicate per m):\n")
for (m in sort(unique(runs$m))) {
  b <- best_replicate(runs, m)
  if (is.null(b)) {
    cat(sprintf("  m=%d : incomplete run set\n", m))
  } else {
    cat(sprintf("  m=%d : %s (llik=%.2f, %d replicate%s)\n",
                m, basename(b$stem), b$llik, b$n_reps,
                if (b$n_reps > 1) "s" else ""))
  }
}
cat("\n")

# -------------------- 4. Pop-order file (from any m=min run) -----------------

min_m <- min(runs$m)
first_ok <- best_replicate(runs, min_m)
if (is.null(first_ok))
  stop("Cannot build pop_order: no complete run set at m=", min_m)

generate_pop_order(first_ok$stem, CONFIG$pop_order)

# -------------------- 5. Faceted tree panel ----------------------------------

save_plot(out("TreeMix_Trees_Facet"),
          width = CONFIG$facet_width, height = CONFIG$facet_tree_h,
          dpi = CONFIG$png_dpi, style = style, {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  treemix_par(style, mar = c(3, 1, 2.2, 2))
  par(mfrow = c(CONFIG$facet_nrow, CONFIG$facet_ncol),
      oma   = c(0, 0, 2.5, 0))

  for (m in CONFIG$m_range) {
    best <- best_replicate(runs, m)
    if (is.null(best)) {
      plot.new(); title(main = paste("Missing: m =", m), col.main = "#999999")
      next
    }
    label <- if (best$n_reps > 1)
      sprintf("m = %d (best of %d)", m, best$n_reps)
    else
      sprintf("m = %d", m)
    tryCatch(
      plot_tree(best$stem, cex = CONFIG$facet_tree_cex, style = style,
                title = label),
      error = function(e) {
        plot.new(); title(main = paste("Error at m =", m), col.main = "#cc3333")
        message("Error for m=", m, ": ", conditionMessage(e))
      })
  }
  mtext("TreeMix: best-likelihood tree per m", outer = TRUE,
        cex = 1.1, font = 2, col = style$tree_col, line = 0.8)
})

# -------------------- 6. Faceted residual panel ------------------------------

save_plot(out("TreeMix_Residuals_Facet"),
          width = CONFIG$facet_width, height = CONFIG$facet_resid_h,
          dpi = CONFIG$png_dpi, style = style, {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  treemix_par(style, mar = c(4.5, 4.5, 2.2, 5))
  par(mfrow = c(CONFIG$facet_nrow, CONFIG$facet_ncol),
      oma   = c(0, 0, 2.5, 0))

  for (m in CONFIG$m_range) {
    best <- best_replicate(runs, m)
    if (is.null(best)) {
      plot.new(); title(main = paste("Missing: m =", m), col.main = "#999999")
      next
    }
    tryCatch(
      plot_resid(best$stem, CONFIG$pop_order, cex = CONFIG$facet_resid_cex,
                 style = style, title = sprintf("m = %d", m)),
      error = function(e) {
        plot.new(); title(main = paste("Error at m =", m), col.main = "#cc3333")
        message("Error for m=", m, ": ", conditionMessage(e))
      })
  }
  mtext("Residual covariance (data − model), in SE units", outer = TRUE,
        cex = 1.1, font = 2, col = style$tree_col, line = 0.8)
})

# -------------------- 7. Main-figure tree & residual for best_m --------------

main_best <- best_replicate(runs, CONFIG$best_m)
if (is.null(main_best))
  stop("No valid replicates found for best_m = ", CONFIG$best_m)

cat(sprintf("Main figure: m=%d, replicate = %s (llik=%.2f, from %d rep%s)\n",
            CONFIG$best_m, basename(main_best$stem),
            main_best$llik, main_best$n_reps,
            if (main_best$n_reps > 1) "s" else ""))

save_plot(out(sprintf("TreeMix_Tree_m%d_best", CONFIG$best_m)),
          width = CONFIG$main_tree_wh[1], height = CONFIG$main_tree_wh[2],
          dpi = CONFIG$png_dpi, style = style, {
  label <- if (main_best$n_reps > 1)
    sprintf("TreeMix — m = %d (best of %d runs)",
            CONFIG$best_m, main_best$n_reps)
  else
    sprintf("TreeMix — m = %d", CONFIG$best_m)
  plot_tree(main_best$stem,
            cex = CONFIG$main_tree_cex, lwd = 1.8, arrow = 0.08,
            style = style, title = label)
})

save_plot(out(sprintf("TreeMix_Residuals_m%d_best", CONFIG$best_m)),
          width = CONFIG$main_resid_wh[1], height = CONFIG$main_resid_wh[2],
          dpi = CONFIG$png_dpi, style = style, {
  plot_resid(main_best$stem, CONFIG$pop_order,
             cex = CONFIG$main_resid_cex, style = style,
             title = sprintf("Residuals — m = %d", CONFIG$best_m))
})

message("\nDone. Outputs in: ", CONFIG$plots_dir)
