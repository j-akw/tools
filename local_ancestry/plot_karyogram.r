#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    input_dir   = NULL,
    sample      = NULL,
    output      = NULL,
    pattern     = "*.msp",
    color_0     = NULL,
    color_1     = NULL,
    color_2     = NULL,
    highlight   = NULL,
    fade_alpha  = 0.12,
    figwidth    = 260,
    figheight   = 200,
    title       = NULL,
    bar_height  = 0.35,
    gap         = 0.08,
    show_title  = TRUE,
    show_legend = TRUE,
    show_chr    = TRUE
  )

  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    nxt <- function() { i <<- i + 1; args[i] }
    if (a %in% c("-d", "--input_dir"))         { opts$input_dir   <- nxt()
    } else if (a %in% c("-s", "--sample"))     { opts$sample      <- nxt()
    } else if (a %in% c("-o", "--output"))     { opts$output      <- nxt()
    } else if (a %in% c("-p", "--pattern"))    { opts$pattern     <- nxt()
    } else if (a == "-0")                      { opts$color_0     <- nxt()
    } else if (a == "-1")                      { opts$color_1     <- nxt()
    } else if (a == "-2")                      { opts$color_2     <- nxt()
    } else if (a %in% c("-hl", "--highlight")) { opts$highlight   <- as.integer(nxt())
    } else if (a == "--fade_alpha")            { opts$fade_alpha  <- as.numeric(nxt())
    } else if (a == "--figwidth")              { opts$figwidth    <- as.numeric(nxt())
    } else if (a == "--figheight")             { opts$figheight   <- as.numeric(nxt())
    } else if (a == "--title")                 { opts$title       <- nxt()
    } else if (a == "--bar_height")            { opts$bar_height  <- as.numeric(nxt())
    } else if (a == "--gap")                   { opts$gap         <- as.numeric(nxt())
    } else if (a == "--no-title")              { opts$show_title  <- FALSE
    } else if (a == "--no-legend")             { opts$show_legend <- FALSE
    } else if (a == "--no-chr")                { opts$show_chr    <- FALSE
    } else if (a %in% c("-h", "--help"))       {
      cat("
Usage: Rscript plot_karyogram.R -d INPUT_DIR -s SAMPLE [options]

Required:
  -d, --input_dir DIR     Directory containing .msp files
  -s, --sample    ID      Sample ID (without .0/.1 suffix)

Ancestry colours:
  -0  HEX                 Colour for ancestry 0  (default: #4A90C4)
  -1  HEX                 Colour for ancestry 1  (default: #E8913A)
  -2  HEX                 Colour for ancestry 2  (default: #6DBE6D)

Highlight:
  -hl, --highlight CODE   Highlight one ancestry; fade the rest
  --fade_alpha     NUM    Opacity for faded ancestries (default: 0.12)

Layout:
  -o, --output   FILE     Output PDF (default: <sample>_karyogram.pdf)
  -p, --pattern  GLOB     Glob for .msp files (default: *.msp)
  --figwidth     NUM      Figure width in mm (default: 260)
  --figheight    NUM      Figure height in mm (default: 200)
  --bar_height   NUM      Haplotype bar height (default: 0.35)
  --gap          NUM      Gap between haplotype pair (default: 0.08)
  --title        TEXT     Custom plot title
  --no-title              Hide plot title
  --no-legend             Hide ancestry legend
  --no-chr                Hide chromosome labels

Examples:
  Rscript plot_karyogram.R -d msp_files/ -s BRH03
  Rscript plot_karyogram.R -d msp_files/ -s BRH03 -0 '#2166AC' -1 '#D6604D'
  Rscript plot_karyogram.R -d msp_files/ -s BRH03 -hl 1 --fade_alpha 0.10
")
      quit(status = 0)
    } else {
      stop(paste("Unknown argument:", a))
    }
    i <- i + 1
  }

  if (is.null(opts$input_dir)) stop("ERROR: --input_dir / -d is required.")
  if (is.null(opts$sample))    stop("ERROR: --sample / -s is required.")
  if (is.null(opts$output))    opts$output <- paste0(opts$sample, "_karyogram.pdf")

  opts
}


read_msp_files <- function(input_dir, pattern, sample_id) {
  files <- sort(Sys.glob(file.path(input_dir, pattern)))
  if (length(files) == 0) stop(paste("ERROR: No files matched", pattern, "in", input_dir))

  pop_names <- NULL
  all_segments <- list()

  for (fpath in files) {
    con <- file(fpath, "r")
    on.exit(close(con), add = TRUE)

    line1 <- readLines(con, n = 1)
    if (is.null(pop_names)) {
      parts <- strsplit(line1, ":")[[1]]
      if (length(parts) >= 2) {
        tokens <- trimws(strsplit(parts[2], "\\s+")[[1]])
        tokens <- tokens[tokens != ""]
        pop_names <- list()
        for (tok in tokens) {
          if (grepl("=", tok)) {
            kv <- strsplit(tok, "=")[[1]]
            pop_names[[as.character(as.integer(kv[2]))]] <- kv[1]
          }
        }
      }
    }

    line2 <- readLines(con, n = 1)
    line2 <- sub("^#", "", line2)
    cols  <- strsplit(line2, "\t")[[1]]

    hap0_col <- paste0(sample_id, ".0")
    hap1_col <- paste0(sample_id, ".1")
    idx0 <- match(hap0_col, cols)
    idx1 <- match(hap1_col, cols)
    if (is.na(idx0) || is.na(idx1)) {
      sample_cols <- sort(unique(sub("\\.[01]$", "", cols[7:length(cols)])))
      prefix  <- strsplit(sample_id, "_")[[1]][1]
      partial <- grep(prefix, sample_cols, value = TRUE, ignore.case = TRUE)
      msg <- paste0(
        "ERROR: Sample '", sample_id, "' not found in ", basename(fpath), ".\n",
        "  Total samples: ", length(sample_cols), "\n"
      )
      if (length(partial) > 0) {
        msg <- paste0(msg,
          "  Partial matches for '", prefix, "': ",
          paste(head(partial, 30), collapse = ", "), "\n")
      } else {
        msg <- paste0(msg,
          "  No partial matches for '", prefix, "'.\n",
          "  First 30 samples: ",
          paste(head(sample_cols, 30), collapse = ", "), "\n")
      }
      stop(msg)
    }

    close(con)
    on.exit(NULL)

    dat <- read.table(fpath, header = FALSE, sep = "\t",
                      skip = 2, stringsAsFactors = FALSE, comment.char = "")

    segs <- data.frame(
      chrom = as.character(dat[, 1]),
      spos  = as.numeric(dat[, 2]),
      epos  = as.numeric(dat[, 3]),
      anc0  = as.integer(dat[, idx0]),
      anc1  = as.integer(dat[, idx1]),
      stringsAsFactors = FALSE
    )
    all_segments[[length(all_segments) + 1]] <- segs
  }

  segments <- do.call(rbind, all_segments)

  chrom_levels <- unique(segments$chrom)
  chrom_nums   <- suppressWarnings(as.integer(chrom_levels))
  chrom_order  <- chrom_levels[order(ifelse(is.na(chrom_nums), 1000, chrom_nums),
                                     chrom_levels)]
  segments$chrom <- factor(segments$chrom, levels = chrom_order)

  list(pop_names = pop_names, segments = segments)
}


plot_karyogram <- function(pop_names, segments, opts) {

  default_pal <- c("0" = "#4A90C4", "1" = "#E8913A", "2" = "#6DBE6D",
                    "3" = "#C75A93", "4" = "#8B6DB0", "5" = "#D4A843")
  pal <- default_pal
  if (!is.null(opts$color_0)) pal["0"] <- opts$color_0
  if (!is.null(opts$color_1)) pal["1"] <- opts$color_1
  if (!is.null(opts$color_2)) pal["2"] <- opts$color_2

  bar_h <- opts$bar_height
  gap   <- opts$gap

  chroms   <- levels(segments$chrom)
  n_chroms <- length(chroms)

  chrom_y <- setNames(rev(seq_len(n_chroms)), chroms)

  rects <- bind_rows(
    segments %>%
      mutate(
        hap       = "hap0",
        ancestry  = as.character(anc0),
        ycentre   = chrom_y[as.character(chrom)],
        ymin      = ycentre + gap / 2,
        ymax      = ycentre + gap / 2 + bar_h,
        spos_mb   = spos / 1e6,
        epos_mb   = epos / 1e6
      ),
    segments %>%
      mutate(
        hap       = "hap1",
        ancestry  = as.character(anc1),
        ycentre   = chrom_y[as.character(chrom)],
        ymin      = ycentre - gap / 2 - bar_h,
        ymax      = ycentre - gap / 2,
        spos_mb   = spos / 1e6,
        epos_mb   = epos / 1e6
      )
  )

  if (!is.null(opts$highlight)) {
    rects$seg_alpha <- ifelse(rects$ancestry == as.character(opts$highlight),
                              1.0, opts$fade_alpha)
  } else {
    rects$seg_alpha <- 1.0
  }

  anc_codes <- sort(unique(rects$ancestry))
  anc_labels <- sapply(anc_codes, function(x) {
    if (!is.null(pop_names[[x]])) pop_names[[x]] else paste("Ancestry", x)
  })
  rects$ancestry <- factor(rects$ancestry, levels = anc_codes, labels = anc_labels)

  pal_named <- setNames(pal[anc_codes], anc_labels)

  y_breaks <- unname(chrom_y)
  y_labels <- paste("Chr", names(chrom_y))

  plot_title <- if (opts$show_title) {
    if (!is.null(opts$title)) opts$title else paste0("Local ancestry \u2014 ", opts$sample)
  } else {
    NULL
  }

  if (opts$show_chr) {
    y_scale <- scale_y_continuous(breaks = y_breaks, labels = y_labels,
                                  expand = expansion(add = 0.6))
    y_text  <- element_text(colour = "black", size = 7)
  } else {
    y_scale <- scale_y_continuous(breaks = NULL, expand = expansion(add = 0.6))
    y_text  <- element_blank()
  }

  legend_pos <- if (opts$show_legend) "right" else "none"

  p <- ggplot(rects) +
    geom_rect(
      aes(xmin = spos_mb, xmax = epos_mb,
          ymin = ymin, ymax = ymax,
          fill = ancestry),
      alpha = rects$seg_alpha,
      colour = NA
    ) +
    scale_fill_manual(values = pal_named, name = "Ancestry") +
    scale_x_continuous(
      name   = "Position (Mb)",
      expand = expansion(mult = c(0, 0.01))
    ) +
    y_scale +
    labs(y = NULL, title = plot_title) +
    theme_classic(base_size = 10) +
    theme(
      axis.line.x       = element_line(colour = "black", linewidth = 0.4),
      axis.line.y       = element_blank(),
      axis.ticks.x      = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.y      = element_blank(),
      axis.text.x       = element_text(colour = "black"),
      axis.text.y       = y_text,
      panel.grid        = element_blank(),
      legend.position   = legend_pos,
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(4, "mm"),
      legend.title      = element_text(size = 8, face = "bold"),
      legend.text       = element_text(size = 7),
      plot.title        = element_text(size = 11, face = "bold", hjust = 0)
    )

  ggsave(opts$output, plot = p,
         device = cairo_pdf,
         width  = opts$figwidth,
         height = opts$figheight,
         units  = "mm", dpi = 600)

  cat("Saved:", opts$output, "\n")
}


main <- function() {
  opts   <- parse_args()
  result <- read_msp_files(opts$input_dir, opts$pattern, opts$sample)
  plot_karyogram(result$pop_names, result$segments, opts)
}

main()
