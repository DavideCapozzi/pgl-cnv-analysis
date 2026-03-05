# =============================================================================
# CNVkit CNV Profile Analysis — ComplexHeatmap Edition (v5)
# =============================================================================
# Two heatmaps are produced:
#   1. Chromosome-arm-level heatmap for ARM_CHROMS, patients clustered by
#      both arms simultaneously (ward.D2 / Euclidean on log2-ratio values).
#   2. Whole-genome heatmap using per-chromosome weighted median log2-ratio.
#
# Changes from v4:
#   - Colour scale is now data-driven and asymmetric around 0 (diploid = white).
#     Breakpoints are set at (p02, p10, 0, p90, p98) of the observed log2-ratio
#     distribution, independently calibrated for deletion and gain arms.
#     This approach mirrors GISTIC2 / cBioPortal conventions and maximises
#     chromatic discrimination across the actual data range.
#   - Legend ticks are the exact colorRamp2 breakpoints — no floating-point
#     artefacts possible by construction.
#   - Column split and ward.D2 clustering unchanged from v4.
#
# Input  : *_call.cns files inside sample sub-directories of OUT_DIR.
#          Expected naming: <OUT_DIR>/<sample>/<sample>.call.cns
#          Required columns: chromosome, start, end, log2, cn
# Output : PNG files written to VIZ_DIR.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

# =============================================================================
# CONFIGURATION — edit only this block
# =============================================================================

OUT_DIR <- "D:/CNVkit/tumor/tumor_newout/dynamic_flat_reference_output"
VIZ_DIR <- file.path(OUT_DIR, "visualizations")

SAMPLE_DIR_PATTERN <- "-t"
EXCLUDED_DIRS      <- c("visualizations", "multi_sample", "nfr", "nfrmulti_sample")

ARM_CHROMS <- c("chr1", "chr3", "chr7", "chr11", "chr17", "chr22")

CENTROMERES_HG38 <- c(
  chr1  = 122026460,
  chr3  =  90772459,
  chr7  =  61053696,
  chr11 =  51078349,
  chr17 =  22813680,
  chr22 =  15979963
)

# Minimum number of samples per visible column cluster.
# Singleton clusters are merged into their nearest dendrogram neighbour.
MIN_CLUSTER_SIZE <- 2

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

if (!dir.exists(VIZ_DIR)) {
  dir.create(VIZ_DIR, recursive = TRUE)
  message("Created output directory: ", VIZ_DIR)
}

# =============================================================================
# 1.  DATA IMPORT
# =============================================================================

import_segments <- function() {
  all_subdirs <- list.dirs(OUT_DIR, full.names = FALSE, recursive = FALSE)
  
  samples <- all_subdirs[
    grepl(SAMPLE_DIR_PATTERN, all_subdirs, fixed = TRUE) &
      !(all_subdirs %in% EXCLUDED_DIRS)
  ]
  
  if (length(samples) == 0) {
    stop(
      "No sample directories found in: ", OUT_DIR,
      "\n  Pattern : '", SAMPLE_DIR_PATTERN, "'",
      "\n  Found   : ", paste(all_subdirs, collapse = ", ")
    )
  }
  
  message("Found ", length(samples), " sample(s).")
  
  map_dfr(samples, function(sample) {
    file_path <- file.path(OUT_DIR, sample, paste0(sample, ".call.cns"))
    
    if (!file.exists(file_path)) {
      message("  [SKIP] Missing: ", file_path)
      return(NULL)
    }
    
    tryCatch(
      read_tsv(file_path, show_col_types = FALSE) %>%
        filter(!is.na(cn)) %>%
        mutate(
          sample_id  = sample,
          cn         = as.numeric(cn),
          log2       = as.numeric(log2),
          chromosome = as.character(chromosome)
        ),
      error = function(e) {
        message("  [ERROR] ", file_path, ": ", e$message)
        NULL
      }
    )
  })
}

# =============================================================================
# 2.  DIAGNOSTICS
# =============================================================================

run_diagnostics <- function(segs) {
  
  message("\n", strrep("-", 60))
  message("DIAGNOSTICS")
  message(strrep("-", 60))
  
  log2_vals <- segs$log2[is.finite(segs$log2)]
  cn_vals   <- segs$cn[is.finite(segs$cn)]
  
  message("\nlog2-ratio distribution (all segments):")
  message("  min    = ", round(min(log2_vals),   3))
  message("  p02    = ", round(quantile(log2_vals, 0.02), 3))
  message("  p10    = ", round(quantile(log2_vals, 0.10), 3))
  message("  median = ", round(median(log2_vals), 3))
  message("  mean   = ", round(mean(log2_vals),   3))
  message("  p90    = ", round(quantile(log2_vals, 0.90), 3))
  message("  p98    = ", round(quantile(log2_vals, 0.98), 3))
  message("  max    = ", round(max(log2_vals),    3))
  
  message("\nCopy-number distribution (all segments):")
  print(table(segs$cn))
  
  message("\nPer-chromosome median log2-ratio:")
  chrom_summary <- segs %>%
    filter(grepl("^chr([1-9]|1[0-9]|2[012])$", chromosome)) %>%
    group_by(chromosome) %>%
    summarise(
      median_log2 = round(median(log2, na.rm = TRUE), 3),
      sd_log2     = round(sd(log2,     na.rm = TRUE), 3),
      n_segs      = n(),
      .groups     = "drop"
    ) %>%
    mutate(chrom_num = as.integer(sub("chr", "", chromosome))) %>%
    arrange(chrom_num)
  print(as.data.frame(chrom_summary))
  
  # Asymmetric scale: calibrate deletion and gain arms independently.
  # Breakpoints mirror GISTIC2 / cBioPortal conventions:
  #   (p02, p10, 0, p90, p98) — 0 = diploid is always the central anchor.
  # Values are rounded to 3 decimal places to avoid floating-point noise
  # in colorRamp2 while preserving meaningful precision for the legend.
  p02 <- round(as.numeric(quantile(log2_vals, 0.02)), 3)
  p10 <- round(as.numeric(quantile(log2_vals, 0.10)), 3)
  p90 <- round(as.numeric(quantile(log2_vals, 0.90)), 3)
  p98 <- round(as.numeric(quantile(log2_vals, 0.98)), 3)
  
  # Guard: if p10 >= 0 or p90 <= 0 (extremely skewed cohort), fall back to
  # symmetric mid-points so colorRamp2 never receives non-monotone breaks.
  if (p10 >= 0) p10 <- round(p02 / 2, 3)
  if (p90 <= 0) p90 <- round(p98 / 2, 3)
  
  scale_params <- list(p02 = p02, p10 = p10, p90 = p90, p98 = p98)
  
  message("\nAsymmetric colour scale breakpoints:")
  message("  deletion extreme  (p02) = ", p02)
  message("  deletion mid      (p10) = ", p10)
  message("  diploid anchor         =  0.000")
  message("  gain mid          (p90) = ", p90)
  message("  gain extreme      (p98) = ", p98)
  
  diag_path <- file.path(VIZ_DIR, "diagnostic_log2_distribution.png")
  png(diag_path, width = 1400, height = 600, res = 130)
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
  
  hist(log2_vals,
       breaks = 80,
       col    = "#4393C3",
       border = "white",
       main   = "log2-ratio distribution (all segments)",
       xlab   = "log2-ratio",
       ylab   = "Frequency")
  abline(v = c(p02, p10, 0, p90, p98),
         col = c("red", "orange", "black", "orange", "red"),
         lty = c(2, 3, 1, 3, 2), lwd = 2)
  legend("topright",
         legend = c("0 (diploid)", "p10 / p90", "p02 / p98"),
         col    = c("black", "orange", "red"),
         lty    = c(1, 3, 2), lwd = 2, cex = 0.8)
  
  barplot(table(cn_vals),
          col    = "#D6604D",
          border = "white",
          main   = "Integer copy-number distribution",
          xlab   = "Copy Number",
          ylab   = "Frequency")
  
  dev.off()
  message("\n  [SAVED] ", diag_path)
  message(strrep("-", 60), "\n")
  
  invisible(scale_params)
}

# =============================================================================
# 3.  MATRIX BUILDERS
# =============================================================================

build_genome_matrix <- function(segs) {
  segs %>%
    filter(grepl("^chr([1-9]|1[0-9]|2[012])$", chromosome),
           is.finite(log2)) %>%
    mutate(
      chrom_num  = as.integer(sub("chr", "", chromosome)),
      seg_length = end - start
    ) %>%
    filter(chrom_num >= 1, chrom_num <= 22) %>%
    group_by(chromosome, chrom_num, sample_id) %>%
    summarise(
      wmed_log2 = {
        o  <- order(log2)
        wt <- seg_length[o]
        cs <- cumsum(wt)
        log2[o][which(cs >= sum(wt) / 2)[1]]
      },
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols     = c(chromosome, chrom_num),
      names_from  = sample_id,
      values_from = wmed_log2,
      values_fill = 0
    ) %>%
    arrange(chrom_num) %>%
    select(-chrom_num) %>%
    column_to_rownames("chromosome") %>%
    as.matrix()
}

assign_arm <- function(start_val, end_val, centromere) {
  if (end_val   <= centromere) return("p")
  if (start_val >= centromere) return("q")
  p_len <- centromere - start_val
  q_len <- end_val    - centromere
  if (p_len >= q_len) "p" else "q"
}

build_arm_matrix <- function(segs) {
  arm_rows <- map_dfr(ARM_CHROMS, function(chrom) {
    centromere <- CENTROMERES_HG38[[chrom]]
    
    segs %>%
      filter(chromosome == chrom, is.finite(log2)) %>%
      mutate(seg_length = end - start) %>%
      rowwise() %>%
      mutate(arm = assign_arm(start, end, centromere)) %>%
      ungroup() %>%
      group_by(arm, sample_id) %>%
      summarise(
        wmean_log2 = sum(log2 * seg_length) / sum(seg_length),
        .groups    = "drop"
      ) %>%
      mutate(arm_label = paste0(chrom, arm))
  })
  
  expected_arms <- unlist(lapply(ARM_CHROMS, function(ch) paste0(ch, c("p", "q"))))
  
  mat <- arm_rows %>%
    pivot_wider(
      id_cols     = arm_label,
      names_from  = sample_id,
      values_from = wmean_log2,
      values_fill = 0
    ) %>%
    column_to_rownames("arm_label") %>%
    as.matrix()
  
  row_order <- intersect(expected_arms, rownames(mat))
  mat[row_order, , drop = FALSE]
}

# =============================================================================
# 4.  COLOUR SCALE  (asymmetric, data-driven, 0 = white anchor)
#
# Breakpoints: (p02, p10, 0, p90, p98) — calibrated independently for the
# deletion and gain arms of the observed log2-ratio distribution.
# This mirrors GISTIC2 / cBioPortal conventions and ensures full chromatic
# range is used for both arms regardless of cohort-level asymmetry.
# Legend ticks are the exact breakpoints — no floating-point artefacts.
# =============================================================================

make_log2_color_fun <- function(sp) {
  colorRamp2(
    breaks = c(sp$p02, sp$p10, 0, sp$p90, sp$p98),
    colors = c("#053061",   # deep blue  — strong deletion
               "#4393C3",   # mid blue   — shallow deletion
               "#F7F7F7",   # white      — diploid
               "#D6604D",   # mid red    — shallow gain
               "#67001F")   # deep red   — strong amplification
  )
}

legend_params <- function(sp) {
  breaks <- c(sp$p02, sp$p10, 0, sp$p90, sp$p98)
  labels <- formatC(breaks, format = "f", digits = 2)
  list(
    title         = "log2 ratio",
    at            = breaks,
    labels        = labels,
    legend_height = unit(5, "cm"),
    title_gp      = gpar(fontsize = 9, fontface = "bold"),
    labels_gp     = gpar(fontsize = 8),
    color_bar     = "continuous"
  )
}

# =============================================================================
# 5.  CLUSTERING  +  COLUMN SPLIT
# =============================================================================

ward_col_clust <- function(mat) {
  hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
}

# Determine the optimal integer k for column_split.
#
# When cluster_columns receives an hclust / dendrogram object, ComplexHeatmap
# only accepts column_split as a *single integer* — not as a factor or vector.
# The function therefore returns k as an integer scalar.
#
# Strategy: start from k_max = floor(n_samples * k_max_frac) and step down
# until every cluster produced by cutree(hc, k) contains >= min_size samples.
# This prevents visually distracting singleton columns while preserving the
# hierarchical structure drawn by the dendrogram.
derive_column_split <- function(hc, n_samples,
                                min_size   = MIN_CLUSTER_SIZE,
                                k_max_frac = 0.25) {
  
  k_upper <- max(2L, floor(n_samples * k_max_frac))
  
  for (k in seq(k_upper, 2L, by = -1L)) {
    tbl <- table(cutree(hc, k = k))
    if (all(tbl >= min_size)) {
      message(sprintf(
        "  Column split: k = %d  (sizes: %s)",
        k, paste(sort(as.integer(tbl)), collapse = ", ")
      ))
      return(k)
    }
  }
  
  message("  Column split: k = 2 (fallback — minimum k with no singletons)")
  2L
}

# =============================================================================
# 6.  SAVE HELPER
# =============================================================================

save_heatmap <- function(ht, filename, width_px = 1800, height_px = 950,
                         res = 150) {
  out_path <- file.path(VIZ_DIR, filename)
  png(out_path, width = width_px, height = height_px, res = res,
      bg = "white")
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legends          = FALSE,
       padding                = unit(c(8, 8, 8, 8), "mm"))
  dev.off()
  message("  [SAVED] ", out_path)
}

# =============================================================================
# 7.  HEATMAP 1 — Whole-genome log2-ratio
# =============================================================================

plot_genome_heatmap <- function(mat, sp) {
  col_fun   <- make_log2_color_fun(sp)
  col_clust <- ward_col_clust(mat)
  col_split <- derive_column_split(col_clust, ncol(mat))
  
  ht <- Heatmap(
    mat,
    name               = "log2 ratio",
    col                = col_fun,
    cluster_rows       = FALSE,
    cluster_columns    = col_clust,
    column_split       = col_split,
    column_gap         = unit(3, "mm"),
    cluster_column_slices = FALSE,
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    row_names_gp       = gpar(fontsize = 9),
    column_names_gp    = gpar(fontsize = 7),
    column_names_rot   = 50,
    column_title       = NULL,          # suppress per-slice titles
    column_title_gp    = gpar(fontsize = 0),
    row_title          = "Chromosome",
    row_title_gp       = gpar(fontsize = 10, fontface = "bold"),
    border             = TRUE,
    border_gp          = gpar(col = "grey30", lwd = 0.8),
    rect_gp            = gpar(col = "grey90", lwd = 0.3),
    row_names_side     = "left",
    heatmap_legend_param = legend_params(sp)
  )
  
  # Wrap in a titled layout to restore the main title above the dendrogram.
  out_path <- file.path(VIZ_DIR, "cnv_genome_heatmap.png")
  png(out_path, width = 1800, height = 950, res = 150, bg = "white")
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legends          = FALSE,
       padding                = unit(c(8, 8, 8, 8), "mm"),
       column_title = "Whole-Genome CNV Profile  \u2014  Weighted Median log2-Ratio per Chromosome",
       column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  dev.off()
  message("  [SAVED] ", out_path)
}

# =============================================================================
# 8.  HEATMAP 2 — Chromosome-arm level
# =============================================================================

plot_arm_heatmap <- function(mat, sp) {
  col_fun   <- make_log2_color_fun(sp)
  col_clust <- ward_col_clust(mat)
  col_split <- derive_column_split(col_clust, ncol(mat))
  
  chrom_of_arm  <- sub("[pq]$", "", rownames(mat))
  chrom_factor  <- factor(chrom_of_arm, levels = ARM_CHROMS)
  compact_labels <- sub("chr", "", rownames(mat))
  
  ht <- Heatmap(
    mat,
    name               = "log2 ratio",
    col                = col_fun,
    cluster_rows       = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns    = col_clust,
    column_split       = col_split,
    column_gap         = unit(3, "mm"),
    cluster_column_slices = FALSE,
    row_split          = chrom_factor,
    row_gap            = unit(3, "mm"),
    row_title_gp       = gpar(fontsize = 9, fontface = "bold"),
    row_title_rot      = 0,
    show_row_names     = TRUE,
    row_labels         = compact_labels,
    show_column_names  = TRUE,
    row_names_gp       = gpar(fontsize = 9, fontface = "bold"),
    column_names_gp    = gpar(fontsize = 7),
    column_names_rot   = 50,
    column_title       = NULL,
    column_title_gp    = gpar(fontsize = 0),
    border             = TRUE,
    border_gp          = gpar(col = "grey30", lwd = 0.8),
    rect_gp            = gpar(col = "grey90", lwd = 0.4),
    heatmap_legend_param = legend_params(sp)
  )
  
  out_path <- file.path(VIZ_DIR, "cnv_arm_heatmap.png")
  png(out_path, width = 1600, height = 900, res = 150, bg = "white")
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legends          = FALSE,
       padding                = unit(c(8, 8, 8, 8), "mm"),
       column_title = paste0(
         "Chromosome-Arm CNV Profile  \u2014  ",
         paste(ARM_CHROMS, collapse = "  \u00b7  ")
       ),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  dev.off()
  message("  [SAVED] ", out_path)
}

# =============================================================================
# 9.  MAIN
# =============================================================================

main <- function() {
  message("\n", strrep("=", 70))
  message("CNVkit — ComplexHeatmap CNV Analysis  (v4, log2-ratio)")
  message(strrep("=", 70))
  
  message("\n[1/5] Importing CNV segments ...")
  segs <- import_segments()
  
  if (is.null(segs) || nrow(segs) == 0) stop("No segment data imported.")
  
  n_samples <- length(unique(segs$sample_id))
  message("  Imported ", nrow(segs), " segments from ", n_samples, " sample(s).")
  
  message("\n[2/5] Running diagnostics ...")
  sp <- run_diagnostics(segs)
  
  message("\n[3/5] Building log2-ratio matrices ...")
  
  genome_mat <- build_genome_matrix(segs)
  message("  Genome matrix : ",
          nrow(genome_mat), " chromosomes x ", ncol(genome_mat), " samples")
  message("  log2 range    : [",
          round(min(genome_mat), 3), ", ", round(max(genome_mat), 3), "]")
  
  arm_mat <- build_arm_matrix(segs)
  message("  Arm matrix    : ",
          nrow(arm_mat), " arms x ", ncol(arm_mat), " samples")
  message("  log2 range    : [",
          round(min(arm_mat), 3), ", ", round(max(arm_mat), 3), "]")
  
  message("\n[4/5] Plotting whole-genome heatmap ...")
  plot_genome_heatmap(genome_mat, sp)
  
  message("\n[5/5] Plotting chromosome-arm heatmap ...")
  plot_arm_heatmap(arm_mat, sp)
  
  message("\n", strrep("=", 70))
  message("Done. Output in: ", VIZ_DIR)
  message(strrep("=", 70), "\n")
}

main()