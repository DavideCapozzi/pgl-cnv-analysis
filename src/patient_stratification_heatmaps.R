# =============================================================================
# CNVkit CNV Profile Analysis — ComplexHeatmap Edition (v3)
# =============================================================================
# Two heatmaps are produced:
#   1. Chromosome-arm-level heatmap for ARM_CHROMS, patients clustered by
#      both arms simultaneously (ward.D2 / Manhattan on log2-ratio values).
#   2. Whole-genome heatmap using per-chromosome median log2-ratio.
#
# A diagnostic section runs first to reveal data distribution and guide
# colour-scale decisions before committing to the final plots.
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

SAMPLE_DIR_PATTERN <- "-t"          # fixed string present in every sample folder
EXCLUDED_DIRS      <- c("visualizations", "multi_sample", "nfr", "nfrmulti_sample")

# Chromosomes for the arm-level heatmap
ARM_CHROMS <- c("chr1", "chr3", "chr7", "chr11", "chr17", "chr22")

# hg38 centromere midpoints (bp) — UCSC cytoBand table, hg38 release 2023.
CENTROMERES_HG38 <- c(
  chr1  = 122026460,
  chr3  =  90772459,
  chr7  =  61053696,
  chr11 =  51078349,
  chr17 =  22813680,
  chr22 =  15979963
)

# Colour-scale saturation for log2-ratio plots.
# Values beyond ±LOG2_SAT are clipped to the extreme colours.
# Will be auto-calibrated from data in run_diagnostics() — override here if needed.
LOG2_SAT <- NULL    # NULL = auto-detect from 2nd/98th percentile of the data

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
#     Prints a full statistical summary of log2-ratio and CN distributions.
#     Saves a side-by-side diagnostic PNG so you can judge the colour scale.
# =============================================================================

run_diagnostics <- function(segs) {
  
  message("\n", strrep("-", 60))
  message("DIAGNOSTICS")
  message(strrep("-", 60))
  
  # --- Global log2 distribution -------------------------------------------
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
  
  # --- Per-chromosome median log2 summary ---------------------------------
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
  
  # --- Recommend LOG2_SAT --------------------------------------------------
  p02 <- quantile(log2_vals, 0.02)
  p98 <- quantile(log2_vals, 0.98)
  recommended_sat <- round(max(abs(p02), abs(p98)), 2)
  message("\nRecommended LOG2_SAT (max of |p02|, |p98|): ", recommended_sat)
  
  # Set global LOG2_SAT if still NULL
  if (is.null(LOG2_SAT)) {
    LOG2_SAT <<- recommended_sat
    message("  -> LOG2_SAT set automatically to ", LOG2_SAT)
  } else {
    message("  -> LOG2_SAT already set to ", LOG2_SAT, " (override active)")
  }
  
  # --- Save diagnostic histogram PNG --------------------------------------
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
  abline(v = c(-LOG2_SAT, 0, LOG2_SAT),
         col = c("red", "black", "red"),
         lty = c(2, 1, 2), lwd = 2)
  legend("topright",
         legend = c("0 (diploid)", paste0("+/-", LOG2_SAT, " sat.")),
         col    = c("black", "red"),
         lty    = c(1, 2), lwd = 2, cex = 0.8)
  
  barplot(table(cn_vals),
          col    = "#D6604D",
          border = "white",
          main   = "Integer copy-number distribution",
          xlab   = "Copy Number",
          ylab   = "Frequency")
  
  dev.off()
  message("\n  [SAVED] ", diag_path)
  message(strrep("-", 60), "\n")
  
  invisible(LOG2_SAT)
}

# =============================================================================
# 3.  MATRIX BUILDERS  (log2-ratio based)
# =============================================================================

# Autosomes chr1-chr22: weighted median log2 per chromosome per sample.
# Weighting by segment length gives more reliable estimates than simple median.
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
      # Weighted median: segments weighted by their bp length
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
      values_fill = 0          # 0 = diploid in log2-ratio space
    ) %>%
    arrange(chrom_num) %>%
    select(-chrom_num) %>%
    column_to_rownames("chromosome") %>%
    as.matrix()
}

# Arm-level: weighted mean log2 per arm per sample.
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
# 4.  COLOUR SCALE  (log2-ratio, symmetric around 0)
# =============================================================================

make_log2_color_fun <- function(sat = 1.0) {
  # sat = saturation point: values beyond ±sat clipped to extreme colour
  colorRamp2(
    breaks = c(-sat, -sat / 2, 0, sat / 2, sat),
    colors = c("#053061",   # deep blue  — strong deletion
               "#4393C3",   # mid blue   — shallow deletion
               "#F7F7F7",   # white      — diploid (log2 = 0)
               "#D6604D",   # mid red    — shallow gain
               "#67001F")   # deep red   — strong amplification
  )
}

# =============================================================================
# 5.  CLUSTERING
# =============================================================================

ward_col_clust <- function(mat) {
  hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
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

plot_genome_heatmap <- function(mat, sat) {
  col_fun   <- make_log2_color_fun(sat)
  col_clust <- ward_col_clust(mat)
  
  legend_breaks <- round(seq(-sat, sat, length.out = 5), 2)
  
  ht <- Heatmap(
    mat,
    name               = "log2 ratio",
    col                = col_fun,
    cluster_rows       = FALSE,
    cluster_columns    = col_clust,
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    row_names_gp       = gpar(fontsize = 9),
    column_names_gp    = gpar(fontsize = 7),
    column_names_rot   = 50,
    column_title       = "Whole-Genome CNV Profile  \u2014  Weighted Median log2-Ratio per Chromosome",
    column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    row_title          = "Chromosome",
    row_title_gp       = gpar(fontsize = 10, fontface = "bold"),
    border             = TRUE,
    border_gp          = gpar(col = "grey30", lwd = 0.8),
    rect_gp            = gpar(col = "grey90", lwd = 0.3),
    row_names_side     = "left",
    heatmap_legend_param = list(
      title          = "log2 ratio",
      at             = legend_breaks,
      labels         = as.character(legend_breaks),
      legend_height  = unit(5, "cm"),
      title_gp       = gpar(fontsize = 9, fontface = "bold"),
      labels_gp      = gpar(fontsize = 8),
      color_bar      = "continuous"
    )
  )
  
  save_heatmap(ht, "cnv_genome_heatmap.png",
               width_px = 1800, height_px = 950)
}

# =============================================================================
# 8.  HEATMAP 2 — Chromosome-arm level (row_split replaces row_gap)
# =============================================================================

plot_arm_heatmap <- function(mat, sat) {
  col_fun   <- make_log2_color_fun(sat)
  col_clust <- ward_col_clust(mat)
  
  # Chromosome grouping for row_split (drives the inter-chromosome gap)
  # Factor with levels in genomic order guarantees correct row ordering.
  chrom_of_arm  <- sub("[pq]$", "", rownames(mat))
  chrom_factor  <- factor(chrom_of_arm, levels = ARM_CHROMS)
  
  # Compact axis labels: "chr1p" -> "1p"
  compact_labels <- sub("chr", "", rownames(mat))
  
  
  legend_breaks <- round(seq(-sat, sat, length.out = 5), 2)
  
  ht <- Heatmap(
    mat,
    name               = "log2 ratio",
    col                = col_fun,
    cluster_rows       = FALSE,        # order governed by row_split factor
    cluster_row_slices = FALSE,        # keep chromosome order fixed
    cluster_columns    = col_clust,
    row_split          = chrom_factor, # one slice per chromosome → proper gap
    row_gap            = unit(3, "mm"),# single scalar: gap between all slices
    row_title_gp       = gpar(fontsize = 9, fontface = "bold"),
    row_title_rot      = 0,
    show_row_names     = TRUE,
    row_labels         = compact_labels,
    show_column_names  = TRUE,
    row_names_gp       = gpar(fontsize = 9, fontface = "bold"),
    column_names_gp    = gpar(fontsize = 7),
    column_names_rot   = 50,
    column_title       = paste0(
      "Chromosome-Arm CNV Profile  \u2014  ",
      paste(ARM_CHROMS, collapse = "  \u00b7  ")
    ),
    column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    border             = TRUE,
    border_gp          = gpar(col = "grey30", lwd = 0.8),
    rect_gp            = gpar(col = "grey90", lwd = 0.4),
    heatmap_legend_param = list(
      title          = "log2 ratio",
      at             = legend_breaks,
      labels         = as.character(legend_breaks),
      legend_height  = unit(5, "cm"),
      title_gp       = gpar(fontsize = 9, fontface = "bold"),
      labels_gp      = gpar(fontsize = 8),
      color_bar      = "continuous"
    )
  )
  
  save_heatmap(ht, "cnv_arm_heatmap.png",
               width_px = 1600, height_px = 900)
}

# =============================================================================
# 9.  MAIN
# =============================================================================

main <- function() {
  message("\n", strrep("=", 70))
  message("CNVkit — ComplexHeatmap CNV Analysis  (v3, log2-ratio)")
  message(strrep("=", 70))
  
  # 1. Import ----------------------------------------------------------------
  message("\n[1/5] Importing CNV segments ...")
  segs <- import_segments()
  
  if (is.null(segs) || nrow(segs) == 0) {
    stop("No segment data imported.")
  }
  
  n_samples <- length(unique(segs$sample_id))
  message("  Imported ", nrow(segs), " segments from ", n_samples, " sample(s).")
  
  # 2. Diagnostics -----------------------------------------------------------
  message("\n[2/5] Running diagnostics ...")
  sat <- run_diagnostics(segs)   # also sets LOG2_SAT globally
  
  # 3. Build matrices --------------------------------------------------------
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
  
  # 4. Whole-genome heatmap --------------------------------------------------
  message("\n[4/5] Plotting whole-genome heatmap ...")
  plot_genome_heatmap(genome_mat, sat)
  
  # 5. Arm-level heatmap -----------------------------------------------------
  message("\n[5/5] Plotting chromosome-arm heatmap ...")
  plot_arm_heatmap(arm_mat, sat)
  
  message("\n", strrep("=", 70))
  message("Done. Output in: ", VIZ_DIR)
  message(strrep("=", 70), "\n")
}

main()