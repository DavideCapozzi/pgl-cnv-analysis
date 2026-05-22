# =============================================================================
# CNVkit CNV Profile Analysis — ComplexHeatmap Edition (v6)
# =============================================================================
# Two heatmaps are produced:
#   1. Whole-genome heatmap using per-chromosome-arm weighted median log2-ratio.
#   2. Chromosome-arm-level heatmap for ARM_CHROMS, patients clustered by
#      both arms simultaneously (ward.D2 / Euclidean on log2-ratio values).
#
# Changes from v4/v5:
#   - Heatmap 1 now also stratifies the entire genome by chromosome arms (p/q)
#     using a complete set of hg38 centromere coordinates for all autosomes.
#   - Colour scale is data-driven and asymmetric around 0 (diploid = white).
#     Breakpoints are set at (p02, p10, 0, p90, p98) of the observed log2-ratio
#     distribution, independently calibrated for deletion and gain arms.
#     This approach mirrors GISTIC2 / cBioPortal conventions and maximises
#     chromatic discrimination across the actual data range.
#   - Legend ticks are the exact colorRamp2 breakpoints.
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
  library(readxl)
})

# =============================================================================
# CONFIGURATION — edit only this block
# =============================================================================

OUT_DIR <- "/mnt/d/CNVkit/tumor/tumor_newout/dynamic_flat_reference_output"
VIZ_DIR <- file.path(OUT_DIR, "visualizations")

SAMPLE_DIR_PATTERN <- "-t"
EXCLUDED_DIRS      <- c("visualizations", "multi_sample", "nfr", "nfrmulti_sample")

ARM_CHROMS <- c("chr1", "chr3", "chr7", "chr11", "chr17", "chr22")

META_FILE      <- file.path(OUT_DIR, "cohort_metadata.tsv")
CLINICAL_FILE  <- "/home/davidec/projects/pgl-cnv-analysis/results/clinca WES PGL.xlsx"

# Expanded to cover all autosomes for whole-genome arm splitting (hg38)
CENTROMERES_HG38 <- c(
  chr1  = 122026460,
  chr2  =  92326171,
  chr3  =  90772459,
  chr4  =  49660117,
  chr5  =  46405641,
  chr6  =  58830166,
  chr7  =  61053696,
  chr8  =  43838887,
  chr9  =  47367679,
  chr10 =  39258202,
  chr11 =  51078349,
  chr12 =  34856694,
  chr13 =  16000000,
  chr14 =  16000000,
  chr15 =  17000000,
  chr16 =  35335801,
  chr17 =  22813680,
  chr18 =  15460696,
  chr19 =  24681782,
  chr20 =  25963282,
  chr21 =  11288129,
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
  
  p02 <- round(as.numeric(quantile(log2_vals, 0.02)), 3)
  p10 <- round(as.numeric(quantile(log2_vals, 0.10)), 3)
  p90 <- round(as.numeric(quantile(log2_vals, 0.90)), 3)
  p98 <- round(as.numeric(quantile(log2_vals, 0.98)), 3)
  
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

assign_arm <- function(start_val, end_val, centromere) {
  if (end_val   <= centromere) return("p")
  if (start_val >= centromere) return("q")
  p_len <- centromere - start_val
  q_len <- end_val    - centromere
  if (p_len >= q_len) "p" else "q"
}

# Builds a whole-genome matrix stratified by p/q arms for all autosomes
build_genome_arm_matrix <- function(segs) {
  all_chroms <- paste0("chr", 1:22)
  
  arm_rows <- map_dfr(all_chroms, function(chrom) {
    centromere <- CENTROMERES_HG38[[chrom]]
    
    segs %>%
      filter(chromosome == chrom, is.finite(log2)) %>%
      mutate(seg_length = end - start) %>%
      rowwise() %>%
      mutate(arm = assign_arm(start, end, centromere)) %>%
      ungroup() %>%
      group_by(arm, sample_id) %>%
      summarise(
        val_log2 = {
          o  <- order(log2)
          wt <- seg_length[o]
          cs <- cumsum(wt)
          log2[o][which(cs >= sum(wt) / 2)[1]]
        },
        .groups = "drop"
      ) %>%
      mutate(arm_label = paste0(chrom, arm))
  })
  
  expected_arms <- unlist(lapply(all_chroms, function(ch) paste0(ch, c("p", "q"))))
  
  mat <- arm_rows %>%
    pivot_wider(
      id_cols     = arm_label,
      names_from  = sample_id,
      values_from = val_log2,
      values_fill = 0
    ) %>%
    column_to_rownames("arm_label") %>%
    as.matrix()
  
  row_order <- intersect(expected_arms, rownames(mat))
  mat[row_order, , drop = FALSE]
}

# Builds the targeted matrix stratified by p/q arms for ARM_CHROMS
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
# 4.  COLOUR SCALE
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
# 4.5. CLINICAL METADATA ANNOTATION
# =============================================================================

build_clinical_annotation <- function(mat, meta_path) {
  if (!file.exists(meta_path)) {
    return(NULL)
  }
  
  meta_df <- read_tsv(meta_path, show_col_types = FALSE)
  target_samples <- colnames(mat)
  meta_aligned <- meta_df[match(target_samples, meta_df$sample_id), ]
  
  meta_aligned$sex[is.na(meta_aligned$sex)] <- "Unknown"
  meta_aligned$tumor_location[is.na(meta_aligned$tumor_location)] <- "Unknown"
  
  # High-contrast clinical palette (Teal & Purple) orthogonal to RdBu heatmap
  sex_colors <- c(
    "Male"         = "#008080", # Teal
    "Female"       = "#800080", # Purple
    "Female_LOH_X" = "#BA55D3", # Medium Orchid
    "Unknown"      = "#BDBDBD", # Light Grey
    "Error"        = "#000000"  # Black
  )
  
  # Custom safe palette avoiding Teal, Purple, Blue, and Deep Red
  base_loc_colors <- c(
    "#FF7F00", # Orange
    "#4DAF4A", # Grass Green (Pure green, no blue tint)
    "#E7298A", # Magenta
    "#E6AB02", # Mustard Gold
    "#A65628", # Saddle Brown
    "#B2DF8A", # Light Green
    "#FB9A99"  # Soft Coral
  )
  
  unique_locs <- unique(meta_aligned$tumor_location)
  num_locs <- length(unique_locs)
  
  if (num_locs <= length(base_loc_colors)) {
    loc_palette <- base_loc_colors[1:num_locs]
  } else {
    loc_palette <- colorRampPalette(base_loc_colors)(num_locs)
  }
  
  loc_colors <- setNames(loc_palette, unique_locs)
  
  HeatmapAnnotation(
    Sex      = meta_aligned$sex,
    Location = meta_aligned$tumor_location,
    col      = list(Sex = sex_colors, Location = loc_colors),
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size   = unit(4, "mm"),
    show_legend        = c(TRUE, TRUE)
  )
}
# Extends build_clinical_annotation with SDHB IHC + germline variant tracks
build_extended_annotation <- function(mat, meta_path, clinical_path) {
  if (!file.exists(meta_path) || !file.exists(clinical_path)) {
    return(build_clinical_annotation(mat, meta_path))
  }

  meta_df <- read_tsv(meta_path, show_col_types = FALSE)
  clin_df <- read_excel(clinical_path)
  colnames(clin_df) <- c("sample_base", "sdhb_ihc", "germline_variant")

  # Strip "-t" suffix to join with Excel IDs (PTJ62-t → PTJ62)
  target_samples <- colnames(mat)
  base_ids       <- sub("-t$", "", target_samples)

  meta_aligned <- meta_df[match(target_samples, meta_df$sample_id), ]
  meta_aligned$sex[is.na(meta_aligned$sex)]                     <- "Unknown"
  meta_aligned$tumor_location[is.na(meta_aligned$tumor_location)] <- "Unknown"

  clin_aligned <- clin_df[match(base_ids, clin_df$sample_base), ]
  sdhb_vals    <- ifelse(is.na(clin_aligned$sdhb_ihc),
                         "Unknown",
                         as.character(clin_aligned$sdhb_ihc))

  sex_colors <- c(
    "Male"         = "#008080",
    "Female"       = "#800080",
    "Female_LOH_X" = "#BA55D3",
    "Unknown"      = "#BDBDBD",
    "Error"        = "#000000"
  )

  base_loc_colors <- c(
    "#FF7F00", "#4DAF4A", "#E7298A", "#E6AB02",
    "#A65628", "#B2DF8A", "#FB9A99"
  )
  unique_locs <- unique(meta_aligned$tumor_location)
  num_locs    <- length(unique_locs)
  loc_palette <- if (num_locs <= length(base_loc_colors)) {
    base_loc_colors[1:num_locs]
  } else {
    colorRampPalette(base_loc_colors)(num_locs)
  }
  loc_colors <- setNames(loc_palette, unique_locs)

  sdhb_colors <- c(
    "pos"     = "#D62728",  # Red — SDHB loss detected
    "neg"     = "#1F77B4",  # Blue — SDHB intact
    "Unknown" = "#BDBDBD"
  )

  HeatmapAnnotation(
    Sex      = meta_aligned$sex,
    Location = meta_aligned$tumor_location,
    SDHB_IHC = sdhb_vals,
    col = list(
      Sex      = sex_colors,
      Location = loc_colors,
      SDHB_IHC = sdhb_colors
    ),
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size   = unit(4, "mm"),
    show_legend        = c(TRUE, TRUE, TRUE)
  )
}

# =============================================================================
# 5.  CLUSTERING  +  COLUMN SPLIT
# =============================================================================

ward_col_clust <- function(mat) {
  hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
}

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
# 7.  HEATMAP 1 — Whole-genome log2-ratio (Arm Level)
# =============================================================================

plot_genome_heatmap <- function(mat, sp) {
  col_fun   <- make_log2_color_fun(sp)
  col_clust <- ward_col_clust(mat)
  col_split <- derive_column_split(col_clust, ncol(mat))
  
  top_anno <- build_clinical_annotation(mat, META_FILE)
  
  chrom_of_arm  <- sub("[pq]$", "", rownames(mat))
  chrom_factor  <- factor(chrom_of_arm, levels = paste0("chr", 1:22))
  compact_labels <- sub("chr", "", rownames(mat))
  
  ht_args <- list(
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
    row_gap            = unit(1, "mm"),
    show_row_names     = TRUE,
    row_labels         = compact_labels,
    show_column_names  = TRUE,
    row_names_gp       = gpar(fontsize = 8),
    column_names_gp    = gpar(fontsize = 7),
    column_names_rot   = 50,
    column_title       = NULL,          
    column_title_gp    = gpar(fontsize = 0),
    row_title_gp       = gpar(fontsize = 9, fontface = "bold"),
    row_title_rot      = 0,
    border             = TRUE,
    border_gp          = gpar(col = "grey30", lwd = 0.8),
    rect_gp            = gpar(col = "grey90", lwd = 0.3),
    row_names_side     = "left",
    heatmap_legend_param = legend_params(sp)
  )
  
  if (!is.null(top_anno)) {
    ht_args$top_annotation <- top_anno
  }
  
  ht <- do.call(Heatmap, ht_args)
  
  out_path <- file.path(VIZ_DIR, "cnv_genome_heatmap.png")
  png(out_path, width = 1800, height = 1100, res = 150, bg = "white")
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legends          = FALSE,
       padding                = unit(c(8, 8, 8, 8), "mm"),
       column_title = "Whole-Genome CNV Profile  \u2014  Weighted Median log2-Ratio per Chromosome Arm",
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
  
  top_anno <- build_clinical_annotation(mat, META_FILE)
  
  chrom_of_arm  <- sub("[pq]$", "", rownames(mat))
  chrom_factor  <- factor(chrom_of_arm, levels = ARM_CHROMS)
  compact_labels <- sub("chr", "", rownames(mat))
  
  ht_args <- list(
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
  
  if (!is.null(top_anno)) {
    ht_args$top_annotation <- top_anno
  }
  
  ht <- do.call(Heatmap, ht_args)
  
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
# 9a. HEATMAP 3 — Whole-genome + clinical annotation (SDHB IHC + Germline)
# =============================================================================

plot_genome_heatmap_clinical <- function(mat, sp) {
  col_fun   <- make_log2_color_fun(sp)
  col_clust <- ward_col_clust(mat)
  col_split <- derive_column_split(col_clust, ncol(mat))

  top_anno <- build_extended_annotation(mat, META_FILE, CLINICAL_FILE)

  chrom_of_arm   <- sub("[pq]$", "", rownames(mat))
  chrom_factor   <- factor(chrom_of_arm, levels = paste0("chr", 1:22))
  compact_labels <- sub("chr", "", rownames(mat))

  ht_args <- list(
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
    row_gap            = unit(1, "mm"),
    show_row_names     = TRUE,
    row_labels         = compact_labels,
    show_column_names  = TRUE,
    row_names_gp       = gpar(fontsize = 8),
    column_names_gp    = gpar(fontsize = 7),
    column_names_rot   = 50,
    column_title       = NULL,
    column_title_gp    = gpar(fontsize = 0),
    row_title_gp       = gpar(fontsize = 9, fontface = "bold"),
    row_title_rot      = 0,
    border             = TRUE,
    border_gp          = gpar(col = "grey30", lwd = 0.8),
    rect_gp            = gpar(col = "grey90", lwd = 0.3),
    row_names_side     = "left",
    heatmap_legend_param = legend_params(sp)
  )

  if (!is.null(top_anno)) ht_args$top_annotation <- top_anno

  ht <- do.call(Heatmap, ht_args)

  out_path <- file.path(VIZ_DIR, "cnv_genome_heatmap_clinical.png")
  png(out_path, width = 1800, height = 1200, res = 150, bg = "white")
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legends          = FALSE,
       padding                = unit(c(8, 8, 8, 8), "mm"),
       column_title    = "Whole-Genome CNV Profile — SDHB IHC",
       column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  dev.off()
  message("  [SAVED] ", out_path)
}

# =============================================================================
# 9b. HEATMAP 4 — Chromosome-arm + clinical annotation (SDHB IHC + Germline)
# =============================================================================

plot_arm_heatmap_clinical <- function(mat, sp) {
  col_fun   <- make_log2_color_fun(sp)
  col_clust <- ward_col_clust(mat)
  col_split <- derive_column_split(col_clust, ncol(mat))

  top_anno <- build_extended_annotation(mat, META_FILE, CLINICAL_FILE)

  chrom_of_arm   <- sub("[pq]$", "", rownames(mat))
  chrom_factor   <- factor(chrom_of_arm, levels = ARM_CHROMS)
  compact_labels <- sub("chr", "", rownames(mat))

  ht_args <- list(
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

  if (!is.null(top_anno)) ht_args$top_annotation <- top_anno

  ht <- do.call(Heatmap, ht_args)

  out_path <- file.path(VIZ_DIR, "cnv_arm_heatmap_clinical.png")
  png(out_path, width = 1600, height = 1000, res = 150, bg = "white")
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legends          = FALSE,
       padding                = unit(c(8, 8, 8, 8), "mm"),
       column_title    = paste0(
         "Chromosome-Arm CNV Profile — SDHB IHC  ·  ",
         paste(ARM_CHROMS, collapse = "  ·  ")
       ),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  dev.off()
  message("  [SAVED] ", out_path)
}

# =============================================================================
# 9c. HEATMAP 5 — Arm level, columns split by combined clinical group
#     Single annotation bar: germline variant (SDHD/SDHB) takes priority,
#     then SDHB IHC pos/neg, then Unknown for samples with no data.
# =============================================================================

plot_arm_heatmap_clinical_grouped <- function(mat, sp) {
  if (!file.exists(META_FILE) || !file.exists(CLINICAL_FILE)) {
    message("  [SKIP] Missing META_FILE or CLINICAL_FILE for grouped heatmap")
    return(invisible(NULL))
  }

  clin_df <- read_excel(CLINICAL_FILE)
  colnames(clin_df) <- c("sample_base", "sdhb_ihc", "germline_variant")

  target_samples <- colnames(mat)
  base_ids       <- sub("-t$", "", target_samples)
  clin_aligned   <- clin_df[match(base_ids, clin_df$sample_base), ]

  # Combined label based on SDHB IHC only
  clinical_group <- sapply(clin_aligned$sdhb_ihc, function(sdhb) {
    sdhb <- if (is.na(sdhb)) "Unknown" else as.character(sdhb)
    if (sdhb == "pos") "SDHB-IHC pos"
    else if (sdhb == "neg") "SDHB-IHC neg"
    else "Unknown"
  })

  all_group_colors <- c(
    "SDHB-IHC pos" = "#D62728",
    "SDHB-IHC neg" = "#1F77B4",
    "Unknown"      = "#BDBDBD"
  )
  group_colors <- all_group_colors[names(all_group_colors) %in% unique(clinical_group)]

  top_anno <- HeatmapAnnotation(
    Clinical_Group   = clinical_group,
    col              = list(Clinical_Group = group_colors),
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size   = unit(5, "mm"),
    show_legend        = TRUE
  )

  col_fun      <- make_log2_color_fun(sp)
  group_factor <- factor(clinical_group,
                         levels = c("SDHB-IHC pos", "SDHB-IHC neg", "Unknown"))

  chrom_of_arm   <- sub("[pq]$", "", rownames(mat))
  chrom_factor   <- factor(chrom_of_arm, levels = ARM_CHROMS)
  compact_labels <- sub("chr", "", rownames(mat))

  ht <- Heatmap(
    mat,
    name                  = "log2 ratio",
    col                   = col_fun,
    cluster_rows          = FALSE,
    cluster_row_slices    = FALSE,
    cluster_columns       = TRUE,
    clustering_method_columns = "ward.D2",
    clustering_distance_columns = "euclidean",
    column_split          = group_factor,
    column_gap            = unit(4, "mm"),
    cluster_column_slices = FALSE,
    row_split             = chrom_factor,
    row_gap               = unit(3, "mm"),
    row_title_gp          = gpar(fontsize = 9, fontface = "bold"),
    row_title_rot         = 0,
    show_row_names        = TRUE,
    row_labels            = compact_labels,
    show_column_names     = TRUE,
    row_names_gp          = gpar(fontsize = 9, fontface = "bold"),
    column_names_gp       = gpar(fontsize = 7),
    column_names_rot      = 50,
    column_title_gp       = gpar(fontsize = 9, fontface = "bold"),
    border                = TRUE,
    border_gp             = gpar(col = "grey30", lwd = 0.8),
    rect_gp               = gpar(col = "grey90", lwd = 0.4),
    top_annotation        = top_anno,
    heatmap_legend_param  = legend_params(sp)
  )

  out_path <- file.path(VIZ_DIR, "cnv_arm_heatmap_clinical_grouped.png")
  png(out_path, width = 1800, height = 1000, res = 150, bg = "white")
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legends          = FALSE,
       padding                = unit(c(8, 8, 8, 8), "mm"),
       column_title    = paste0(
         "Chromosome-Arm CNV Profile — Grouped by Clinical Status  ·  ",
         paste(ARM_CHROMS, collapse = "  ·  ")
       ),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  dev.off()
  message("  [SAVED] ", out_path)
}

# =============================================================================
# 10.  ORCHESTRATOR / MAIN EXECUTION
# =============================================================================

main <- function() {
  message("=== Starting CNV Profile Analysis Pipeline ===")
  
  # 1. Import Data
  segs <- import_segments()
  if (is.null(segs) || nrow(segs) == 0) {
    stop("Execution halted: No valid segments imported.")
  }
  
  # 2. Diagnostics & Color Mapping
  sp <- run_diagnostics(segs)
  
  # 3. Whole Genome Heatmap
  message("Building whole-genome matrix...")
  mat_genome <- build_genome_arm_matrix(segs)
  message("Plotting whole-genome heatmap...")
  plot_genome_heatmap(mat_genome, sp)
  
  message("\nBuilding targeted arm matrix...")
  mat_arm <- build_arm_matrix(segs)
  message("Plotting targeted arm heatmap...")
  plot_arm_heatmap(mat_arm, sp)

  message("\nPlotting whole-genome heatmap with clinical annotation...")
  plot_genome_heatmap_clinical(mat_genome, sp)

  message("Plotting targeted arm heatmap with clinical annotation...")
  plot_arm_heatmap_clinical(mat_arm, sp)

  message("Plotting arm heatmap grouped by clinical status...")
  plot_arm_heatmap_clinical_grouped(mat_arm, sp)

  message("\n=== Pipeline Execution Completed Successfully ===")
}

# Execute the script
main()