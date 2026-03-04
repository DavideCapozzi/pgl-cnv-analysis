#!/bin/bash
set -euo pipefail

# =====================================================================
# CNVkit Module: Visualization Suite (v4.2)
# Description: Generates Cohort Heatmaps and Single Sample plots.
# Resolves WSL permission issues and purely filters autosomes for heatmap.
# =====================================================================

# --- CONFIGURATION PATHS ---
OUT_DIR="/mnt/d/CNVkit/tumor/tumor_newout/dynamic_flat_reference_output"
VIZ_DIR="${OUT_DIR}/visualizations"
TMP_CNS_DIR="${VIZ_DIR}/tmp_cns_copies"

mkdir -p "$VIZ_DIR"

# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================

log_info() { echo -e "[INFO] $(date '+%H:%M:%S') - $1"; }
log_error() { echo -e "[ERROR] $(date '+%H:%M:%S') - $1" >&2; }

cleanup() {
    [[ -d "$TMP_CNS_DIR" ]] && { 
        log_info "Cleaning up temporary copy directory..."
        rm -rf "$TMP_CNS_DIR"
    } || true
}
trap cleanup EXIT

# =====================================================================
# COHORT HEATMAPS
# =====================================================================

generate_cohort_heatmaps() {
    log_info "Preparing files for Cohort Heatmaps..."
    mkdir -p "$TMP_CNS_DIR"
    
    local cns_list=$(find "$OUT_DIR" -name "*.call.cns" | sort)
    
    [[ -z "$cns_list" ]] && {
        log_error "No .call.cns files found in $OUT_DIR. Cannot generate heatmaps."
        exit 1
    }

    local clean_cns_list=""
    for file in $cns_list; do
        local base_name=$(basename "$file")
        local clean_name="${base_name/.call.cns/.cns}"
        local copy_path="${TMP_CNS_DIR}/${clean_name}"
        
        # Safely filter autosomes (chr1-chr22) via awk and keep the header (NR==1)
        # This replaces cp and perfectly drops non-autosomes without CNVkit -c errors
        awk 'NR==1 || $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "$file" > "$copy_path"
        
        clean_cns_list="$clean_cns_list $copy_path"
    done

    local count=$(echo "$cns_list" | wc -l)
    log_info "Found $count samples. Generating Global Heatmap (Autosomes Only)..."

    # 1. Global Heatmap (Autosomes Only)
    # The -c flag is removed because the files are already pre-filtered for autosomes
    cnvkit.py heatmap $clean_cns_list -d --output "${VIZ_DIR}/cohort_heatmap_global_autosomes.pdf" \
        && log_info "Global autosomal heatmap created." \
        || log_error "Failed to create global heatmap."

    # 2. PGL Specific Chromosomes
    log_info "Generating targeted heatmaps for PGL loci..."
    for chr in chr1 chr3 chr7 chr11 chr17 chr22; do
         cnvkit.py heatmap $clean_cns_list -d -c "$chr" --output "${VIZ_DIR}/cohort_heatmap_${chr}.pdf" \
            && log_info "Heatmap for $chr created." \
            || log_error "Failed to create heatmap for $chr."
    done
}

# =====================================================================
# SINGLE SAMPLE PLOTS
# =====================================================================

generate_single_plots() {
    log_info "Generating Single Sample Plots (Scatter & Diagram)..."

    find "$OUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort | while read sample_dir; do
        
        local sample_name=$(basename "$sample_dir")
        local cnr_file="${sample_dir}/${sample_name}.cnr"
        local call_file="${sample_dir}/${sample_name}.call.cns"
        local cns_file="${sample_dir}/${sample_name}.cns"
        
        local seg_file="$cns_file"
        [[ -f "$call_file" ]] && seg_file="$call_file"

        [[ -f "$cnr_file" && -f "$seg_file" ]] && {
            
            local scatter_out="${sample_dir}/${sample_name}-scatter.pdf"
            [[ -f "$scatter_out" ]] || {
                cnvkit.py scatter "$cnr_file" -s "$seg_file" -o "$scatter_out" \
                && log_info "Scatter plot created for $sample_name" \
                || log_error "Failed scatter plot for $sample_name"
            }

            local diagram_out="${sample_dir}/${sample_name}-diagram.pdf"
            [[ -f "$diagram_out" ]] || {
                cnvkit.py diagram "$cnr_file" -s "$seg_file" -o "$diagram_out" \
                && log_info "Diagram created for $sample_name" \
                || log_error "Failed diagram for $sample_name"
            }
            
        } || log_error "Missing CNR or CNS for $sample_name. Skipping plots."
    done
}

# =====================================================================
# MAIN EXECUTION
# =====================================================================

log_info "=== STARTING VISUALIZATION SUITE ==="

generate_cohort_heatmaps
generate_single_plots

log_info "=== VISUALIZATION COMPLETED ==="
