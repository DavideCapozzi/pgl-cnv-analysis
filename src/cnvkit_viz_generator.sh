#!/bin/bash

# =====================================================================
# CNVkit Module: Visualization Suite
# Description: Generates Cohort Heatmaps and Single Sample plots.
# Input: Scans OUT_DIR for existing .cns and .cnr files.
# =====================================================================

# --- CONFIGURATION PATHS ---
# Must match the output dir of the processor script
OUT_DIR="/mnt/d/CNVkit/tumor/tumor_test_output"
VIZ_DIR="${OUT_DIR}/visualizations"

mkdir -p "$VIZ_DIR"

# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================

log_info() { echo -e "[INFO] $(date '+%H:%M:%S') - $1"; }
log_error() { echo -e "[ERROR] $(date '+%H:%M:%S') - $1" >&2; }

# =====================================================================
# COHORT HEATMAPS
# =====================================================================

generate_cohort_heatmaps() {
    log_info "Collecting files for Cohort Heatmaps..."
    
    # Find all .call.cns files (preferred) or .cns files
    # We use 'find' to get full paths
    local cns_list=$(find "$OUT_DIR" -name "*.call.cns" | sort)
    
    # Fallback if no .call.cns found
    [[ -z "$cns_list" ]] && cns_list=$(find "$OUT_DIR" -name "*.cns" ! -name "*.call.cns" ! -name "*.bintest.cns" | sort)

    [[ -z "$cns_list" ]] && {
        log_error "No .cns or .call.cns files found in $OUT_DIR. Cannot generate heatmaps."
        return 1
    }

    # Count samples
    local count=$(echo "$cns_list" | wc -l)
    log_info "Found $count samples. Generating Global Heatmap..."

    # 1. Global Heatmap (All Chromosomes)
    cnvkit.py heatmap $cns_list \
        -d \
        --output "${VIZ_DIR}/cohort_heatmap_global.pdf" \
    && log_info "Global heatmap created." \
    || log_error "Failed to create global heatmap."

    # 2. PGL Specific Chromosomes (1, 11, 22)
    # Allows detailed view of typical PGL losses
    for chr in chr1 chr11 chr22; do
         cnvkit.py heatmap $cns_list \
            -d \
            -c "$chr" \
            --output "${VIZ_DIR}/cohort_heatmap_${chr}.pdf" \
         && log_info "Heatmap for $chr created."
    done
}

# =====================================================================
# SINGLE SAMPLE PLOTS
# =====================================================================

generate_single_plots() {
    log_info "Generating Single Sample Plots (Scatter & Diagram)..."

    # Iterate over subdirectories in OUT_DIR
    find "$OUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort | while read sample_dir; do
        
        local sample_name=$(basename "$sample_dir")
        local cnr_file="${sample_dir}/${sample_name}.cnr"
        local cns_file="${sample_dir}/${sample_name}.cns"
        local call_file="${sample_dir}/${sample_name}.call.cns" # Preferred for diagrams
        
        # Use call.cns if available, otherwise cns
        local seg_file="$cns_file"
        [[ -f "$call_file" ]] && seg_file="$call_file"

        # Check if required files exist
        [[ -f "$cnr_file" && -f "$seg_file" ]] && {
            
            # A. Scatter Plot (Whole Genome)
            # Essential for quality control and visual verification
            local scatter_out="${sample_dir}/${sample_name}-scatter.pdf"
            [[ ! -f "$scatter_out" ]] && {
                cnvkit.py scatter "$cnr_file" -s "$seg_file" -o "$scatter_out" \
                && log_info "Scatter plot created for $sample_name"
            }

            # B. Diagram (Ideogram)
            # Good for summary and clinical reporting
            local diagram_out="${sample_dir}/${sample_name}-diagram.pdf"
            [[ ! -f "$diagram_out" ]] && {
                cnvkit.py diagram "$cnr_file" -s "$seg_file" -o "$diagram_out" \
                && log_info "Diagram created for $sample_name"
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
