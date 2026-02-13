#!/bin/bash

# =====================================================================
# CNVkit Pipeline: Stable Batch Runner (v0.9.10)
# Target: PGL WES - Tumor Only - Flat Reference
# Optimization: Scatter removed to prevent race conditions/missing files
# =====================================================================

# --- CONFIGURATION PATHS ---
INPUT_BASE_DIR="/mnt/d/CNVkit/tumor/PTJ_WES_IDT-30802789"
REF_DIR="/mnt/d/CNVkit/tumor/tumor_targets"
REF_FILE="${REF_DIR}/curated_flat_reference.cnn"
VCF_LIST="${REF_DIR}/vcf_list.txt"


# --- DERIVED PATHS ---
REF_BASENAME=$(basename "${REF_FILE%.*}")
OUT_DIR="/mnt/d/CNVkit/tumor/tumor_newout/${REF_BASENAME}_output"

# Ensure directories exist
mkdir -p "$OUT_DIR"

# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================

log_info() { echo -e "[INFO] $(date '+%H:%M:%S') - $1"; }
log_error() { echo -e "[ERROR] $(date '+%H:%M:%S') - $1" >&2; }

# Extract VCF path from the provided list file using grep
get_vcf_from_list() {
    local sample_name="$1"
    local vcf_list="$2"
    # Search for the sample name in the list and return the first match
    grep "${sample_name}" "$vcf_list" | head -n 1
}

# =====================================================================
# CORE PROCESSING FUNCTION
# =====================================================================

process_sample() {
    local bam_file="$1"
    local sample_name=$(basename "$bam_file" .bam)
    local sample_dir="${OUT_DIR}/${sample_name}"
    
    # Define outputs
    local cnr_file="${sample_dir}/${sample_name}.cnr"
    local cns_file="${sample_dir}/${sample_name}.cns"
    local call_file="${sample_dir}/${sample_name}.call.cns"
    local breaks_file="${sample_dir}/${sample_name}.breaks.txt"
    local genemetrics_file="${sample_dir}/${sample_name}.genemetrics.txt"

    log_info "Processing sample: ${sample_name}"
    mkdir -p "$sample_dir"

    # --- STEP 1: BATCH (Coverage -> Fix -> Segment) ---
    # CRITICAL CHANGE: Removed --scatter and --diagram to prevent file writing errors.
    # Added explicit check [[ -f $cns_file ]] to ensure success.
    
    log_info "Running batch (Coverage + Segmentation)..."
    
    cnvkit.py batch "$bam_file" \
        --reference "$REF_FILE" \
        --output-dir "$sample_dir" \
        --segment-method cbs \
        --drop-low-coverage \
        --processes $(nproc) \
    && [[ -f "$cns_file" ]] \
    && log_info "Batch completed. CNS generated: $cns_file" \
    || { log_error "Batch failed or CNS not created for ${sample_name}"; return 1; }

    # --- STEP 2: CALLING (With VCF Integration) ---
    local vcf_path=$(get_vcf_from_list "$sample_name" "$VCF_LIST")
    
    if [[ -n "$vcf_path" && -f "$vcf_path" ]]; then
        log_info "VCF found: $(basename "$vcf_path"). Running Clonal Call..."
        
        cnvkit.py call "$cns_file" \
            --vcf "$vcf_path" \
            --method clonal \
            --output "$call_file" \
        && log_info "Call (with VCF) completed." \
        || log_error "Call failed despite VCF presence."
    else
        log_info "VCF NOT found in list. Running basic Call..."
        
        cnvkit.py call "$cns_file" \
            --method clonal \
            --output "$call_file" \
        && log_info "Call (basic) completed."
    fi

    # --- STEP 3: DOWNSTREAM (Breaks & Genemetrics) ---
    # Only execute if call file was created successfully
    [[ -f "$call_file" ]] && {
        log_info "Generating metrics..."
        
        # Breaks
        cnvkit.py breaks "$cnr_file" "$cns_file" \
            --min-probes 5 \
            | tee "$breaks_file" > /dev/null
            
        # Genemetrics
        cnvkit.py genemetrics "$cnr_file" \
            -s "$call_file" \
            --drop-low-coverage \
            --output "$genemetrics_file" \
        && log_info "Metrics generated."
    } || log_error "Skipping metrics: Call file missing."

    echo "---------------------------------------------------"
}

# =====================================================================
# HEATMAP MODULE
# =====================================================================

generate_heatmaps() {
    log_info "Generating Cohort Heatmaps..."
    
    local heatmap_dir="${OUT_DIR}/heatmaps"
    mkdir -p "$heatmap_dir"
    
    # Priority: Use .call.cns (ploidy corrected) if available, else .cns
    local cns_list=$(find "$OUT_DIR" -name "*.call.cns")
    
    if [[ -z "$cns_list" ]]; then
        log_error "No .call.cns files found. Trying standard .cns..."
        cns_list=$(find "$OUT_DIR" -name "*.cns" ! -name "*.call.cns" ! -name "*.bintest.cns")
    fi

    [[ -n "$cns_list" ]] && {
        # Global Heatmap
        cnvkit.py heatmap $cns_list -d -o "${heatmap_dir}/cohort_heatmap_all.pdf" && \
        log_info "Global heatmap created."

        # PGL Specific Chromosomes
        for chr in chr1 chr11 chr22; do
             cnvkit.py heatmap $cns_list -d -c "$chr" -o "${heatmap_dir}/cohort_heatmap_${chr}.pdf"
        done
    } || log_error "No CNS files available for heatmap."
}

# =====================================================================
# MAIN EXECUTION
# =====================================================================

log_info "=== STARTING PIPELINE ==="

# Check requirements
[[ ! -f "$REF_FILE" ]] && { log_error "Reference file missing: $REF_FILE"; exit 1; }
[[ ! -f "$VCF_LIST" ]] && { log_error "VCF List missing: $VCF_LIST"; exit 1; }

# Find BAMs and process
find "$INPUT_BASE_DIR" -type f -name "*-t.bam" | sort | while read bam_file; do
    process_sample "$bam_file"
done

# Generate final visualizations
generate_heatmaps

log_info "=== PIPELINE COMPLETED ==="