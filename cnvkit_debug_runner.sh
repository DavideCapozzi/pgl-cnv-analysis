#!/bin/bash

# =====================================================================
# CNVkit Pipeline: Deconstructed Debug Runner
# Target: PGL WES - Tumor Only - Flat Reference
# Focus: Isolating 'fix' crash
# =====================================================================

# --- CONFIGURATION PATHS ---
INPUT_BASE_DIR="/mnt/d/CNVkit/tumor/PTJ_WES_IDT-30802789"
REF_DIR="/mnt/d/CNVkit/tumor/tumor_targets"
REF_FILE="${REF_DIR}/flat_reference.cnn"
VCF_LIST="${REF_DIR}/vcf_list.txt"
OUT_DIR="/mnt/d/CNVkit/tumor/tumor_test_output"

# Ensure directories exist
mkdir -p "$OUT_DIR"

# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================

log_info() { echo -e "[INFO] $(date '+%H:%M:%S') - $1"; }
log_error() { echo -e "[ERROR] $(date '+%H:%M:%S') - $1" >&2; }

get_vcf_from_list() {
    local sample_name="$1"
    local vcf_list="$2"
    grep "${sample_name}" "$vcf_list" | head -n 1
}

# =====================================================================
# PROCESSING FUNCTION (DECONSTRUCTED)
# =====================================================================

process_sample_manual() {
    local bam_file="$1"
    local sample_name=$(basename "$bam_file" .bam)
    local sample_dir="${OUT_DIR}/${sample_name}"
    
    # Define outputs
    local target_cov="${sample_dir}/${sample_name}.targetcoverage.cnn"
    local anti_cov="${sample_dir}/${sample_name}.antitargetcoverage.cnn"
    local cnr_file="${sample_dir}/${sample_name}.cnr"
    local cns_file="${sample_dir}/${sample_name}.cns"
    local call_file="${sample_dir}/${sample_name}.call.cns"

    # Define Temporary BEDs (relying on artifacts from previous batch run)
    # NOTE: In a clean run, you would provide the original target BEDs.
    local target_bed="${sample_dir}/flat_reference.target-tmp.bed"
    local antitarget_bed="${sample_dir}/flat_reference.antitarget-tmp.bed"

    log_info "Processing sample: ${sample_name}"
    mkdir -p "$sample_dir"

    # --- CHECK REQUIREMENTS ---
    if [[ ! -f "$target_bed" || ! -f "$antitarget_bed" ]]; then
        log_error "Missing temporary BED files in ${sample_dir}. Cannot run manual coverage without targets."
        log_error "Run 'cnvkit.py batch' once to generate beds or point to original BEDs."
        return 1
    fi

    # --- STEP 1: COVERAGE (Explicit) ---
    log_info "Step 1: Running Coverage..."
    
    # Target Coverage
    if [[ ! -f "$target_cov" ]]; then
        cnvkit.py coverage "$bam_file" "$target_bed" -o "$target_cov" -p $(nproc) \
        && log_info "Target coverage generated." \
        || { log_error "Target coverage failed."; return 1; }
    else
        log_info "Target coverage exists. Skipping."
    fi

    # Antitarget Coverage
    if [[ ! -f "$anti_cov" ]]; then
        cnvkit.py coverage "$bam_file" "$antitarget_bed" -o "$anti_cov" -p $(nproc) \
        && log_info "Antitarget coverage generated." \
        || { log_error "Antitarget coverage failed."; return 1; }
    else
        log_info "Antitarget coverage exists. Skipping."
    fi

    # --- STEP 2: FIX (The Critical Fail Point) ---
    log_info "Step 2: Running Fix (Target + Antitarget + Flat Ref)..."
    
    # Attempt standard fix first. 
    # If this fails, we will see the specific python error in stdout/stderr.
    cnvkit.py fix "$target_cov" "$anti_cov" "$REF_FILE" \
        --output "$cnr_file" \
    && log_info "Fix completed successfully. CNR created." \
    || { 
        log_error "Fix FAILED. Retrying with --no-gc --no-edge (Flat Ref Fallback)..."
        
        # Fallback for Flat References lacking GC content
        cnvkit.py fix "$target_cov" "$anti_cov" "$REF_FILE" \
            --output "$cnr_file" \
            --no-gc \
            --no-edge \
        && log_info "Fix completed with GC/Edge correction disabled." \
        || { log_error "Fix CRITICAL FAILURE. Check reference compatibility."; return 1; }
    }

    # --- STEP 3: SEGMENTATION ---
    log_info "Step 3: Running Segmentation (CBS)..."
    
    cnvkit.py segment "$cnr_file" \
        --method cbs \
        --drop-low-coverage \
        --output "$cns_file" \
        --processes $(nproc) \
    && log_info "Segmentation completed." \
    || { log_error "Segmentation failed."; return 1; }

    # --- STEP 4: CALLING ---
    local vcf_path=$(get_vcf_from_list "$sample_name" "$VCF_LIST")
    
    if [[ -n "$vcf_path" && -f "$vcf_path" ]]; then
        log_info "VCF found: $(basename "$vcf_path"). Running Clonal Call..."
        cnvkit.py call "$cns_file" --vcf "$vcf_path" --method clonal --output "$call_file"
    else
        log_info "No VCF. Running basic Call..."
        cnvkit.py call "$cns_file" --method clonal --output "$call_file"
    fi

    echo "---------------------------------------------------"
}

# =====================================================================
# MAIN EXECUTION
# =====================================================================

log_info "=== STARTING DEBUG PIPELINE ==="

# Find BAMs and process
find "$INPUT_BASE_DIR" -type f -name "*-t.bam" | sort | while read bam_file; do
    process_sample_manual "$bam_file"
done

log_info "=== PIPELINE COMPLETED ==="
