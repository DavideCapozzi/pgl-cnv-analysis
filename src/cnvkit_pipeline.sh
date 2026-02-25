#!/bin/bash

# =====================================================================
# CNVkit Pipeline: Stable Batch Runner (v0.9.11)
# Target: PGL WES - Tumor Only - Flat Reference
# Optimization: Smart skipping, VCF baf check, timestamp validation
# =====================================================================

# --- CONFIGURATION PATHS ---
INPUT_BASE_DIR="/mnt/d/CNVkit/tumor/PTJ_WES_IDT-30802789"
REF_DIR="/mnt/d/CNVkit/tumor/tumor_targets"
REF_FILE="${REF_DIR}/dynamic_flat_reference.cnn"
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
    grep "${sample_name}" "$vcf_list" | head -n 1
}

# =====================================================================
# CORE PROCESSING FUNCTION
# =====================================================================

process_sample() {
    local bam_file="$1"
    local current_idx="$2"
    local total_bams="$3"
    
    local sample_name=$(basename "$bam_file" .bam)
    local sample_dir="${OUT_DIR}/${sample_name}"
    
    # Define outputs
    local cnr_file="${sample_dir}/${sample_name}.cnr"
    local cns_file="${sample_dir}/${sample_name}.cns"
    local call_file="${sample_dir}/${sample_name}.call.cns"
    local breaks_file="${sample_dir}/${sample_name}.breaks.txt"
    local genemetrics_file="${sample_dir}/${sample_name}.genemetrics.txt"

    log_info "--- Processing sample [${current_idx}/${total_bams}]: ${sample_name} ---"
    mkdir -p "$sample_dir"

    # --- STEP 1: BATCH (Coverage -> Fix -> Segment -> Basic Call) ---
    [[ -f "$cns_file" && -f "$cnr_file" ]] && log_info "Skip Batch: $cns_file and $cnr_file already exist." || {
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
    }

    # --- STEP 2: CALLING (With VCF Integration & BAF Check) ---
    local vcf_path=$(get_vcf_from_list "$sample_name" "$VCF_LIST")
    
    [[ -n "$vcf_path" && -f "$vcf_path" ]] && {
        # VCF exists: Skip ONLY if call_file exists AND has the 'baf' column
        [[ -f "$call_file" ]] && head -n 1 "$call_file" | grep -qw "baf" && log_info "Skip Calling: $call_file (with VCF BAF) already exists." || {
            log_info "VCF found: $(basename "$vcf_path"). Running Clonal Call to integrate VCF..."
            cnvkit.py call "$cns_file" --vcf "$vcf_path" --method clonal --output "$call_file" \
            && log_info "Call (with VCF) completed." \
            || log_error "Call failed despite VCF presence."
        }
    } || {
        # No VCF found: Skip if any call_file exists, otherwise run basic call
        [[ -f "$call_file" ]] && log_info "Skip Calling: $call_file (basic) already exists." || {
            log_info "VCF NOT found in list. Running basic Call..."
            cnvkit.py call "$cns_file" --method clonal --output "$call_file" \
            && log_info "Call (basic) completed." \
            || log_error "Basic Call failed."
        }
    }

    # --- STEP 3: DOWNSTREAM (Breaks & Genemetrics) ---
    [[ ! -f "$call_file" ]] && log_error "Skipping metrics: Call file $call_file missing." || {
        # Skip ONLY if metrics exist AND are newer than the call_file
        [[ -f "$breaks_file" && -f "$genemetrics_file" && "$breaks_file" -nt "$call_file" && "$genemetrics_file" -nt "$call_file" ]] && log_info "Skip Metrics: updated $breaks_file and $genemetrics_file already exist." || {
            log_info "Generating downstream metrics..."
            
            # Breaks
            cnvkit.py breaks "$cnr_file" "$cns_file" --min-probes 5 | tee "$breaks_file" > /dev/null
                
            # Genemetrics
            cnvkit.py genemetrics "$cnr_file" -s "$call_file" --drop-low-coverage --output "$genemetrics_file" \
            && log_info "Metrics successfully generated."
        }
    }

    echo "---------------------------------------------------"
}

# =====================================================================
# MAIN EXECUTION
# =====================================================================

log_info "=== STARTING PIPELINE ==="

# Check requirements
[[ ! -f "$REF_FILE" ]] && { log_error "Reference file missing: $REF_FILE"; exit 1; }
[[ ! -f "$VCF_LIST" ]] && { log_error "VCF List missing: $VCF_LIST"; exit 1; }

# Count total target BAMs for the counter
TOTAL_BAMS=$(find "$INPUT_BASE_DIR" -type f -name "*-t.bam" | wc -l)
[[ "$TOTAL_BAMS" -eq 0 ]] && { log_error "No BAM files found in $INPUT_BASE_DIR"; exit 1; }

log_info "Found $TOTAL_BAMS samples to process."

# Find BAMs, sort them, and process with a counter
CURRENT_IDX=1
find "$INPUT_BASE_DIR" -type f -name "*-t.bam" | sort | while read -r bam_file; do
    process_sample "$bam_file" "$CURRENT_IDX" "$TOTAL_BAMS"
    ((CURRENT_IDX++))
done

log_info "=== PIPELINE COMPLETED ==="