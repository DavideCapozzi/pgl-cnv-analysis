#!/bin/bash

# =====================================================================
# CNVkit Pipeline Runner - FIXED for v0.9.10 & Flat Reference
# =====================================================================

# --- CONFIGURAZIONE PATH (Verifica che siano corretti per il tuo sistema) ---
# Directory base dove si trovano i dati
BASE_DIR="/mnt/d/CNVkit/tumor/PTJ_WES_IDT-30802789"
# Directory dove hai salvato la flat reference creata
TARGETS_DIR="/mnt/d/CNVkit/tumor/tumor_targets"
# Directory di output (verrà creata)
OUT_DIR="/mnt/d/CNVkit/tumor/tumor_test_output"

# File di riferimento (Flat Reference)
REF_FILE="${TARGETS_DIR}/flat_reference.cnn"

# Log files
ERROR_LOG="${OUT_DIR}/pipeline_error.log"

# Creazione directory
mkdir -p "$OUT_DIR"

# =====================================================================
# FUNZIONI
# =====================================================================

process_sample() {
    local bam_file="$1"
    local sample_name=$(basename "$bam_file" | sed 's/.bam//')
    local sample_out_dir="${OUT_DIR}/${sample_name}"
    
    echo "[INFO] Processing sample: ${sample_name}"
    mkdir -p "$sample_out_dir"

    # ---------------------------------------------------------
    # STEP 1: CNVkit Batch (Coverage -> Fix -> Segment)
    # ---------------------------------------------------------
    # NOTA: Rimosse opzioni -t, -a, -f perché incompatibili con -r
    # NOTA: Rimosso --prune (deprecated)
    # NOTA: Aggiunto --drop-low-coverage e --segment-method cbs per Flat Ref
    
    echo "[EXEC] Running batch with Flat Reference..."
    
    cnvkit.py batch "$bam_file" \
        -r "$REF_FILE" \
        --output-dir "$sample_out_dir" \
        --segment-method cbs \
        --drop-low-coverage \
        --scatter \
        --diagram \
        -p $(nproc)

    if [ $? -ne 0 ]; then
        echo "[ERROR] Batch failed for ${sample_name}" | tee -a "$ERROR_LOG"
        return 1
    fi

    # ---------------------------------------------------------
    # STEP 2: Integrazione VCF (Call & BAF)
    # ---------------------------------------------------------
    # Cerchiamo il VCF corrispondente
    # Pattern: Cerca file che iniziano con il nome del sample e contengono "hard-filtered"
    vcf_file=$(find "$BASE_DIR" -type f -name "${sample_name}*.hard-filtered.vcf.gz" | head -n 1)

    # Se non trova .gz, cerca .vcf non compresso
    if [ -z "$vcf_file" ]; then
        vcf_file=$(find "$BASE_DIR" -type f -name "${sample_name}*.hard-filtered.vcf" | head -n 1)
    fi

    local cns_file="${sample_out_dir}/${sample_name}.cns"

    if [[ -f "$vcf_file" && -f "$cns_file" ]]; then
        echo "[INFO] Found VCF: $(basename "$vcf_file"). Running CALL with BAF..."
        
        # CNVkit Call usando il VCF per calcolare la BAF
        cnvkit.py call "$cns_file" \
            --vcf "$vcf_file" \
            --method clonal \
            -o "${sample_out_dir}/${sample_name}_call.cns"
            
        echo "[SUCCESS] Call completed with VCF integration."
    else
        echo "[WARNING] VCF not found for ${sample_name} OR .cns missing. Running call without VCF."
        
        cnvkit.py call "$cns_file" \
            --method clonal \
            -o "${sample_out_dir}/${sample_name}_call.cns"
    fi

    # ---------------------------------------------------------
    # STEP 3: Genemetrics
    # ---------------------------------------------------------
    local cnr_file="${sample_out_dir}/${sample_name}.cnr"
    local call_file="${sample_out_dir}/${sample_name}_call.cns"
    
    if [[ -f "$cnr_file" && -f "$call_file" ]]; then
        echo "[EXEC] Generating genemetrics..."
        cnvkit.py genemetrics "$cnr_file" \
            -s "$call_file" \
            --drop-low-coverage \
            -o "${sample_out_dir}/${sample_name}_genemetrics.txt"
    fi
}

# =====================================================================
# MAIN LOOP
# =====================================================================

echo "=== STARTING PIPELINE ==="
echo "Reference: $REF_FILE"

# Verifica esistenza reference
if [ ! -f "$REF_FILE" ]; then
    echo "[FATAL] Reference file not found at $REF_FILE"
    exit 1
fi

# Loop su tutti i file BAM nella directory dei tumori
# Modifica il pattern "-t*.bam" se i tuoi file hanno nomi diversi (es. "*.bam")
find "$BASE_DIR" -type f -name "*-t*.bam" | sort | while read bam_file; do
    process_sample "$bam_file"
done

echo "=== PIPELINE COMPLETED ==="
