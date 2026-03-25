#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNVkit Minimal Common Region (MCR) & Variant Analyzer
------------------------------------------------------
Purpose: Identifies overlapping genomic regions (MCRs) of concordant copy number 
alterations across multiple patients, and subsequently extracts mutations 
falling within these intervals from associated VCF files.

Architecture:
  1. Recursive directory traversal with strict suffix-stripping matching.
  2. Segment contiguous fusion (prevents technical fragmentation).
  3. O(N log N) Sweep-line algorithm for precise overlap intersections.
  4. Tabix-indexed Pysam queries for O(1) variant extraction.
"""

import os
import sys
import argparse
import logging
from typing import List, Dict, Tuple, Set

import pandas as pd
import pysam

# =============================================================================
# LOGGING SETUP
# =============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger("MCR_Analyzer")

# =============================================================================
# UTILITIES & I/O
# =============================================================================

def normalize_chrom(name: str) -> str:
    """Ensures consistent chromosome naming (chr prefix)."""
    s = str(name).strip()
    return s if s.startswith("chr") else f"chr{s}"

def _extract_core_id(filename: str, file_type: str) -> str:
    """
    Safely extracts the core sample ID by removing known exact pipeline suffixes.
    This prevents matching failures caused by complex DRAGEN/CNVkit suffixes.
    """
    if file_type == "cns":
        suffixes = [".call.cns", ".cns"]
    else:
        # Prioritize hard-filtered as requested
        suffixes = [".hard-filtered.vcf.gz", ".hard-filtered.vcf", ".vcf.gz", ".vcf"]
        
    core_id = filename
    for suff in suffixes:
        if core_id.endswith(suff):
            core_id = core_id[:-len(suff)]
            break
            
    return core_id

def map_sample_files(cnv_dir: str, vcf_dir: str) -> Dict[str, Dict[str, str]]:
    """
    Recursively pairs .call.cns files with corresponding VCF files.
    Matches based on the core sample ID.
    """
    cnv_dir = os.path.abspath(cnv_dir)
    vcf_dir = os.path.abspath(vcf_dir)
    
    sample_map = {}
    
    # 1. Recursive locate CNV files
    for root, _, files in os.walk(cnv_dir):
        for file in files:
            if file.endswith(".call.cns"):
                core_id = _extract_core_id(file, "cns")
                if core_id not in sample_map:
                    sample_map[core_id] = {'cns': None, 'vcf': None, 'sample_name': core_id}
                sample_map[core_id]['cns'] = os.path.join(root, file)
                
    # 2. Recursive locate VCF files
    for root, _, files in os.walk(vcf_dir):
        for file in files:
            # Check for hard-filtered VCFs to avoid capturing raw intermediate files
            if "hard-filtered.vcf" in file:
                core_id = _extract_core_id(file, "vcf")
                if core_id in sample_map:
                    sample_map[core_id]['vcf'] = os.path.join(root, file)

    # 3. Validate mapping
    valid_pairs = {k: v for k, v in sample_map.items() if v['cns'] is not None}
    missing_vcfs = [k for k, v in valid_pairs.items() if v['vcf'] is None]
    
    logger.info(f"Discovered {len(valid_pairs)} CNV profiles.")
    if missing_vcfs:
        logger.warning(f"Missing VCFs for {len(missing_vcfs)} samples. Variant analysis will be skipped for them.")
        logger.debug(f"Unmatched CNV IDs: {', '.join(missing_vcfs[:5])}...")
        
    return valid_pairs

# =============================================================================
# CORE INTERVAL LOGIC
# =============================================================================

def parse_and_merge_cnv(file_path: str, sample_id: str, chroms: List[str]) -> Tuple[List[Dict], List[Dict]]:
    """
    Parses CNVkit output and merges contiguous segments of the same alteration type.
    We isolate Amplifications (CN > 2) and Deletions (CN < 2).
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        logger.error(f"Failed to read {file_path}: {e}")
        return [], []

    df['chromosome'] = df['chromosome'].apply(normalize_chrom)
    df = df[df['chromosome'].isin(chroms)]
    
    amps = df[df['cn'] > 2].copy()
    dels = df[df['cn'] < 2].copy()
    
    def merge_segments(sub_df: pd.DataFrame, alt_type: str) -> List[Dict]:
        if sub_df.empty:
            return []
            
        sub_df = sub_df.sort_values(by=['chromosome', 'start'])
        merged = []
        
        for _, row in sub_df.iterrows():
            if not merged:
                merged.append({
                    'chrom': row['chromosome'], 'start': row['start'], 
                    'end': row['end'], 'sample': sample_id, 'type': alt_type
                })
            else:
                last = merged[-1]
                # Fuse contiguous or overlapping segments
                if row['chromosome'] == last['chrom'] and row['start'] <= last['end'] + 1:
                    last['end'] = max(last['end'], row['end'])
                else:
                    merged.append({
                        'chrom': row['chromosome'], 'start': row['start'], 
                        'end': row['end'], 'sample': sample_id, 'type': alt_type
                    })
        return merged

    return merge_segments(amps, 'AMP'), merge_segments(dels, 'DEL')

def sweep_line_mcr(segments: List[Dict], min_patients: int) -> List[Dict]:
    """
    O(N log N) Sweep-line algorithm to find Minimal Common Regions (MCR).
    """
    if not segments:
        return []
        
    events = []
    for seg in segments:
        events.append((seg['chrom'], seg['start'], 'start', seg['sample']))
        events.append((seg['chrom'], seg['end'], 'end', seg['sample']))
        
    # Order: chrom -> position -> 'end' processes before 'start' at exact same position
    events.sort(key=lambda x: (x[0], x[1], 0 if x[2] == 'end' else 1))
    
    active_samples: Set[str] = set()
    mcrs = []
    current_mcr_start = None
    
    for chrom, pos, event_type, sample in events:
        previous_count = len(active_samples)
        
        if event_type == 'start':
            active_samples.add(sample)
        else:
            active_samples.discard(sample)
            
        current_count = len(active_samples)
        
        # MCR Start: Crossed threshold upwards
        if previous_count < min_patients and current_count >= min_patients:
            current_mcr_start = pos
            
        # MCR End: Crossed threshold downwards
        elif previous_count >= min_patients and current_count < min_patients:
            if current_mcr_start is not None and pos > current_mcr_start:
                mcrs.append({
                    'chrom': chrom,
                    'start': current_mcr_start,
                    'end': pos,
                    'n_patients': previous_count,
                    'patients': ",".join(sorted(list(active_samples.union({sample}))))
                })
            current_mcr_start = None
            
    return mcrs

# =============================================================================
# VARIANT EXTRACTION
# =============================================================================

def extract_variants_in_mcr(mcr: Dict, sample_map: Dict[str, Dict[str, str]], mcr_id: str) -> List[Dict]:
    """
    Uses pysam to extract variants for the patients involved in the specific MCR.
    """
    involved_samples = mcr['patients'].split(",")
    variants_found = []
    
    for sample in involved_samples:
        vcf_path = sample_map.get(sample, {}).get('vcf')
        if not vcf_path or not os.path.exists(vcf_path):
            continue
            
        try:
            vcf = pysam.VariantFile(vcf_path)
        except Exception as e:
            logger.error(f"Failed to open VCF for {sample}: {e}")
            continue
            
        # Fetch using TBI index. Fallback to linear iteration if index is missing.
        try:
            records = vcf.fetch(mcr['chrom'], mcr['start'], mcr['end'])
            for rec in records:
                variants_found.append(_parse_vcf_record(rec, sample, mcr_id, mcr['chrom']))
        except ValueError:
            logger.warning(f"VCF {vcf_path} lacks TBI index. Parsing linearly within contig...")
            vcf.reset()
            for rec in vcf:
                if rec.chrom == mcr['chrom'] and mcr['start'] <= rec.pos <= mcr['end']:
                    variants_found.append(_parse_vcf_record(rec, sample, mcr_id, mcr['chrom']))
        
        vcf.close()
        
    return variants_found

def _parse_vcf_record(rec: pysam.VariantRecord, sample: str, mcr_id: str, chrom: str) -> Dict:
    """Extracts standard genetic metrics from a VCF record."""
    alt = ",".join(rec.alts) if rec.alts else "."
    
    dp = "."
    af = "."
    if rec.samples:
        sample_data = rec.samples[0] 
        dp = sample_data.get('DP', ".")
        af = sample_data.get('AF', ".")
        if isinstance(af, tuple):
            af = ",".join(map(str, af))
            
    return {
        'MCR_ID': mcr_id,
        'Sample': sample,
        'Chromosome': chrom,
        'Position': rec.pos,
        'Ref': rec.ref,
        'Alt': alt,
        'Filter': ";".join(rec.filter.keys()) if rec.filter.keys() else "PASS",
        'Depth': dp,
        'Allele_Freq': af
    }

# =============================================================================
# MAIN ORCHESTRATOR
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Extracts Minimal Common Regions (MCRs) from CNVkit output and maps associated VCF variants.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example usage:
  python src/mcr_variant_analyzer.py \\
    --cnv-dir /mnt/d/CNVkit/tumor/tumor_newout/dynamic_flat_reference_output/ \\
    --vcf-dir /mnt/d/CNVkit/tumor/PTJ_WES_IDT-30802789/ \\
    --out-dir results/ \\
    --chromosomes 1,3,7,11,17,22 \\
    --min-patients 3
        """
    )
    
    parser.add_argument("--cnv-dir", required=True, 
                        help="Recursive root directory containing *.call.cns files.")
    parser.add_argument("--vcf-dir", required=True, 
                        help="Recursive root directory containing *.hard-filtered.vcf.gz files.")
    parser.add_argument("--out-dir", required=True, 
                        help="Output directory for TSV reports.")
    parser.add_argument("--chromosomes", required=True, 
                        help="Comma-separated list of chromosomes to analyze.\nValid formats: '1,3,7' or 'chr1,chr3,chr7'.")
    parser.add_argument("--min-patients", type=int, default=3, 
                        help="Minimum number of samples that must share a CNV alteration to form an MCR (default: 3).")
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)

    # Dynamically add file handler to the logger to save logs in the results folder
    log_file_path = os.path.join(args.out_dir, "mcr_analysis_run.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
    logger.addHandler(file_handler)
    logger.info(f"Initialized file logging at: {log_file_path}")

    target_chroms = [normalize_chrom(c.strip()) for c in args.chromosomes.split(",")]
    
    logger.info("--- Phase 1: File Mapping ---")
    sample_map = map_sample_files(args.cnv_dir, args.vcf_dir)
    if not sample_map:
        logger.error("No valid CNV files found. Exiting.")
        return

    logger.info(f"--- Phase 2: Interval Parsing & MCR Detection (Target: {args.min_patients}+ patients) ---")
    all_amps, all_dels = [], []
    
    for sample, paths in sample_map.items():
        amps, dels = parse_and_merge_cnv(paths['cns'], sample, target_chroms)
        all_amps.extend(amps)
        all_dels.extend(dels)
        
    mcr_amps = sweep_line_mcr(all_amps, args.min_patients)
    mcr_dels = sweep_line_mcr(all_dels, args.min_patients)
    
    for m in mcr_amps: m['type'] = 'AMP'
    for m in mcr_dels: m['type'] = 'DEL'
    
    all_mcrs = mcr_amps + mcr_dels
    logger.info(f"Detected {len(all_mcrs)} Minimal Common Regions ({len(mcr_amps)} AMP, {len(mcr_dels)} DEL).")
    
    if not all_mcrs:
        logger.info("No MCRs met the threshold criteria. Pipeline successfully terminated.")
        return

    logger.info("--- Phase 3: VCF Integration ---")
    mcr_records = []
    variant_records = []
    
    for i, mcr in enumerate(all_mcrs, 1):
        mcr_id = f"MCR_{mcr['chrom']}_{mcr['type']}_{i:03d}"
        
        mcr_records.append({
            'MCR_ID': mcr_id,
            'Chromosome': mcr['chrom'],
            'Start': mcr['start'],
            'End': mcr['end'],
            'Type': mcr['type'],
            'N_Patients': mcr['n_patients'],
            'Patient_IDs': mcr['patients']
        })
        
        variants = extract_variants_in_mcr(mcr, sample_map, mcr_id)
        variant_records.extend(variants)

    logger.info("--- Phase 4: Report Generation ---")
    
    df_mcr = pd.DataFrame(mcr_records)
    out_mcr_path = os.path.join(args.out_dir, "mcr_regions_summary.xlsx")
    df_mcr.to_excel(out_mcr_path, index=False, engine='openpyxl')
    logger.info(f"[+] Saved MCR regions to: {out_mcr_path}")
    
    if variant_records:
        df_var = pd.DataFrame(variant_records)
        out_var_path = os.path.join(args.out_dir, "mcr_patient_variants.xlsx")
        try:
            df_var.to_excel(out_var_path, index=False, engine='openpyxl')
            logger.info(f"[+] Saved {len(variant_records)} variant annotations to: {out_var_path}")
        except Exception as e:
            logger.error(f"Failed to save variants to Excel (possibly exceeded 1M rows). Error: {e}")
    else:
        logger.info("[-] No variants found within the detected MCRs for the involved patients.")
        
    logger.info("Analysis complete.")

if __name__ == "__main__":
    main()