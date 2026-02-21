#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNVkit Dynamic Curated Flat Reference Builder (v2.1)
------------------------------------------------------
Strategy: "Autosomal Global QC -> Dynamic Noise Evaluation -> Adaptive Fallback"
Context:  Paraganglioma (PGL) Tumor-Only WES analysis.

Changes in v2.1:
  1. I/O Resilience: Auto-handles directory vs file paths for report output.
  2. Pairing Logs: Explicit warnings for missing target/antitarget pairs.
  3. Autosomal QC: Global noise evaluation is strictly restricted to chr1-22
     to prevent physiological sex differences (chrX/Y) from artificially 
     inflating the IQR and falsely rejecting male samples.
  4. Optimized Masking: Pre-computed chromosome masks for CPU efficiency.
"""

import os
import sys
import argparse
import logging
import warnings
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
import cnvlib

# =============================================================================
# LOGGING SETUP
# =============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger("DynamicRefBuilder")

# =============================================================================
# ALGORITHMIC THRESHOLDS
# =============================================================================

# 1. GLOBAL SAMPLE QC (Autosomes Only)
MAX_SAMPLE_NOISE = 0.80 

# 2. CHROMOSOME SHIFT & NOISE
MAX_CHROM_SHIFT = 0.25 
MAX_CHROM_NOISE = 0.65 
MIN_VALID_BINS_FRAC = 0.10

# =============================================================================
# UTILITIES & FILE I/O
# =============================================================================

def normalize_chrom(name: str) -> str:
    s = str(name)
    return s if s.startswith("chr") else f"chr{s}"

def is_autosome(chrom: str) -> bool:
    """Checks if a normalized chromosome string is an autosome (chr1-chr22)."""
    valid_autosomes = {f"chr{i}" for i in range(1, 23)}
    return chrom in valid_autosomes

def locate_files(base_dir: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    targets, antitargets = {}, {}
    base_dir = os.path.abspath(base_dir)
    
    logger.info(f"Scanning directory tree: {base_dir}")
    
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".targetcoverage.cnn"):
                sample_id = file.replace(".targetcoverage.cnn", "")
                targets[sample_id] = os.path.join(root, file)
            elif file.endswith(".antitargetcoverage.cnn"):
                sample_id = file.replace(".antitargetcoverage.cnn", "")
                antitargets[sample_id] = os.path.join(root, file)
                
    common_samples = set(targets.keys()).intersection(set(antitargets.keys()))
    
    # Logging the pairing process explicitly
    for sid in common_samples:
        logger.info(f"  [+] Matched Pair Found: {sid}")
        
    missing_targets = set(antitargets.keys()) - common_samples
    for sid in missing_targets:
        logger.warning(f"  [-] Orphan Antitarget (Missing Target): {sid}")
        
    missing_antis = set(targets.keys()) - common_samples
    for sid in missing_antis:
        logger.warning(f"  [-] Orphan Target (Missing Antitarget): {sid}")
        
    return {s: targets[s] for s in common_samples}, {s: antitargets[s] for s in common_samples}

# =============================================================================
# DYNAMIC EVALUATION ENGINE
# =============================================================================

def evaluate_regions(file_dict: Dict[str, str]) -> Tuple[Dict[str, List[str]], pd.DataFrame, List[str]]:
    sample_ids = sorted(list(file_dict.keys()))
    
    # --- STEP 1: LOAD ALL DATA ---
    template_cnv = cnvlib.read(file_dict[sample_ids[0]])
    chromosomes = np.array([normalize_chrom(c) for c in template_cnv.data['chromosome']])
    unique_chroms = np.unique(chromosomes)
    
    # Pre-compute masks for CPU efficiency
    chrom_masks = {chrom: (chromosomes == chrom) for chrom in unique_chroms}
    autosome_mask = np.array([is_autosome(c) for c in chromosomes])
    
    mat_log2 = np.full((len(template_cnv.data), len(sample_ids)), np.nan, dtype=np.float32)
    for i, sid in enumerate(sample_ids):
        mat_log2[:, i] = cnvlib.read(file_dict[sid]).data['log2'].values.astype(np.float32)

    report_data = []
    
    # --- STEP 2: AUTOSOMAL GLOBAL SAMPLE QC ---
    logger.info("Performing Autosomal Global Sample Quality Control...")
    good_sample_indices = []
    good_sample_ids = []
    
    for i, sid in enumerate(sample_ids):
        # Apply mask to evaluate ONLY chr1-chr22 to avoid sex-bias
        sample_data = mat_log2[autosome_mask, i]
        valid_data = sample_data[~np.isnan(sample_data)]
        
        if len(valid_data) == 0:
            report_data.append((sid, "ALL", np.nan, np.nan, "REJECTED_SAMPLE", "Empty Autosomes"))
            continue
            
        global_iqr = float(np.percentile(valid_data, 75) - np.percentile(valid_data, 25))
        
        if global_iqr > MAX_SAMPLE_NOISE:
            logger.warning(f"  [!] Sample {sid} REJECTED globally (Autosomal IQR: {global_iqr:.2f} > {MAX_SAMPLE_NOISE})")
            report_data.append((sid, "ALL", np.nan, global_iqr, "REJECTED_SAMPLE", "Global Noise Too High"))
        else:
            good_sample_indices.append(i)
            good_sample_ids.append(sid)
            
    if not good_sample_indices:
        raise RuntimeError("All samples failed Global QC. Cannot build reference.")

    # --- STEP 3: POOLED BASELINE ---
    logger.info(f"Building baseline from {len(good_sample_indices)} clean samples...")
    clean_mat = mat_log2[:, good_sample_indices]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        pooled_median = np.nanmedian(clean_mat, axis=1)
        
    residuals = clean_mat - pooled_median[:, np.newaxis]
    
    # --- STEP 4: CHROMOSOME-LEVEL EVALUATION ---
    logger.info("Evaluating chromosomal shifts and local noise...")
    inclusion_map = {sid: [] for sid in good_sample_ids}
    
    for idx, sid in enumerate(good_sample_ids):
        for chrom in unique_chroms:
            mask = chrom_masks[chrom]
            chrom_residuals = residuals[mask, idx]
            valid_residuals = chrom_residuals[~np.isnan(chrom_residuals)]
            
            n_valid, n_total = len(valid_residuals), np.sum(mask)
            
            if n_total == 0 or (n_valid / n_total) < MIN_VALID_BINS_FRAC:
                report_data.append((sid, chrom, np.nan, np.nan, "EXCLUDED", "No Valid Data"))
                continue
                
            chrom_shift = float(np.median(valid_residuals))
            chrom_iqr = float(np.percentile(valid_residuals, 75) - np.percentile(valid_residuals, 25))
            
            if abs(chrom_shift) > MAX_CHROM_SHIFT:
                status, reason = "EXCLUDED", f"Shift out of bounds ({chrom_shift:+.2f})"
            elif chrom_iqr > MAX_CHROM_NOISE:
                status, reason = "EXCLUDED", f"Local noise too high (IQR {chrom_iqr:.2f})"
            else:
                status, reason = "INCLUDED", "Clean"
                inclusion_map[sid].append(chrom)
                
            report_data.append((sid, chrom, chrom_shift, chrom_iqr, status, reason))

    return inclusion_map, pd.DataFrame(report_data, columns=["Sample", "Chromosome", "Median_Shift", "IQR_Noise", "Status", "Reason"]), good_sample_ids

# =============================================================================
# MATRIX BUILDER
# =============================================================================

def build_reference_matrix(file_dict: Dict[str, str], inclusion_map: Dict[str, List[str]], good_sample_ids: List[str]) -> object:
    template_cnv = cnvlib.read(file_dict[good_sample_ids[0]])
    chromosomes = np.array([normalize_chrom(c) for c in template_cnv.data['chromosome']])
    
    mat_log2 = np.full((len(template_cnv.data), len(good_sample_ids)), np.nan, dtype=np.float32)
    mat_depth = np.full((len(template_cnv.data), len(good_sample_ids)), np.nan, dtype=np.float32)
    
    for i, sample_id in enumerate(good_sample_ids):
        curr_cnv = cnvlib.read(file_dict[sample_id])
        allowed_chroms = set(inclusion_map.get(sample_id, []))
        
        if not allowed_chroms:
            continue

        keep_mask = np.isin(chromosomes, list(allowed_chroms))
        mat_log2[keep_mask, i] = curr_cnv.data['log2'].values.astype(np.float32)[keep_mask]
        mat_depth[keep_mask, i] = curr_cnv.data['depth'].values.astype(np.float32)[keep_mask]

    logger.info("Calculating robust statistics and enforcing safety clamps...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ref_log2 = np.nanmean(mat_log2, axis=1)
        ref_depth = np.nanmean(mat_depth, axis=1)
        ref_spread = np.nanstd(mat_log2, axis=1, ddof=1) 

    fallback_mask = np.isnan(ref_log2)
    
    # Depth Clamping
    ref_depth = np.maximum(ref_depth, 1e-6)
    ref_depth[np.isnan(ref_depth)] = np.nanmean(ref_depth) if not np.isnan(np.nanmean(ref_depth)) else 1.0

    # Spread Clamping (75th percentile for conservative weighting)
    safe_spread_val = np.nanpercentile(ref_spread, 75) if not np.isnan(np.nanpercentile(ref_spread, 75)) else 0.1
    ref_spread[np.isnan(ref_spread)] = safe_spread_val
    ref_spread = np.maximum(ref_spread, 0.001)

    ref_log2[fallback_mask] = 0.0

    final_cnv = template_cnv.copy()
    final_cnv.data['log2'] = ref_log2
    final_cnv.data['depth'] = ref_depth
    final_cnv.data['spread'] = ref_spread
    final_cnv.data['weight'] = 1.0 / (ref_spread ** 2)
    
    return final_cnv

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="CNVkit Dynamic Curated Reference Builder")
    parser.add_argument("-i", "--input", required=True, help="Path to directory containing .cnn files")
    parser.add_argument("-o", "--output", required=True, help="Output path for the curated reference.cnn")
    parser.add_argument("-r", "--report", required=False, help="Optional output path for TSV report (file or directory)")
    args = parser.parse_args()
    
    logger.info("--- Phase 1: Locating Files ---")
    targets_map, antitargets_map = locate_files(args.input)
    
    logger.info("--- Phase 2: Dynamic Evaluation ---")
    inclusion_map, report_df, good_ids = evaluate_regions(targets_map)
    
    if args.report:
        report_path = args.report
        # Auto-handle if user passes a directory instead of a file
        if os.path.isdir(report_path):
            report_path = os.path.join(report_path, "evaluation_report.tsv")
        
        report_dir = os.path.dirname(os.path.abspath(report_path))
        if report_dir and not os.path.exists(report_dir):
            os.makedirs(report_dir)
            
        report_df.to_csv(report_path, sep='\t', index=False)
        logger.info(f"Evaluation report saved to: {report_path}")

    logger.info("--- Phase 3: Building Reference Matrices ---")
    ref_targets = build_reference_matrix(targets_map, inclusion_map, good_ids)
    ref_antitargets = build_reference_matrix(antitargets_map, inclusion_map, good_ids)

    logger.info("--- Phase 4: Merging and Saving ---")
    ref_targets.add(ref_antitargets)
    ref_targets.sort()
    
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    logger.info(f"Writing output to: {args.output}")
    ref_targets.data.to_csv(args.output, sep='\t', index=False, float_format='%.6g')
    logger.info("Pipeline completed successfully.")

if __name__ == "__main__":
    main()