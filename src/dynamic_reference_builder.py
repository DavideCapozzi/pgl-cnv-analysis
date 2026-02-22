#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNVkit Dynamic Curated Flat Reference Builder (v3.1)
------------------------------------------------------
Strategy: "Autosomal Global QC -> Biological Segment Masking -> Robust Consensus"
Context:  Paraganglioma (PGL) Tumor-Only WES analysis.

Core Algorithm (v3.1):
  Overcomes the "Cohort Flattening" paradox. Computes a 75th-percentile crude baseline
  to anchor to the diploid state, applies a modular rolling median to isolate broad 
  somatic CNVs from bin-level technical noise, and actively masks pathological 
  segments (NaN) before computing the final technical baseline and DDOF=1 spread.
"""

import os
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

# 2. CHROMOSOME EVALUATION
MAX_CHROM_SHIFT = 0.25 
MAX_CHROM_NOISE = 0.65 
MIN_VALID_BINS_FRAC = 0.10

# 3. BIOLOGICAL SEGMENT MASKING
BIOLOGICAL_SMOOTHING_WINDOW = 50  # Bins (captures broad segmental/arm-level events)
MAX_BIOLOGICAL_DEVIATION = 0.30   # log2 threshold to mask out tumor segments

# =============================================================================
# UTILITIES & FILE I/O
# =============================================================================

def normalize_chrom(name: str) -> str:
    s = str(name)
    return s if s.startswith("chr") else f"chr{s}"

def is_autosome(chrom: str) -> bool:
    valid_autosomes = {f"chr{i}" for i in range(1, 23)}
    return chrom in valid_autosomes

def locate_files(base_dir: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    targets, antitargets = {}, {}
    base_dir = os.path.abspath(base_dir)
    
    logger.info(f"Scanning directory tree: {base_dir}")
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".targetcoverage.cnn"):
                sid = file.replace(".targetcoverage.cnn", "")
                targets[sid] = os.path.join(root, file)
            elif file.endswith(".antitargetcoverage.cnn"):
                sid = file.replace(".antitargetcoverage.cnn", "")
                antitargets[sid] = os.path.join(root, file)
                
    common_samples = set(targets.keys()).intersection(set(antitargets.keys()))
    

    for sid in common_samples:
        logger.info(f"  [+] Matched Pair Found: {sid}")
        
    for sid in (set(antitargets.keys()) - common_samples):
        logger.warning(f"  [-] Orphan Antitarget (Missing Target): {sid}")
        
    for sid in (set(targets.keys()) - common_samples):
        logger.warning(f"  [-] Orphan Target (Missing Antitarget): {sid}")
        
    return {s: targets[s] for s in common_samples}, {s: antitargets[s] for s in common_samples}

# =============================================================================
# CORE ALGORITHMS
# =============================================================================

def _apply_biological_mask(mat_log2: np.ndarray, base_track: np.ndarray, 
                           window: int, threshold: float) -> np.ndarray:
    """
    Isolates and masks broad biological events using a rolling median low-pass filter.
    Returns a boolean mask of positions to replace with NaN to protect the baseline.
    """
    tumor_mask = np.zeros_like(mat_log2, dtype=bool)
    
    for i in range(mat_log2.shape[1]):
        sample_track = mat_log2[:, i]
        
        # Skip if the track is mostly empty
        if np.isnan(sample_track).mean() > 0.90:
            continue
            
        residuals = sample_track - base_track
        
        # Pandas rolling median isolates segment-level biological shifts
        smoothed_signal = pd.Series(residuals).rolling(
            window=window, center=True, min_periods=1
        ).median().values
        
        tumor_mask[:, i] = np.abs(smoothed_signal) > threshold
        
    return tumor_mask

# =============================================================================
# DYNAMIC EVALUATION ENGINE
# =============================================================================

def evaluate_regions(file_dict: Dict[str, str]) -> Tuple[Dict[str, List[str]], pd.DataFrame, List[str]]:
    sample_ids = sorted(list(file_dict.keys()))
    template_cnv = cnvlib.read(file_dict[sample_ids[0]])
    chromosomes = np.array([normalize_chrom(c) for c in template_cnv.data['chromosome']])
    unique_chroms = np.unique(chromosomes)
    
    chrom_masks = {chrom: (chromosomes == chrom) for chrom in unique_chroms}
    autosome_mask = np.array([is_autosome(c) for c in chromosomes])
    
    mat_log2 = np.full((len(template_cnv.data), len(sample_ids)), np.nan, dtype=np.float32)
    for i, sid in enumerate(sample_ids):
        mat_log2[:, i] = cnvlib.read(file_dict[sid]).data['log2'].values.astype(np.float32)

    report_data = []
    
    # --- STEP 1: AUTOSOMAL GLOBAL SAMPLE QC ---
    logger.info("Performing Autosomal Global Sample Quality Control...")
    good_sample_indices, good_sample_ids = [], []
    
    for i, sid in enumerate(sample_ids):
        sample_data = mat_log2[autosome_mask, i]
        valid_data = sample_data[~np.isnan(sample_data)]
        
        if len(valid_data) == 0:
            report_data.append((sid, "ALL", np.nan, np.nan, "REJECTED_SAMPLE", "Empty Autosomes"))
            continue
            
        global_iqr = float(np.percentile(valid_data, 75) - np.percentile(valid_data, 25))
        
        if global_iqr > MAX_SAMPLE_NOISE:
            logger.warning(f"  [!] Sample {sid} REJECTED globally (Autosomal IQR: {global_iqr:.2f})")
            report_data.append((sid, "ALL", np.nan, global_iqr, "REJECTED_SAMPLE", "Global Noise Too High"))
        else:
            good_sample_indices.append(i)
            good_sample_ids.append(sid)
            
    if not good_sample_indices:
        raise RuntimeError("All samples failed Global QC. Cannot build reference.")

    clean_mat = mat_log2[:, good_sample_indices]
    
    # --- STEP 2: ROBUST CHROMOSOME-LEVEL EVALUATION ---
    logger.info("Evaluating chromosomal shifts against crude 75th-percentile baseline...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        crude_baseline = np.nanpercentile(clean_mat, 75, axis=1)
        
    residuals = clean_mat - crude_baseline[:, np.newaxis]
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

    report_df = pd.DataFrame(report_data, columns=["Sample", "Chromosome", "Median_Shift", "IQR_Noise", "Status", "Reason"])
    return inclusion_map, report_df, good_sample_ids

# =============================================================================
# MATRIX BUILDER WITH BIOLOGICAL MASKING
# =============================================================================

def build_reference_matrix(file_dict: Dict[str, str], inclusion_map: Dict[str, List[str]], good_sample_ids: List[str]) -> object:
    template_cnv = cnvlib.read(file_dict[good_sample_ids[0]])
    chromosomes = np.array([normalize_chrom(c) for c in template_cnv.data['chromosome']])
    unique_chroms = np.unique(chromosomes)
    chrom_masks = {chrom: (chromosomes == chrom) for chrom in unique_chroms}
    
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

    # --- BIOLOGICAL SEGMENT MASKING ---
    logger.info("Applying Biological Segment Masking to isolate technical noise...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        crude_baseline = np.nanpercentile(mat_log2, 75, axis=1)
        
    masked_segments_count = 0
    sample_mask_counts = {sid: 0 for sid in good_sample_ids}
    
    for chrom in unique_chroms:
        mask = chrom_masks[chrom]
        if not np.any(mask): 
            continue
            
        chrom_log2 = mat_log2[mask, :]
        chrom_base = crude_baseline[mask]
        
        tumor_mask = _apply_biological_mask(
            mat_log2=chrom_log2, 
            base_track=chrom_base, 
            window=BIOLOGICAL_SMOOTHING_WINDOW, 
            threshold=MAX_BIOLOGICAL_DEVIATION
        )
        
        if np.any(tumor_mask):
            mat_log2[mask, :][tumor_mask] = np.nan
            masked_segments_count += np.sum(tumor_mask)
            
            # Track masking per sample for reporting
            for i, sid in enumerate(good_sample_ids):
                sample_mask_counts[sid] += np.sum(tumor_mask[:, i])
                
    logger.info(f"Total bins masked across cohort: {masked_segments_count}")
    for sid, count in sample_mask_counts.items():
        if count > 0:
            logger.debug(f"  [-] {sid}: {count} bins masked")

    # --- FINAL CONSENSUS & SAFETY CLAMPING ---
    logger.info("Calculating robust statistics and enforcing safety clamps...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ref_log2 = np.nanmean(mat_log2, axis=1)
        ref_depth = np.nanmean(mat_depth, axis=1)
        # CRITICAL: ddof=1 prevents infinite weights when n=1. 
        ref_spread = np.nanstd(mat_log2, axis=1, ddof=1) 

    fallback_mask = np.isnan(ref_log2)
    
    # Depth Clamping: Prevents division by zero or negative depths downstream
    ref_depth = np.maximum(ref_depth, 1e-6)
    ref_depth[np.isnan(ref_depth)] = np.nanmean(ref_depth) if not np.isnan(np.nanmean(ref_depth)) else 1.0

    # Spread Clamping: Use 75th percentile for bins lacking variance data (conservative weighting)
    safe_spread_val = np.nanpercentile(ref_spread, 75) if not np.isnan(np.nanpercentile(ref_spread, 75)) else 0.1
    ref_spread[np.isnan(ref_spread)] = safe_spread_val
    # Hard clamp to prevent w = 1/0^2
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
    # Added verbosity flag to see debug info (like per-sample masking)
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose/debug logging")
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info("--- Phase 1: Locating Files ---")
    targets_map, antitargets_map = locate_files(args.input)
    
    logger.info("--- Phase 2: Dynamic Evaluation ---")
    inclusion_map, report_df, good_ids = evaluate_regions(targets_map)
    
    # --- TERMINAL SUMMARY ---
    total_samples = len(targets_map)
    passed_global_qc = len(good_ids)
    failed_global_qc = total_samples - passed_global_qc
    
    logger.info("=== EVALUATION SUMMARY ===")
    logger.info(f"Global Autosomal QC: {passed_global_qc}/{total_samples} samples passed.")
    if failed_global_qc > 0:
        logger.warning(f"Global Autosomal QC: {failed_global_qc} samples entirely rejected.")

    exclusion_df = report_df[report_df["Status"] == "EXCLUDED"]
    if not exclusion_df.empty:
        chrom_exclusions = exclusion_df[exclusion_df["Reason"] != "Global Noise Too High"]
        if not chrom_exclusions.empty:
            summary_stats = chrom_exclusions.groupby("Chromosome").size().sort_values(ascending=False)
            logger.info("Local Chromosome Exclusions (across valid samples):")
            for chrom, count in summary_stats.items():
                logger.info(f"  [-] {chrom}: Excluded in {count} sample(s)")
        else:
            logger.info("Local Chromosome Exclusions: None detected. All valid targets included.")
    else:
         logger.info("Local Chromosome Exclusions: None detected.")
    logger.info("==========================")
    
    if args.report:
        report_path = args.report
        if os.path.isdir(report_path):
            report_path = os.path.join(report_path, "evaluation_report.tsv")
        
        report_dir = os.path.dirname(os.path.abspath(report_path))
        if report_dir and not os.path.exists(report_dir):
            os.makedirs(report_dir)
            
        report_df.to_csv(report_path, sep='\t', index=False)
        logger.info(f"Evaluation report saved to: {report_path}")

    logger.info("--- Phase 3: Building Reference Matrices ---")
    logger.info("[TARGETS]")
    ref_targets = build_reference_matrix(targets_map, inclusion_map, good_ids)
    logger.info("[ANTITARGETS]")
    ref_antitargets = build_reference_matrix(antitargets_map, inclusion_map, good_ids)

    logger.info("--- Phase 4: Merging and Saving ---")
    ref_targets.add(ref_antitargets)
    ref_targets.sort()
    
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    logger.info(f"Writing output via Pandas to avoid TabIO corruptions: {args.output}")
    ref_targets.data.to_csv(args.output, sep='\t', index=False, float_format='%.6g')
    logger.info("Pipeline completed successfully.")

if __name__ == "__main__":
    main()