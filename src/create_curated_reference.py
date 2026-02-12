#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNVkit Curated Flat Reference Builder (Final Optimized Version)
---------------------------------------------------------------
Strategy: "Signal Isolation with Flat Fallback"
Context:  Paraganglioma (PGL) Tumor-Only WES analysis.

Logic:
  1.  Loads target/antitarget coverage files (.cnn).
  2.  Validates that all files share the exact same genomic bins (preventing mismatch).
  3.  Applies an "Inclusive Explicit" mask: data is only used if the chromosome 
      is explicitly listed in SAMPLE_CHROMOSOME_MAP.
  4.  Computes weighted statistics (Log2, Spread, Depth) for the reference.
  5.  Fallback: If a bin has 0 valid samples, forces Log2=0 (Flat) and infers Spread.

Environment: numpy<2.0, pandas<2.0, cnvkit 0.9.10
"""

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
import cnvlib

# =============================================================================
# LOGGING SETUP
# =============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(asctime)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger("RefBuilder")

# =============================================================================
# CONFIGURATION: SAMPLE INCLUSION MAP
# =============================================================================
# STRICT RULE: Keys must match the sample name (folder/filename prefix).
# VALUES: List of chromosomes deemed "diploid/flat" by visual inspection.

SAMPLE_CHROMOSOME_MAP = {
    # --- Selected samples and chromosomes to filter noise, based on flat reference output ---
     
    "PT15-t": [
        "chr1", "chr2", "chr3", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
        "chr11", "chr12", "chr13", "chr15", "chr16", "chr17", "chr18", 
        "chr19", "chr20", "chr21", "chr22"
    ],

    "PTJ50-t": [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        "chr10", "chr11", "chr12", "chr13", "chr15", "chr16", "chr17", "chr18", 
        "chr19", "chr20", "chr21", "chr22"
    ],

    "PT20-t": [
        "chr2", "chr3", "chr4", "chr6", "chr8", "chr9", 
        "chr10", "chr12", "chr13", "chr15", "chr16", 
        "chr18", "chr19", "chr20", "chr21"
    ],

    "PT32-t": [
        "chr2", "chr4", "chr5", "chr7", "chr8", "chr9", 
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
          "chr20", "chr21"
    ],

    "PT61-t": [
        "chr2", "chr3", "chr4", "chr5", "chr6", "chr9", "chr10", 
        "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
        "chr19", "chr20", "chr21"
    ],

    "PT162-t": [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        "chr10", "chr12", "chr13", "chr16",  
        "chr20", "chr21", 
    ],
    
    "PTJ147-t": [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
        "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"
    ]
}

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def normalize_chrom(name):
    """Ensures consistent 'chr' prefix for comparisons."""
    s = str(name)
    return s if s.startswith("chr") else f"chr{s}"

def validate_compatibility(template_df, new_df, sample_name):
    """
    Critical Check: Ensures the new sample uses the exact same target kit
    as the template. Checks length and coordinate boundaries.
    """
    if len(template_df) != len(new_df):
        raise ValueError(f"Bin count mismatch: Template has {len(template_df)}, {sample_name} has {len(new_df)}")
    
    # Check alignment of the first and last bin to catch shifts
    # Using .iloc for integer-location based indexing (safe in pandas < 2.0)
    if (template_df.iloc[0]['start'] != new_df.iloc[0]['start']) or \
       (template_df.iloc[-1]['end'] != new_df.iloc[-1]['end']):
        raise ValueError(f"Coordinate mismatch in {sample_name}. Input files must define identical bins.")

def build_reference_matrix(file_dict, inclusion_map):
    """
    Vectorized construction of the reference profile.
    
    Args:
        file_dict (dict): Map of {sample_id: filepath}
        inclusion_map (dict): Map of {sample_id: [allowed_chromosomes]}
        
    Returns:
        cnvlib.CopyNumArray: The constructed reference object.
    """
    if not file_dict:
        raise ValueError("No input files provided to build_reference_matrix.")

    sample_ids = list(file_dict.keys())
    
    # 1. Load Template (First Sample)
    # We use this to preserve GC, RepeatMasker, Gene names, etc.
    template_cnv = cnvlib.read(file_dict[sample_ids[0]])
    template_df = template_cnv.data
    
    n_bins = len(template_df)
    n_samples = len(sample_ids)
    
    logger.info(f"Initializing matrix: {n_bins} bins x {n_samples} samples")

    # 2. Pre-allocate NumPy arrays for speed (float32 saves RAM)
    # We use NaN to represent "Masked/Ignored" data
    mat_log2 = np.full((n_bins, n_samples), np.nan, dtype=np.float32)
    mat_depth = np.full((n_bins, n_samples), np.nan, dtype=np.float32)
    
    # Pre-calculate chromosome vector for the template (used for masking)
    # Convert to standard format once
    template_chroms = np.array([normalize_chrom(c) for c in template_df['chromosome']])

    # 3. Vectorized Loading & Masking
    valid_samples_count = 0
    
    for i, sample_id in enumerate(sample_ids):
        try:
            # Load file
            curr_cnv = cnvlib.read(file_dict[sample_id])
            
            # Validation
            validate_compatibility(template_df, curr_cnv.data, sample_id)
            
            # Get allowed chromosomes
            allowed = set([normalize_chrom(c) for c in inclusion_map.get(sample_id, [])])
            
            if not allowed:
                logger.warning(f"Sample {sample_id} has no allowed chromosomes in map. Skipping.")
                continue

            # --- THE OPTIMIZATION ---
            # Create a boolean mask for the whole genome at once
            # True = Keep, False = Mask
            keep_mask = np.isin(template_chroms, list(allowed))
            
            # Extract data
            raw_log2 = curr_cnv.data['log2'].values
            raw_depth = curr_cnv.data['depth'].values
            
            # Inject into matrix (only where mask is True)
            # We assume the file is sorted identically (validated above)
            mat_log2[keep_mask, i] = raw_log2[keep_mask]
            mat_depth[keep_mask, i] = raw_depth[keep_mask]
            
            valid_samples_count += 1
            
        except Exception as e:
            logger.error(f"Error processing {sample_id}: {e}")

    if valid_samples_count == 0:
        raise RuntimeError("No valid samples were processed. Check your paths and Inclusion Map.")

    # 4. Statistical Reduction (Row-wise)
    logger.info("Computing reference statistics...")
    
    # np.nanmean / np.nanstd ignores the NaNs (masked regions)
    # This automatically implements the "Signal Isolation"
    with np.errstate(invalid='ignore'): # Suppress warnings for all-NaN rows
        ref_log2 = np.nanmean(mat_log2, axis=1)
        ref_depth = np.nanmean(mat_depth, axis=1)
        ref_spread = np.nanstd(mat_log2, axis=1)

    # 5. Flat Fallback Logic
    # Identify bins where ALL samples were masked (Result is NaN)
    fallback_mask = np.isnan(ref_log2)
    n_fallback = np.sum(fallback_mask)
    
    if n_fallback > 0:
        percent = (n_fallback / n_bins) * 100
        logger.warning(f"Fallback triggered for {n_fallback} bins ({percent:.2f}%).")
        
        # A. Set Log2 to 0.0 (Neutral/Flat)
        ref_log2[fallback_mask] = 0.0
        
        # B. Set Depth to Global Mean (Avoids 0 or 1 bias)
        global_mean_depth = np.nanmean(ref_depth)
        ref_depth[fallback_mask] = global_mean_depth if not np.isnan(global_mean_depth) else 1.0
        
        # C. Set Spread (Crucial for Segmentation Weighting)
        # We use the global average spread for fallback regions
        global_mean_spread = np.nanmean(ref_spread)
        safe_spread = global_mean_spread if not np.isnan(global_mean_spread) else 0.1
        ref_spread[fallback_mask] = safe_spread
        # Also fix any single-sample bins that yielded NaN std
        ref_spread[np.isnan(ref_spread)] = safe_spread

    # 6. Construct Final Object
    final_cnv = template_cnv.copy()
    final_cnv.data['log2'] = ref_log2
    final_cnv.data['depth'] = ref_depth
    final_cnv.data['spread'] = ref_spread
    
    # Recalculate Weights
    # CNVkit relies on 'weight' to penalize noisy bins.
    # Formula approximation: weight ~ 1 / variance
    # We clip spread to avoid division by zero
    clipped_spread = np.maximum(ref_spread, 1e-4)
    final_cnv.data['weight'] = 1.0 / (clipped_spread ** 2)
    
    return final_cnv

def locate_files(base_dir, sample_ids):
    """Finds target/antitarget files for the given samples."""
    targets = {}
    antitargets = {}
    
    # Normalize path
    base_dir = os.path.abspath(base_dir)
    
    for sid in sample_ids:
        # Search strategy: Check explicit subfolder, then flat dir
        # Expected pattern: {sid}.targetcoverage.cnn
        
        # 1. Subfolder check
        t_path = os.path.join(base_dir, sid, f"{sid}.targetcoverage.cnn")
        a_path = os.path.join(base_dir, sid, f"{sid}.antitargetcoverage.cnn")
        
        if not (os.path.exists(t_path) and os.path.exists(a_path)):
            # 2. Flat check
            t_path = os.path.join(base_dir, f"{sid}.targetcoverage.cnn")
            a_path = os.path.join(base_dir, f"{sid}.antitargetcoverage.cnn")
        
        if os.path.exists(t_path) and os.path.exists(a_path):
            targets[sid] = t_path
            antitargets[sid] = a_path
        else:
            logger.warning(f"Files not found for sample: {sid}")
            
    return targets, antitargets

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="CNVkit Curated Flat Reference Builder")
    parser.add_argument("-i", "--input", required=True, help="Path to directory containing .cnn files")
    parser.add_argument("-o", "--output", required=True, help="Output path for the curated reference.cnn")
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input):
        logger.error(f"Input directory does not exist: {args.input}")
        sys.exit(1)

    # 1. Locate Files based on the Config Map
    # We only look for samples defined in the MAP
    samples_to_find = list(SAMPLE_CHROMOSOME_MAP.keys())
    targets_map, antitargets_map = locate_files(args.input, samples_to_find)
    
    if not targets_map:
        logger.critical("No matching coverage files found. Check inputs and keys in script.")
        sys.exit(1)

    # 2. Process Targets
    logger.info("--- Processing Targets ---")
    try:
        ref_targets = build_reference_matrix(targets_map, SAMPLE_CHROMOSOME_MAP)
    except Exception as e:
        logger.critical(f"Failed during Target processing: {e}")
        sys.exit(1)

    # 3. Process Antitargets
    logger.info("--- Processing Antitargets ---")
    try:
        ref_antitargets = build_reference_matrix(antitargets_map, SAMPLE_CHROMOSOME_MAP)
    except Exception as e:
        logger.critical(f"Failed during Antitarget processing: {e}")
        sys.exit(1)

    # 4. Merge and Save
    logger.info("--- Merging and Saving ---")
    ref_targets.add(ref_antitargets)
    ref_targets.sort()
    
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    ref_targets.write(args.output)
    
    logger.info("-" * 50)
    logger.info(f"SUCCESS. Curated Reference saved to:\n{args.output}")
    logger.info("-" * 50)

if __name__ == "__main__":
    main()