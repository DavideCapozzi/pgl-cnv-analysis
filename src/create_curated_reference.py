#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNVkit Curated Flat Reference Builder (Stable Release)
------------------------------------------------------
Strategy: "Signal Isolation with Flat Fallback"
Context:  Paraganglioma (PGL) Tumor-Only WES analysis.

Changes in this version:
  1. FIXED ImportError: Changed 'from cnvlib import tabio' to 'import cnvlib.tabio' 
     to correctly load the submodule in cnvkit 0.9.10 environments.
  2. FIXED AttributeError: Uses 'tabio.write()' for saving output.
  3. OPTIMIZATION: Suppresses expected RuntimeWarnings for NaN slices.
  4. CONFIG: Updated sample map (PTJ62-t).

Environment: numpy<2.0, pandas<2.0, cnvkit 0.9.10
"""

import os
import sys
import argparse
import logging
import warnings
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
# NOTE: Updated keys based on user correction (PTJ62-t instead of PT162-t)

SAMPLE_CHROMOSOME_MAP = {
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

    "PTJ62-t": [
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
    as the template.
    """
    if len(template_df) != len(new_df):
        raise ValueError(f"Bin count mismatch: Template has {len(template_df)}, {sample_name} has {len(new_df)}")
    
    # Check alignment of the first and last bin to catch shifts
    if (template_df.iloc[0]['start'] != new_df.iloc[0]['start']) or \
       (template_df.iloc[-1]['end'] != new_df.iloc[-1]['end']):
        raise ValueError(f"Coordinate mismatch in {sample_name}. Input files must define identical bins.")

def build_reference_matrix(file_dict, inclusion_map):
    """
    Vectorized construction of the reference profile.
    Uses masks to include only clean chromosomes per sample.
    """
    if not file_dict:
        raise ValueError("No input files provided to build_reference_matrix.")

    sample_ids = list(file_dict.keys())
    
    # 1. Load Template (First Sample)
    template_cnv = cnvlib.read(file_dict[sample_ids[0]])
    template_df = template_cnv.data
    
    n_bins = len(template_df)
    n_samples = len(sample_ids)
    
    logger.info(f"Initializing matrix: {n_bins} bins x {n_samples} samples")

    # 2. Pre-allocate NumPy arrays (float32 saves RAM)
    mat_log2 = np.full((n_bins, n_samples), np.nan, dtype=np.float32)
    mat_depth = np.full((n_bins, n_samples), np.nan, dtype=np.float32)
    
    # Pre-calculate chromosome vector for the template
    template_chroms = np.array([normalize_chrom(c) for c in template_df['chromosome']])

    # 3. Vectorized Loading & Masking
    valid_samples_count = 0
    
    for i, sample_id in enumerate(sample_ids):
        try:
            curr_cnv = cnvlib.read(file_dict[sample_id])
            validate_compatibility(template_df, curr_cnv.data, sample_id)
            
            allowed = set([normalize_chrom(c) for c in inclusion_map.get(sample_id, [])])
            if not allowed:
                logger.warning(f"Sample {sample_id} has no allowed chromosomes in map. Skipping.")
                continue

            # Create boolean mask for rows belonging to allowed chromosomes
            keep_mask = np.isin(template_chroms, list(allowed))
            
            raw_log2 = curr_cnv.data['log2'].values
            raw_depth = curr_cnv.data['depth'].values
            
            # Fill the matrix only where mask is True (other rows remain NaN)
            mat_log2[keep_mask, i] = raw_log2[keep_mask]
            mat_depth[keep_mask, i] = raw_depth[keep_mask]
            
            valid_samples_count += 1
            
        except Exception as e:
            logger.error(f"Error processing {sample_id}: {e}")

    if valid_samples_count == 0:
        raise RuntimeError("No valid samples were processed. Check your paths and Inclusion Map.")

    # 4. Statistical Reduction (Row-wise)
    logger.info("Computing reference statistics...")
    
    # Use context manager to strictly suppress "Mean of empty slice" warnings.
    # Empty slices are EXPECTED here (regions where no sample is clean) and handled below.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        
        ref_log2 = np.nanmean(mat_log2, axis=1)
        ref_depth = np.nanmean(mat_depth, axis=1)
        ref_spread = np.nanstd(mat_log2, axis=1)

    # 5. Flat Fallback Logic
    # Identify bins where ALL samples were masked (NaN)
    fallback_mask = np.isnan(ref_log2)
    n_fallback = np.sum(fallback_mask)
    
    if n_fallback > 0:
        percent = (n_fallback / n_bins) * 100
        logger.warning(f"Fallback triggered for {n_fallback} bins ({percent:.2f}%).")
        
        # A. Set Log2 to 0.0 (Neutral/Flat)
        ref_log2[fallback_mask] = 0.0
        
        # B. Set Depth to Global Mean (estimated from valid regions)
        global_mean_depth = np.nanmean(ref_depth)
        ref_depth[fallback_mask] = global_mean_depth if not np.isnan(global_mean_depth) else 1.0
        
        # C. Set Spread to Global Mean Spread (estimated from valid regions)
        global_mean_spread = np.nanmean(ref_spread)
        safe_spread = global_mean_spread if not np.isnan(global_mean_spread) else 0.1
        
        # Apply spread to fallback regions and any other NaNs
        ref_spread[fallback_mask] = safe_spread
        ref_spread[np.isnan(ref_spread)] = safe_spread

    # 6. Construct Final Object
    # We copy the template structure and replace data columns
    final_cnv = template_cnv.copy()
    final_cnv.data['log2'] = ref_log2
    final_cnv.data['depth'] = ref_depth
    final_cnv.data['spread'] = ref_spread
    
    # Recalculate Weights (Inverse variance weighting: weight ~ 1 / spread^2)
    # Clip spread to avoid division by zero
    clipped_spread = np.maximum(ref_spread, 1e-4)
    final_cnv.data['weight'] = 1.0 / (clipped_spread ** 2)
    
    return final_cnv

def locate_files(base_dir, sample_ids):
    """
    Finds target/antitarget files for the given samples.
    Performs robust path checking for flat or nested structures.
    """
    targets = {}
    antitargets = {}
    
    base_dir = os.path.abspath(base_dir)
    
    for sid in sample_ids:
        # Paths to check (Explicit list for debug logging)
        paths_checked = []

        # 1. Subfolder check: base_dir/SAMPLE/SAMPLE.targetcoverage.cnn
        t_path_1 = os.path.join(base_dir, sid, f"{sid}.targetcoverage.cnn")
        a_path_1 = os.path.join(base_dir, sid, f"{sid}.antitargetcoverage.cnn")
        paths_checked.append(f"Subfolder: {t_path_1}")

        # 2. Flat check: base_dir/SAMPLE.targetcoverage.cnn
        t_path_2 = os.path.join(base_dir, f"{sid}.targetcoverage.cnn")
        a_path_2 = os.path.join(base_dir, f"{sid}.antitargetcoverage.cnn")
        paths_checked.append(f"Flat Dir:  {t_path_2}")
        
        if os.path.exists(t_path_1) and os.path.exists(a_path_1):
            targets[sid] = t_path_1
            antitargets[sid] = a_path_1
        elif os.path.exists(t_path_2) and os.path.exists(a_path_2):
            targets[sid] = t_path_2
            antitargets[sid] = a_path_2
        else:
            logger.warning(f"Files not found for sample: {sid}")
            logger.debug(f"Checked paths:\n  " + "\n  ".join(paths_checked))
            
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

    # 1. Locate Files
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
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # 3. Process Antitargets
    logger.info("--- Processing Antitargets ---")
    try:
        ref_antitargets = build_reference_matrix(antitargets_map, SAMPLE_CHROMOSOME_MAP)
    except Exception as e:
        logger.critical(f"Failed during Antitarget processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # 4. Merge and Save
    logger.info("--- Merging and Saving ---")
    
    # In CNVkit, we add antitargets to targets to make the full reference set
    ref_targets.add(ref_antitargets)
    ref_targets.sort()
    
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # WRITING STEP
    try:
        logger.info(f"Writing curated reference (via Pandas) to: {args.output}")
        
        ref_targets.data.to_csv(
            args.output, 
            sep='\t',        
            index=False,    
            float_format='%.6g' 
        )
        
        logger.info("-" * 50)
        logger.info("SUCCESS. Curated Reference created.")
        logger.info("-" * 50)
    except Exception as e:
        logger.critical(f"Error saving file via Pandas: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()