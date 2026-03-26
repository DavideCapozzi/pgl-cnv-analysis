#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clinical Metadata Extractor for CNVkit Heatmaps
------------------------------------------------
Calculates expected biological sex from target coverage depth
and extracts tumor location via filename regex matching.
"""

import os
import re
import argparse
import logging
import pandas as pd
import numpy as np

# =============================================================================
# LOGGING SETUP
# =============================================================================
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger("MetadataExtractor")

# =============================================================================
# THRESHOLDS FOR SEX DETERMINATION (WES TARGETED)
# =============================================================================
MIN_Y_RATIO = 0.10   # Threshold to call Male (overcoming WES sparse Y baits)
MIN_X_RATIO = 0.75   # Threshold to differentiate normal Female vs Female with X-loss

def calculate_weighted_depth(df: pd.DataFrame, chrom_list: list) -> float:
    """Calculates length-weighted mean depth for a set of chromosomes."""
    sub_df = df[df['norm_chrom'].isin(chrom_list)].copy()
    if sub_df.empty or sub_df['depth'].sum() == 0:
        return 0.0
    
    lengths = sub_df['end'] - sub_df['start']
    weights = lengths / lengths.sum()
    return np.average(sub_df['depth'], weights=weights)

def infer_biological_sex(df: pd.DataFrame) -> str:
    """Infers sex using orthogonal X and Y depth normalization against autosomes."""
    # Standardize chromosome names to remove 'chr' prefix safely
    df['norm_chrom'] = df['chromosome'].astype(str).str.replace('chr', '', regex=False)
    
    autosomes = [str(i) for i in range(1, 23)]
    x_chroms = ['X']
    y_chroms = ['Y']
    
    auto_depth = calculate_weighted_depth(df, autosomes)
    if auto_depth == 0:
        return "Unknown"
        
    x_depth = calculate_weighted_depth(df, x_chroms)
    y_depth = calculate_weighted_depth(df, y_chroms)
    
    norm_x = x_depth / auto_depth
    norm_y = y_depth / auto_depth
    
    if norm_y > MIN_Y_RATIO:
        return "Male"
    elif norm_x > MIN_X_RATIO:
        return "Female"
    else:
        # Y is missing/noise, but X is lower than expected for XX.
        return "Female_LOH_X"

def extract_metadata(target_dir: str, output_path: str):
    """Scans for .targetcoverage.cnn files, extracts metadata, and saves to TSV."""
    target_dir = os.path.abspath(target_dir)
    records = []
    
    logger.info(f"Scanning directory for .targetcoverage.cnn files: {target_dir}")
    
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.endswith(".targetcoverage.cnn"):
                filepath = os.path.join(root, file)
                
                # Core ID must perfectly match the R script sample_id logic
                sample_id = file.replace(".targetcoverage.cnn", "")
                
                # Regex modified to capture leading alphabetical characters 
                # (e.g., extracts 'PTJ' from 'PTJ62-t', 'PV' from 'PV5-t')
                location_match = re.search(r'^([A-Za-zA-Z]+)', sample_id)
                tumor_loc = location_match.group(1).upper() if location_match else "Unknown"
                
                try:
                    df = pd.read_csv(filepath, sep='\t')
                    sex_call = infer_biological_sex(df)
                except Exception as e:
                    logger.error(f"Failed parsing {file}: {e}")
                    sex_call = "Error"
                
                records.append({
                    "sample_id": sample_id,
                    "sex": sex_call,
                    "tumor_location": tumor_loc
                })
                
    if not records:
        logger.warning("No .targetcoverage.cnn files found.")
        return
        
    df_meta = pd.DataFrame(records)
    
    # Ensure output directory exists
    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    df_meta.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Processed {len(records)} samples.")
    logger.info(f"Metadata exported successfully to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Extracts Clinical Metadata from CNVkit output")
    parser.add_argument("-i", "--input", required=True, help="Root directory containing sample subdirectories (.targetcoverage.cnn)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path (e.g., cohort_metadata.tsv)")
    args = parser.parse_args()
    
    extract_metadata(args.input, args.output)

if __name__ == "__main__":
    main()