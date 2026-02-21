#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNVkit Reference File Validator
-------------------------------
Performs rigorous structural and logical checks on a .cnn reference file.
Ensures compatibility with CNVkit segmentation and calling algorithms.

Usage:
    python validate_cnv_reference.py --input /path/to/reference.cnn

Author: Pipeline Optimization Team
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np

def print_header(title):
    print(f"\n{'=' * 60}")
    print(f" {title.upper()}")
    print(f"{'=' * 60}")

def check_status(name, passed, message=""):
    status = "PASS" if passed else "FAIL"
    # padding for alignment
    print(f"[{status:^4}] {name:<35} {message}")
    return passed

def validate_reference(file_path):
    results = []
    
    # 1. LOAD FILE
    print_header("Loading & Structural Check")
    if not os.path.exists(file_path):
        print(f"[FAIL] File not found: {file_path}")
        return False
        
    try:
        # CNVkit uses tab-separated files
        df = pd.read_csv(file_path, sep='\t')
        passed = True
        msg = f"Loaded {len(df)} bins"
    except Exception as e:
        passed = False
        msg = str(e)
    
    results.append(("File Load", passed))
    check_status("File Readability", passed, msg)
    if not passed: return False

    # 2. COLUMN VALIDATION
    required_cols = ['chromosome', 'start', 'end', 'log2', 'depth', 'spread', 'weight']
    missing_cols = [c for c in required_cols if c not in df.columns]
    
    passed = len(missing_cols) == 0
    msg = f"Missing: {missing_cols}" if not passed else "All required columns present"
    results.append(("Column Structure", passed))
    check_status("Header Structure", passed, msg)
    if not passed: return False # Critical failure

    # 3. NUMERICAL INTEGRITY (NaN/Inf)
    print_header("Numerical Integrity")
    
    # Check for NaNs
    nan_counts = df[required_cols].isna().sum()
    total_nans = nan_counts.sum()
    passed_nan = total_nans == 0
    msg_nan = f"Found {total_nans} NaNs" if not passed_nan else "No NaNs detected"
    results.append(("NaN Check", passed_nan))
    check_status("NaN Values", passed_nan, msg_nan)

    # Check for Infinite values
    numeric_cols = ['log2', 'depth', 'spread', 'weight']
    inf_counts = np.isinf(df[numeric_cols]).sum().sum()
    passed_inf = inf_counts == 0
    msg_inf = f"Found {inf_counts} Infinite values" if not passed_inf else "No Infinite values detected"
    results.append(("Inf Check", passed_inf))
    check_status("Infinite Values", passed_inf, msg_inf)

    # 4. COORDINATE LOGIC
    print_header("Genomic Logic")
    
    # Start < End
    invalid_coords = df[df['start'] >= df['end']]
    passed_coords = len(invalid_coords) == 0
    msg_coords = f"{len(invalid_coords)} bins with start >= end" if not passed_coords else "Coordinate order valid"
    results.append(("Coordinates", passed_coords))
    check_status("Coordinate Order", passed_coords, msg_coords)

    # Positive coordinates
    negative_coords = df[(df['start'] < 0) | (df['end'] < 0)]
    passed_neg = len(negative_coords) == 0
    msg_neg = f"{len(negative_coords)} bins with negative coords" if not passed_neg else "Non-negative coordinates"
    results.append(("Positivity", passed_neg))
    check_status("Coordinate Positivity", passed_neg, msg_neg)

    # 5. STATISTICAL DISTRIBUTION
    print_header("Statistical Distribution")

    # Depth > 0 (Warn only, 0 depth is technically possible but suspicious for reference)
    zero_depth = df[df['depth'] <= 0]
    passed_depth = len(zero_depth) == 0
    msg_depth = f"{len(zero_depth)} bins with <= 0 depth" if not passed_depth else "All depths positive"
    check_status("Depth Validity", passed_depth, msg_depth)
    # Note: We don't fail the pipeline for zero depth, but it is flagged.

    # Spread > 0 (Critical for weighting)
    # If spread is 0 or very close to 0, weights become infinite.
    low_spread = df[df['spread'] < 1e-5]
    passed_spread = len(low_spread) == 0
    msg_spread = f"{len(low_spread)} bins with near-zero spread" if not passed_spread else "Spread values valid"
    results.append(("Spread Validity", passed_spread))
    check_status("Spread Validity", passed_spread, msg_spread)

    # 6. FLAT REFERENCE LOGIC VERIFICATION
    print_header("Flat Reference Verification")
    
    # We expect a significant portion of bins to be exactly 0.0 (Fallback)
    # But not ALL bins (unless it's a completely flat ref, which is also valid but suspicious)
    
    flat_bins = df[df['log2'] == 0.0]
    n_flat = len(flat_bins)
    pct_flat = (n_flat / len(df)) * 100
    
    passed_flat = n_flat > 0
    msg_flat = f"{n_flat} bins ({pct_flat:.2f}%) are exactly 0.0"
    results.append(("Fallback Logic", passed_flat))
    check_status("Flat Fallback Detection", passed_flat, msg_flat)

    # 7. SUMMARY TABLE
    print_header("Validation Summary")
    
    print(f"{'CHECK NAME':<30} | {'STATUS':<10}")
    print("-" * 45)
    
    all_passed = True
    for name, result in results:
        status_str = "PASS" if result else "FAIL"
        print(f"{name:<30} | {status_str:<10}")
        if not result:
            all_passed = False

    print("-" * 45)
    
    if all_passed:
        print("\n[CONCLUSION] File is VALID and ready for CNVkit.")
        return True
    else:
        print("\n[CONCLUSION] File contains CRITICAL ERRORS. Do not use.")
        return False

def main():
    parser = argparse.ArgumentParser(description="Validate CNVkit .cnn reference file.")
    parser.add_argument("-i", "--input", required=True, help="Path to .cnn file")
    
    args = parser.parse_args()
    
    success = validate_reference(args.input)
    
    if success:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
