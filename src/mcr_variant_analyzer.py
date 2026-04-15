#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNVkit Minimal Common Region (MCR) & Variant Analyzer (hg38 Edition)
--------------------------------------------------------------------
Purpose: Identifies overlapping genomic regions (MCRs) of concordant copy number 
alterations across multiple patients, extracts variants, and maps genes.
Excludes centromeric regions (hg38) to reduce mapping artifacts.
"""

import os
import sys
import argparse
import logging
import time
import requests
from typing import List, Dict, Tuple, Set

import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg') # Safe backend for headless server environments
import matplotlib.pyplot as plt

# =============================================================================
# LOGGING SETUP
# =============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger("MCR_Analyzer")

# =============================================================================
# HG38 CENTROMERE COORDINATES (NCBI GRC)
# =============================================================================
# Coordinates mapping based on GRCh38.p14 patches.
HG38_CENTROMERES = {
    'chr1': (122026460, 125184587),
    'chr2': (92188146, 94090557),
    'chr3': (90772459, 93655574),
    'chr4': (49708101, 51743951),
    'chr5': (46485901, 50059807),
    'chr6': (58553889, 59829934),
    'chr7': (58169654, 60828234),
    'chr8': (44033745, 45877265),
    'chr9': (43236168, 45518558),
    'chr10': (39686683, 41593521),
    'chr11': (51078349, 54425074),
    'chr12': (34769408, 37185252),
    'chr13': (16000001, 18051248),
    'chr14': (16000001, 18173523),
    'chr15': (17000001, 19725254),
    'chr16': (36311159, 38280682),
    'chr17': (22813680, 26885980),
    'chr18': (15460900, 20861206),
    'chr19': (24498981, 27190874),
    'chr20': (26436233, 30038348),
    'chr21': (10864561, 12915808),
    'chr22': (12954789, 15054318),
    'chrX': (58605580, 62412542),
    'chrY': (10316945, 10544039)
}

# =============================================================================
# UTILITIES & I/O
# =============================================================================

def normalize_chrom(name: str) -> str:
    s = str(name).strip()
    return s if s.startswith("chr") else f"chr{s}"

def chrom_sort_key(chrom: str) -> int:
    c = chrom.lower().replace("chr", "")
    if c.isdigit(): return int(c)
    if c in ['x']: return 23
    if c in ['y']: return 24
    if c in ['m', 'mt']: return 25
    return 99

def _extract_core_id(filename: str, file_type: str) -> str:
    suffixes = [".call.cns", ".cns"] if file_type == "cns" else [".hard-filtered.vcf.gz", ".hard-filtered.vcf", ".vcf.gz", ".vcf"]
    core_id = filename
    for suff in suffixes:
        if core_id.endswith(suff):
            core_id = core_id[:-len(suff)]
            break
    return core_id

def map_sample_files(cnv_dir: str, vcf_dir: str) -> Dict[str, Dict[str, str]]:
    cnv_dir, vcf_dir = os.path.abspath(cnv_dir), os.path.abspath(vcf_dir)
    sample_map = {}
    
    for root, _, files in os.walk(cnv_dir):
        for file in files:
            if file.endswith(".call.cns"):
                core_id = _extract_core_id(file, "cns")
                if core_id not in sample_map:
                    sample_map[core_id] = {'cns': None, 'vcf': None, 'sample_name': core_id}
                sample_map[core_id]['cns'] = os.path.join(root, file)
                
    for root, _, files in os.walk(vcf_dir):
        for file in files:
            if any(ext in file for ext in [".hard-filtered.vcf", ".vcf.gz", ".vcf"]):
                core_id = _extract_core_id(file, "vcf")
                if core_id in sample_map:
                    sample_map[core_id]['vcf'] = os.path.join(root, file)

    return {k: v for k, v in sample_map.items() if v['cns'] is not None}

# =============================================================================
# API GENE ANNOTATION (ENSEMBL REST)
# =============================================================================

def fetch_genes_for_region(chrom: str, start: int, end: int) -> str:
    clean_chrom = chrom.replace("chr", "")
    max_chunk = 4000000 
    genes_found = set()
    
    intervals = []
    curr_start = start
    while curr_start <= end:
        curr_end = min(curr_start + max_chunk, end)
        intervals.append((curr_start, curr_end))
        curr_start = curr_end + 1

    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}

    for c_start, c_end in intervals:
        ext = f"/overlap/region/human/{clean_chrom}:{c_start}-{c_end}?feature=gene"
        try:
            r = requests.get(server + ext, headers=headers, timeout=10)
            if r.status_code == 200:
                for gene in r.json():
                    name = gene.get('external_name')
                    if name: genes_found.add(name)
            time.sleep(0.1) 
        except Exception:
            continue

    return ",".join(sorted(list(genes_found))) if genes_found else "None_or_API_Error"

# =============================================================================
# CORE INTERVAL LOGIC & CENTROMERE FILTERING
# =============================================================================

def filter_centromeres(segments: List[Dict]) -> List[Dict]:
    """
    Substracts hg38 centromere coordinates from patient segments.
    Splits segments into two if they span across the entire centromere.
    """
    filtered = []
    for seg in segments:
        chrom = seg['chrom']
        start = seg['start']
        end = seg['end']
        
        if chrom not in HG38_CENTROMERES:
            filtered.append(seg)
            continue
            
        cen_start, cen_end = HG38_CENTROMERES[chrom]
        
        if end < cen_start or start > cen_end:
            filtered.append(seg)
        elif start >= cen_start and end <= cen_end:
            continue 
        elif start < cen_start and end <= cen_end:
            seg['end'] = cen_start - 1
            filtered.append(seg)
        elif start >= cen_start and end > cen_end:
            seg['start'] = cen_end + 1
            filtered.append(seg)
        elif start < cen_start and end > cen_end:
            seg1 = seg.copy()
            seg1['end'] = cen_start - 1
            filtered.append(seg1)
            
            seg2 = seg.copy()
            seg2['start'] = cen_end + 1
            filtered.append(seg2)
            
    # Protect against any invalid zero/negative length segments
    return [s for s in filtered if s['start'] <= s['end']]

def parse_and_merge_cnv(file_path: str, sample_id: str, chroms: List[str], min_log2: float) -> Tuple[List[Dict], List[Dict]]:
    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        logger.error(f"Failed to read {file_path}: {e}")
        return [], []

    df['chromosome'] = df['chromosome'].apply(normalize_chrom)
    df = df[df['chromosome'].isin(chroms)]
    
    if 'log2' in df.columns:
        df = df[df['log2'].abs() >= min_log2]
    
    amps = df[df['cn'] > 2].copy()
    dels = df[df['cn'] < 2].copy()
    
    def merge_segments(sub_df: pd.DataFrame, alt_type: str) -> List[Dict]:
        if sub_df.empty: return []
        sub_df = sub_df.sort_values(by=['chromosome', 'start'])
        merged = []
        for _, row in sub_df.iterrows():
            if not merged:
                merged.append({'chrom': row['chromosome'], 'start': row['start'], 'end': row['end'], 'sample': sample_id, 'type': alt_type})
            else:
                last = merged[-1]
                if row['chromosome'] == last['chrom'] and row['start'] <= last['end'] + 1:
                    last['end'] = max(last['end'], row['end'])
                else:
                    merged.append({'chrom': row['chromosome'], 'start': row['start'], 'end': row['end'], 'sample': sample_id, 'type': alt_type})
        return merged

    return merge_segments(amps, 'AMP'), merge_segments(dels, 'DEL')

def sweep_line_mcr(segments: List[Dict], min_patients: int) -> List[Dict]:
    if not segments: return []
        
    mcrs = []
    chrom_groups = {}
    for seg in segments:
        chrom_groups.setdefault(seg['chrom'], []).append(seg)
        
    for chrom, chrom_segments in chrom_groups.items():
        events = []
        for seg in chrom_segments:
            events.append((seg['start'], 'start', seg['sample']))
            events.append((seg['end'], 'end', seg['sample']))
            
        events.sort(key=lambda x: (x[0], 0 if x[1] == 'start' else 1))
        
        active_samples: Set[str] = set()
        cumulative_samples: Set[str] = set()
        current_mcr_start = None
        
        for pos, event_type, sample in events:
            if event_type == 'start':
                active_samples.add(sample)
                if current_mcr_start is not None:
                    cumulative_samples.add(sample)
            else:
                active_samples.discard(sample)
                
            current_count = len(active_samples)
            
            if current_count >= min_patients and current_mcr_start is None:
                current_mcr_start = pos
                cumulative_samples = set(active_samples)
                
            elif current_count < min_patients and current_mcr_start is not None:
                if pos > current_mcr_start:
                    mcrs.append({
                        'chrom': chrom,
                        'start': current_mcr_start,
                        'end': pos,
                        'n_patients': len(cumulative_samples),
                        'patients': ",".join(sorted(list(cumulative_samples)))
                    })
                current_mcr_start = None
                cumulative_samples = set()
                
    return mcrs

# =============================================================================
# VARIANT EXTRACTION & FILTERING
# =============================================================================

def extract_variants_in_mcr(mcr: Dict, sample_map: Dict, mcr_id: str, min_dp: int, min_af: float) -> List[Dict]:
    involved_samples = mcr['patients'].split(",")
    variants_found = []
    
    for sample in involved_samples:
        vcf_path = sample_map.get(sample, {}).get('vcf')
        if not vcf_path or not os.path.exists(vcf_path): continue
            
        try:
            vcf = pysam.VariantFile(vcf_path)
            
            # Dynamically resolve contig naming mismatch between CNVkit (chr) and VCF
            query_chrom = mcr['chrom']
            if query_chrom not in vcf.header.contigs:
                alt_chrom = query_chrom.replace("chr", "")
                if alt_chrom in vcf.header.contigs:
                    query_chrom = alt_chrom
                    
            records = vcf.fetch(query_chrom, mcr['start'], mcr['end'])
            for rec in records:
                var_dict = _parse_vcf_record(rec, sample, mcr_id, mcr['chrom'], min_dp, min_af)
                if var_dict: variants_found.append(var_dict)
        except ValueError:
            # Fallback for unindexed VCFs, applying contig normalization
            try:
                vcf.reset()
                alt_chrom_fb = mcr['chrom'].replace("chr", "")
                for rec in vcf:
                    if rec.chrom in [mcr['chrom'], alt_chrom_fb] and mcr['start'] <= rec.pos <= mcr['end']:
                        var_dict = _parse_vcf_record(rec, sample, mcr_id, mcr['chrom'], min_dp, min_af)
                        if var_dict: variants_found.append(var_dict)
            except Exception as inner_e:
                logger.warning(f"Failed linear scan for {sample} in {mcr_id}: {inner_e}")
        except Exception as e:
            logger.warning(f"Failed to process VCF for {sample}: {e}")
        finally:
            vcf.close()
        
    return variants_found

def _parse_vcf_record(rec: pysam.VariantRecord, sample: str, mcr_id: str, chrom: str, min_dp: int, min_af: float) -> Dict:
    dp_val, max_af_val = 0, 0.0
    if rec.samples:
        sample_data = rec.samples[0] 
        dp_val = sample_data.get('DP', 0) or 0
        af_raw = sample_data.get('AF', 0.0)
        
        # Robust handling for multiallelic variants (tuple) to avoid missing driver mutations
        if isinstance(af_raw, tuple) and af_raw:
            valid_afs = [float(v) for v in af_raw if v is not None]
            max_af_val = max(valid_afs) if valid_afs else 0.0
        else:
            max_af_val = float(af_raw) if af_raw is not None else 0.0
            
    if dp_val < min_dp or max_af_val < min_af: return {}
    
    return {
        'MCR_ID': mcr_id, 'Sample': sample, 'Chromosome': chrom,
        'Position': rec.pos, 'Ref': rec.ref, 'Alt': ",".join(rec.alts) if rec.alts else ".",
        'Filter': ";".join(rec.filter.keys()) if rec.filter.keys() else "PASS",
        'Depth': dp_val, 'Allele_Freq': round(max_af_val, 4)
    }

# =============================================================================
# REPORTING & VISUALIZATION
# =============================================================================

def generate_mcr_barplot(df_mcr: pd.DataFrame, out_dir: str):
    """Generates a bar plot representing patient counts for each MCR ordered by genome."""
    if df_mcr.empty: return
    
    # Create an internal numeric key for perfect sorting along the genome
    df_plot = df_mcr.copy()
    df_plot['chrom_num'] = df_plot['Chromosome'].str.replace('chr', '').replace({'X': 23, 'Y': 24, 'M': 25, 'MT': 25}).astype(int)
    df_plot = df_plot.sort_values(['chrom_num', 'Start_MB'])
    
    plt.figure(figsize=(16, 8))
    bars = plt.bar(
        df_plot['MCR_ID'], 
        df_plot['N_Patients'], 
        color=['#d62728' if t == 'AMP' else '#1f77b4' for t in df_plot['Type']]
    )
    
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.ylabel('Number of Patients', fontsize=12)
    plt.title('MCR Patient Involvement across hg38 Genome', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Legend handling
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#d62728', label='Amplification (AMP)'),
                       Patch(facecolor='#1f77b4', label='Deletion (DEL)')]
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plot_path = os.path.join(out_dir, "mcr_patient_counts_hg38.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    logger.info(f"[+] Saved MCR distribution plot to: {plot_path}")

# =============================================================================
# MAIN ORCHESTRATOR
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Extracts MCRs (Centromere excluded), maps genes, and filters variants.")
    parser.add_argument("--cnv-dir", required=True, help="Recursive root directory containing *.call.cns files.")
    parser.add_argument("--vcf-dir", required=True, help="Recursive root directory containing *.vcf or *.vcf.gz files.")
    parser.add_argument("--out-dir", required=True, help="Output directory for reports.")
    parser.add_argument("--chromosomes", required=True, help="Comma-separated list of chromosomes (e.g., 1,2,3).")
    parser.add_argument("--min-patients", type=int, default=3, help="Min samples to form an MCR (default: 3).")
    parser.add_argument("--min-log2", type=float, default=0.25, help="Min absolute log2 ratio (default: 0.25).")
    parser.add_argument("--min-dp", type=int, default=20, help="Min Depth (DP) for variant (default: 20).")
    parser.add_argument("--min-af", type=float, default=0.05, help="Min Allele Frequency (AF) for variant (default: 0.05).")
    parser.add_argument("--annotate-genes", action="store_true", help="Enable Ensembl API querying.")
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    
    log_file_path = os.path.join(args.out_dir, "mcr_analysis_run.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
    logger.addHandler(file_handler)
    
    target_chroms = [normalize_chrom(c.strip()) for c in args.chromosomes.split(",")]
    
    logger.info("--- Phase 1: File Mapping ---")
    sample_map = map_sample_files(args.cnv_dir, args.vcf_dir)
    if not sample_map:
        logger.error("No valid CNV files found. Exiting.")
        return

    logger.info(f"--- Phase 2: Interval Parsing, Centromere Filtering & MCR Detection ({args.min_patients}+ patients) ---")
    all_amps, all_dels = [], []
    
    for sample, paths in sample_map.items():
        amps, dels = parse_and_merge_cnv(paths['cns'], sample, target_chroms, args.min_log2)
        all_amps.extend(amps)
        all_dels.extend(dels)
        
    # Apply Centromere Filter before Sweep-line
    all_amps = filter_centromeres(all_amps)
    all_dels = filter_centromeres(all_dels)
        
    mcr_amps = sweep_line_mcr(all_amps, args.min_patients)
    mcr_dels = sweep_line_mcr(all_dels, args.min_patients)
    
    for m in mcr_amps: m['type'] = 'AMP'
    for m in mcr_dels: m['type'] = 'DEL'
    
    all_mcrs = mcr_amps + mcr_dels
    all_mcrs.sort(key=lambda x: (chrom_sort_key(x['chrom']), x['start']))
    
    logger.info(f"Detected {len(all_mcrs)} Minimal Common Regions ({len(mcr_amps)} AMP, {len(mcr_dels)} DEL).")
    if not all_mcrs: return

    logger.info(f"--- Phase 3: Gene Annotation & Variant Extraction (DP>={args.min_dp}, AF>={args.min_af}) ---")
    mcr_records, variant_records = [], []
    
    for i, mcr in enumerate(all_mcrs, 1):
        mcr_id = f"MCR_{mcr['chrom']}_{mcr['type']}_{i:03d}"
        annotated_genes = fetch_genes_for_region(mcr['chrom'], mcr['start'], mcr['end']) if args.annotate_genes else "Skipped"
        
        mcr_records.append({
            'MCR_ID': mcr_id,
            'Chromosome': mcr['chrom'],
            'Start_MB': round(mcr['start'] / 1_000_000.0, 6),
            'End_MB': round(mcr['end'] / 1_000_000.0, 6),
            'Length_MB': round((mcr['end'] - mcr['start']) / 1_000_000.0, 6),
            'Type': mcr['type'],
            'N_Patients': mcr['n_patients'],
            'Patient_IDs': mcr['patients'],
            'Gene_Symbols': annotated_genes
        })
        
        variant_records.extend(extract_variants_in_mcr(mcr, sample_map, mcr_id, args.min_dp, args.min_af))

    logger.info("--- Phase 4: Report Generation & Plotting ---")
    df_mcr = pd.DataFrame(mcr_records)
    out_mcr_path = os.path.join(args.out_dir, "mcr_regions_summary_hg38.xlsx")
    df_mcr.to_excel(out_mcr_path, index=False, engine='openpyxl')
    logger.info(f"[+] Saved MCR regions to: {out_mcr_path}")
    
    # Generate Plot
    generate_mcr_barplot(df_mcr, args.out_dir)
    
    if variant_records:
        df_var = pd.DataFrame(variant_records)
        out_var_path = os.path.join(args.out_dir, "mcr_patient_variants_hg38.xlsx")
        try:
            df_var.to_excel(out_var_path, index=False, engine='openpyxl')
            logger.info(f"[+] Saved {len(variant_records)} filtered variant annotations to: {out_var_path}")
        except Exception as e:
            logger.error(f"Failed to save variants. Error: {e}")
            
    logger.info("Analysis complete.")

if __name__ == "__main__":
    main()