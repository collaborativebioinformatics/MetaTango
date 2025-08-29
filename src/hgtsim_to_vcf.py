#!/usr/bin/env python3
"""
Convert HgtSIM insertion report to VCF 4.2 format for use as ground truth

python hgtsim_to_vcf.py --input sim3_Step_2_insertion_report.txt --output sim3_ground_truth.vcf --gff $(ls GFFs/*.gff)


"""

import sys
import argparse
from datetime import datetime
import os

def parse_insertion_report(file_path):
    """Parse HgtSIM insertion report file"""
    insertions = []
    
    with open(file_path, 'r') as f:
        # Skip header line
        header = f.readline().strip()
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                recipient_genome = fields[0]
                recipient_contig = fields[1] 
                break_position = fields[2]
                inserted_sequence = fields[3]
                
                # Skip lines with None values
                if break_position == "None" or inserted_sequence == "None":
                    continue
                
                try:
                    pos = int(break_position)
                    insertions.append({
                        'chrom': recipient_contig,
                        'pos': pos,
                        'inserted_seq_id': inserted_sequence,
                        'genome': recipient_genome
                    })
                except ValueError:
                    continue
    
    return insertions

def parse_gff_files(gff_file_paths):
    """Parse GFF files to extract gene/feature lengths"""
    feature_lengths = {}
    
    for gff_file in gff_file_paths:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    start = int(fields[3])
                    end = int(fields[4])
                    attributes = fields[8]
                    
                    length = end - start + 1
                    
                    # Extract gene/protein IDs from attributes
                    for attr in attributes.split(';'):
                        if attr.startswith('ID=') or attr.startswith('locus_tag=') or attr.startswith('protein_id='):
                            feature_id = attr.split('=')[1]
                            feature_lengths[feature_id] = length
    
    return feature_lengths

def get_insertion_length(seq_id, feature_lengths):
    """Get actual insertion length from GFF data or estimate if not found"""
    # Try exact match first
    if seq_id in feature_lengths:
        return feature_lengths[seq_id]
    
    # Try partial matches for complex IDs
    for feature_id, length in feature_lengths.items():
        if seq_id in feature_id or feature_id in seq_id:
            return length
    
    # Fallback to estimation if not found in GFF
    if "WP_" in seq_id or "RS" in seq_id:
        return 1000  # Gene estimate
    else:
        return 500   # Default estimate

def generate_vcf_header(reference_contigs=None):
    """Generate VCF 4.2 header for HgtSIM ground truth"""
    current_date = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    
    header_lines = [
        "##fileformat=VCFv4.2",
        "##source=HgtSIM_ground_truth",
        f'##command="python hgtsim_to_vcf.py --input Step_2_insertion_report.txt --output ground_truth.vcf"',
        f'##fileDate="{current_date}"',
    ]
    
    # Add contig lines if provided
    if reference_contigs:
        for contig in reference_contigs:
            header_lines.append(f"##contig=<ID={contig}>")
    
    # Add standard VCF metadata
    header_lines.extend([
        "##ALT=<ID=INS,Description=\"Insertion\">",
        "##ALT=<ID=DEL,Description=\"Deletion\">", 
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=INV,Description=\"Inversion\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Structural variation with precise breakpoints\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variation\">",
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting the structural variation\">",
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">",
        "##INFO=<ID=HGTSIM_SEQ,Number=1,Type=String,Description=\"HgtSIM inserted sequence identifier\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    ])
    
    return "\n".join(header_lines)

def create_insertion_vcf_record(insertion, record_id, feature_lengths):
    """Create VCF record for HgtSIM insertion"""
    chrom = insertion['chrom']
    pos = insertion['pos']
    ref = "N"  # Reference is unknown at insertion site
    alt = "<INS>"
    qual = "60"
    filter_field = "PASS"
    
    # Get actual insertion length from GFF data
    svlen = get_insertion_length(insertion['inserted_seq_id'], feature_lengths)
    end_pos = pos  # For insertions, END = POS
    
    # Create INFO field
    info_parts = [
        "PRECISE",
        "SVTYPE=INS",
        f"SVLEN={svlen}",
        f"END={end_pos}",
        "SUPPORT=100",  # High confidence since it's ground truth
        "AF=1.000",     # Assume homozygous insertion
        f"HGTSIM_SEQ={insertion['inserted_seq_id']}"
    ]
    info = ";".join(info_parts)
    
    # FORMAT and sample fields
    format_field = "GT:GQ"
    sample_field = "1/1:60"  # Homozygous insertion, high quality
    
    vcf_id = f"HgtSIM.INS.{record_id}"
    
    return f"{chrom}\t{pos}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{sample_field}"

def convert_hgtsim_to_vcf(input_file, output_file, gff_files=None):
    """Main conversion function"""
    
    # Parse GFF files if provided
    feature_lengths = {}
    if gff_files:
        feature_lengths = parse_gff_files(gff_files)
        print(f"Parsed {len(feature_lengths)} features from {len(gff_files)} GFF files")
    
    # Parse HgtSIM insertion report
    insertions = parse_insertion_report(input_file)
    
    if not insertions:
        print("No valid insertions found in input file")
        return
    
    # Get unique contigs for header
    contigs = sorted(set(ins['chrom'] for ins in insertions))
    
    # Write VCF output
    with open(output_file, 'w') as f:
        # Write header
        f.write(generate_vcf_header(contigs) + "\n")
        
        # Write insertion records
        for i, insertion in enumerate(insertions, 1):
            vcf_record = create_insertion_vcf_record(insertion, i, feature_lengths)
            f.write(vcf_record + "\n")
    
    print(f"Converted {len(insertions)} insertions from {input_file} to {output_file}")
    print(f"Contigs found: {', '.join(contigs)}")
    
    # Report on length estimation
    found_in_gff = sum(1 for ins in insertions if get_insertion_length(ins['inserted_seq_id'], feature_lengths) != 1000)
    print(f"Found actual lengths for {found_in_gff}/{len(insertions)} insertions in GFF files")

def main():
    parser = argparse.ArgumentParser(description="Convert HgtSIM insertion report to VCF 4.2 format")
    parser.add_argument("--input", "-i", required=True, 
                       help="Input HgtSIM Step_2_insertion_report.txt file")
    parser.add_argument("--output", "-o", required=True, 
                       help="Output VCF file for ground truth")
    parser.add_argument("--gff", "-g", nargs='+', 
                       help="GFF files from source genomes to get actual feature lengths")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found")
        sys.exit(1)
    
    # Check GFF files exist if provided
    if args.gff:
        for gff_file in args.gff:
            if not os.path.exists(gff_file):
                print(f"Error: GFF file {gff_file} not found")
                sys.exit(1)
    
    convert_hgtsim_to_vcf(args.input, args.output, args.gff)

if __name__ == "__main__":
    main()
