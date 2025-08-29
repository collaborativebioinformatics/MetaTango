#!/usr/bin/env python3
import argparse
import subprocess
import sys

#example:- 
## truvari bench -b sim3_ground_truth.vcf.gz  -c t3.vcf.gz -o sniffles_t3_vs_sim3/ -f all_refs.fasta --sizemin 50 --pctseq 0 --pctsize 0.5 --typeignore --pick multi

def run_truvari_bench(base_vcf, comp_vcf, output_dir, ref_fasta):
    cmd = [
        "truvari", "bench",
        "-b", base_vcf,
        "-c", comp_vcf,
        "-o", output_dir,
        "-f", ref_fasta,
        "--sizemin", "50",
        "--pctseq", "0",
        "--pctsize", "0.5",
        "--typeignore",
        "--pick", "multi"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print("Truvari bench completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running Truvari bench: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'truvari' command not found. Ensure Truvari is installed and in your PATH.", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run Truvari bench for SV comparison")
    parser.add_argument("-b", "--base_vcf", required=True, help="Base VCF file (ground truth)")
    parser.add_argument("-c", "--comp_vcf", required=True, help="Comparison VCF file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    parser.add_argument("-f", "--ref_fasta", required=True, help="Reference FASTA file")
    
    args = parser.parse_args()
    
    run_truvari_bench(args.base_vcf, args.comp_vcf, args.output_dir, args.ref_fasta)

if __name__ == "__main__":
    main()
