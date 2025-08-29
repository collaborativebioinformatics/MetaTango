
#!/usr/bin/env python3

import sys
import os
import subprocess

def is_fasta(file_path):
    return file_path.endswith((".fa", ".fasta", ".fna"))

def is_fastq(file_path):
    return file_path.endswith((".fq", ".fastq"))

def is_bam(file_path):
    return file_path.endswith(".bam")

def run_cmd(cmd):
    print(f"[CMD] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <reference.fa> <input.bam|input.fastq|input.fasta> <Sp. name>")
        sys.exit(1)

    ref = sys.argv[1]
    inp = sys.argv[2]
    sp = sys.argv[3]
    proc = "6"

    # setting process dirs
    run_cmd(["mkdir", "-p", "scratch_metatango_%s"%sp])
    run_cmd(["mkdir", "-p", "results"])
    # Check first input: reference.fasta
    if not is_fasta(ref):
        print(f"Error: Reference {ref} is not a FASTA file (.fa/.fasta/.fna)")
        sys.exit(1)

    # Handle 2nd input file
    query_file = None
    if is_bam(inp):
        query_file = inp + ".fasta"
        print(f"Converting BAM -> FASTA: {inp} -> {query_file}")
        #run_cmd(["samtools", "view", inp, "| awk '{print ">"$1"\n"$10}' > ", query_file])
        cmd = f"samtools view {inp} | awk '{{print \">\"$1\"\\n\"$10}}' > {query_file}"
        subprocess.run(cmd, shell=True, check=True)

    elif is_fasta(inp) or is_fastq(inp):
        query_file = inp
    else:
        print(f"Error: Input {inp} must be BAM, FASTA, or FASTQ")
        sys.exit(1)

    # Running minimap2 for alignment
    output_paf = f"results/{sp}_to_ref.paf"
    output_bam = f"results/{sp}_to_ref.sorted.bam"

    # temps
    sam_out = "scratch_metatango_%s/tmp.1"%sp
    bam_temp = "scratch_metatango_%s/tmp.2"%sp

    options = ["-ax", "asm20", "-t", proc, "--secondary=yes", "-N", "50", "-L", "--MD"]

    # initial paf generation
    run_cmd(["minimap2", ref, query_file, "-o", output_paf]) # could be helpful for downstream Sv visual if needed

    run_cmd(["minimap2"] + options + [ref, query_file, "-o", sam_out])
    run_cmd(["samtools", "view", "-bS", sam_out, "-o", bam_temp])
    run_cmd(["samtools", "sort", "-o", output_bam, bam_temp])
    run_cmd(["samtools", "index", output_bam])

    os.remove(sam_out)
    os.remove(bam_temp)

    sniffles_vcf = f"results/{sp}_sv_calls.vcf"
    # running sniffles
    snif_options = ["--minsupport", "3", "--mapq", "10", "--minsvlen", "50", "--genotype-ploidy", "1", "--allow-overwrite"]
    run_cmd(["sniffles"] + snif_options + ["-t", proc, "--input", output_bam, "--reference", ref, "--vcf", sniffles_vcf])

if __name__ == "__main__":
    main()

