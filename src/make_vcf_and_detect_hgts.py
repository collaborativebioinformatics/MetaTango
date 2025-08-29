import sys
import subprocess, random
import pandas as pd
import numpy as np
import os
import glob
import argparse


def run_rhea_with_fastas(fasta_files, rhea_flags=None):
    if len(fasta_files) < 2:
        sys.exit("Error: Please provide at least 2 FASTA/FASTQ files.")
    cmd = ["python", "./rhea.py"] + fasta_files
    if rhea_flags:
        cmd += rhea_flags  # append user-provided flags

    print(f"Running Rhea: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def process_structural_variant_file(sv_file, folder):
    os.makedirs(folder, exist_ok=True)
    out_file = os.path.join(folder, os.path.basename(sv_file))
    df = pd.read_csv(sv_file, sep='\t', header=0)
    NEW_HEADERS = ["node", "node_id", "node_length", "c0-SV_type", "c0-neighbor1", "c0-neighbor2", "c0-replacement"]
    df.columns = NEW_HEADERS
    df.to_csv(out_file, sep='\t', index=False)
    return out_file


def rhea_output_to_df(sv_file):
    df = pd.read_csv(sv_file, sep='\t')
    df = df.replace(np.nan, "Unspecified")
    vcf_df = pd.DataFrame({
        'CHROM': df['node'],
        'POS': df['node'],
        'ID': '.',
        'REF': 'Reference_sequence',
        'ALT': df['c0-SV_type'],
        'QUAL': '.',
        'FILTER': 'PASS',
        'INFO': "SVTYPE=" + df['c0-SV_type'].astype(str) + ";" +
                "SVLEN=" + df['node_length'].astype(str)
    })
    return df, vcf_df


def write_vcf(vcf_df, output_file):
    with open(output_file, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=Rhea\n")
        f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
        f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n")
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf_df.to_csv(f, sep="\t", index=False, header=False)


def extract_edges(gfa, out_fa):
    with open(gfa) as f, open(out_fa, "w") as out:
        for line in f:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                edge_id, seq = parts[1], parts[2]
                out.write(f">{edge_id}\n{seq}\n")


def integrate_hits_into_vcf(vcf_df, hits_out):

    hits_df = pd.read_csv(hits_out, sep="\t")

    hits_df["edge_id"] = hits_df["edge_id"].astype(str)
    vcf_df["CHROM"] = vcf_df["CHROM"].astype(str)

    merged_df = vcf_df.copy()

    # Replace CHROM and POS with ref_id and start if a hit exists
    merged_df = merged_df.merge(
        hits_df[["edge_id", "ref_id", "start"]],
        how="left",
        left_on="CHROM",
        right_on="edge_id"
    )

    # If hit not found, replace CHROM, POS with no_hits and not_found
    merged_df["CHROM"] = merged_df["ref_id"].fillna("no_hits")
    merged_df["POS"] = merged_df["start"].apply(
        lambda x: str(int(x)) if pd.notna(x) else "not_found"
    )

    merged_df = merged_df.drop(columns=["edge_id", "ref_id", "start"])

    # Drop unspecified SVs
    merged_df = merged_df[merged_df["ALT"] != "Unspecified"]

    return merged_df


def find_hgts(df, edge_hits_tsv):

    edge_hits = pd.read_csv(edge_hits_tsv, sep="\t")

    # Ensure consistent dtypes
    df["node"] = df["node"].astype(str)
    edge_hits["edge_id"] = edge_hits["edge_id"].astype(str)

    # Map node -> ref_id (via edge_id)
    df = df.merge(
        edge_hits[["edge_id", "ref_id", "start", "end", "mapq"]],
        how="left",
        left_on="node",
        right_on="edge_id"
    )

    # Build number -> ref_id mapping for neighbor/replacement
    edge_hits["edge_num"] = edge_hits["edge_id"].str.replace("edge_", "").astype(int)
    num_to_ref = dict(zip(edge_hits["edge_num"], edge_hits["ref_id"]))

    def map_neighbor(x):
        try:
            # Convert floats to int
            n = int(float(x))
            return num_to_ref.get(n, x)
        except (ValueError, TypeError):
            return x  # leave "Unspecified" etc. untouched

    for col in ["c0-neighbor1", "c0-neighbor2", "c0-replacement"]:
        if col in df.columns:
            df[col] = df[col].apply(map_neighbor)
            # Cleanup: replace leftover floats with "no_hit"
            df[col] = df[col].apply(
                lambda v: "no_hit" if isinstance(v, float) else v
            )

    # Filters
    df = df[df['c0-SV_type'] != 'Unspecified']
    df = df[df['c0-SV_type'] != 'tandem duplication gain']
    df = df[df['edge_id'].notna()]

    # Define HGT logic
    def is_hgt(row):
        svtype = row["c0-SV_type"]

        if svtype in ["deletion", "insertion"]:
            vals = [row.get("c0-neighbor1"), row.get("c0-neighbor2"), row.get("ref_id")]
            vals = [v for v in vals if v != "no_hit"]
            return len(set(vals)) > 1  # at least one different

        elif svtype in ["complex deletion", "complex insertion"]:
            vals = [
                row.get("c0-neighbor1"),
                row.get("c0-neighbor2"),
                row.get("c0-replacement"),
                row.get("ref_id"),
            ]
            vals = [v for v in vals if v != "no_hit"]
            return len(set(vals)) > 1  # at least one different

        else:
            return False

    df["hgt"] = df.apply(is_hgt, axis=1)

    # Keep only HGT rows
    df = df[df["hgt"]]
    df = df.rename(columns={"ref_id": "SV_source"})
    df = df.rename(columns={"start": "SV_start_position"})
    df["SV_start_position"] = df["SV_start_position"].astype("Int64")
    df = df[["node", "node_id", "node_length", "c0-SV_type", "c0-neighbor1", "c0-neighbor2", "c0-replacement", "SV_source", "SV_start_position"]]
    df = df.reset_index(drop=True)
    df.to_csv("detected_HGTs.tsv", sep="\t")

    return df

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run Rhea pipeline with user-provided refs folder.")
    parser.add_argument("--refs_folder", required=True, help="Path to folder with reference FASTA files (.fna/.fasta)")
    parser.add_argument("--rhea_flags", nargs="*", help="Optional flags to pass to Rhea")
    parser.add_argument("inputs", nargs="*", help="Input FASTA/FASTQ files for Rhea")
    args = parser.parse_args()

    if args.inputs:
        run_rhea_with_fastas(args.inputs, args.rhea_flags)

    # Prepare reference concatenation
    refs_folder = args.refs_folder
    refs_fa = "refs_combined.fa"
    refs_index = "refs.mmi"
    with open(refs_fa, "w") as outfile:
        for fname in os.listdir(refs_folder):
            if fname.endswith((".fna", ".fasta")):
                with open(os.path.join(refs_folder, fname)) as infile:
                    outfile.write(infile.read())
    subprocess.run(["minimap2", "-d", refs_index, refs_fa], check=True)

    # Loop over all structural variant TSV files
    sv_files = sorted(glob.glob("rhea_results/structual_variants-c*.tsv"))
    for sv_file in sv_files:
        x = sv_file.split("structual_variants-c")[-1].split(".tsv")[0]
        folder = f"T{x}"
        out_sv_file = process_structural_variant_file(sv_file, folder)  # replace headers of rhea output file
        os.chdir(folder)  # enter folder to run all downstream steps

        # Extract edges
        extract_edges("../rhea_results/metaflye/assembly_graph.gfa", "edges.fa")

        # Align edges to references
        sam_file = "edges.sam"
        bam_file = "edges.bam"
        sorted_bam = "edges.sorted.bam"
        hits_out = "edges_hits.tsv"
        subprocess.run(["minimap2", '-a', "../" + refs_index, "edges.fa"], stdout=open(sam_file, "w"), check=True)
        subprocess.run(["samtools", "view", "-bS", sam_file, "-o", bam_file], check=True)
        subprocess.run(["samtools", "sort", "-o", sorted_bam, bam_file], check=True)
        subprocess.run(["samtools", "index", sorted_bam], check=True)

        # Parse SAM/BAM to get best hits
        hits = {}
        view = subprocess.Popen(["samtools", "view", "-F", "0x904", sorted_bam],
                                stdout=subprocess.PIPE, text=True)
        for line in view.stdout:
            fields = line.strip().split("\t")
            qname, rname, pos, mapq, seq = fields[0], fields[2], int(fields[3]), int(fields[4]), fields[9]
            end = pos + len(seq) - 1
            if qname not in hits:
                hits[qname] = []
            hits[qname].append((mapq, rname, pos, end))

        best_hits = {}
        for qname, alignments in hits.items():
            max_mapq = max(a[0] for a in alignments)
            top_hits = [a for a in alignments if a[0] == max_mapq]
            best_hits[qname] = random.choice(top_hits)

        with open(hits_out, "w") as out:
            out.write("edge_id\tref_id\tstart\tend\tmapq\n")
            for qname, (mapq, rname, pos, end) in best_hits.items():
                out.write(f"{qname}\t{rname}\t{pos}\t{end}\t{mapq}\n")

        # Generate VCF and detect HGTs
        df, vcf_df = rhea_output_to_df(os.path.basename(out_sv_file))
        vcf_df = integrate_hits_into_vcf(vcf_df, hits_out)
        find_hgts(df, hits_out)
        write_vcf(vcf_df, "rhea_output_with_edge_hits.vcf")

        os.chdir("..")  # return to parent folder
