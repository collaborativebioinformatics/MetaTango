import random
from pathlib import Path
from Bio import SeqIO
from BCBio import GFF

# Folders with FASTA + GFF files
fasta_dir = Path("genomes_fasta")
gff_dir = Path("genomes_gff")
output_file = "selected_genes.fasta"

random.seed(42)  # reproducibility

records_out = []

def extract_cds_features(feature, rec):
    """Recursively extract CDS features (with sequence) from a GFF feature tree."""
    genes = []
    if feature.type == "CDS":
        seq = feature.extract(rec.seq)
        gene_id = None
        if "protein_id" in feature.qualifiers:
            gene_id = feature.qualifiers["protein_id"][0]
        elif "ID" in feature.qualifiers:
            gene_id = feature.qualifiers["ID"][0]
        elif "locus_tag" in feature.qualifiers:
            gene_id = feature.qualifiers["locus_tag"][0]
        if not gene_id:
            gene_id = f"gene_{len(genes)+1}"
        genes.append((gene_id, seq))
    # recurse into subfeatures
    for sub in feature.sub_features:
        genes.extend(extract_cds_features(sub, rec))
    return genes


for fasta_file in fasta_dir.glob("*.fna"):
    species_name = fasta_file.stem
    accession_with_version = species_name.split("_")[0] + "_" + species_name.split("_")[1]
    gff_file = gff_dir / f"{accession_with_version}.gff"

    if not gff_file.exists():
        print(f"⚠️ No GFF for {fasta_file.name} (expected {gff_file.name}), skipping...")
        continue

    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    genes = []
    with open(gff_file) as gff_handle:
        for rec in GFF.parse(gff_handle, base_dict=fasta_dict):
            for feature in rec.features:
                genes.extend(extract_cds_features(feature, rec))

    if len(genes) == 0:
        print(f"⚠️ No CDS found in {gff_file.name}, skipping...")
        continue

    chosen = random.sample(genes, min(10, len(genes)))
    for i, (gene_id, seq) in enumerate(chosen, 1):
        # clean up gene naming
        clean_gene_id = gene_id.replace("cds-", "")
        rec = SeqIO.SeqRecord(
            seq,
            id=f"{accession_with_version}_{clean_gene_id}",
            description=""  # no trailing genomic_geneX
        )
        records_out.append(rec)

SeqIO.write(records_out, output_file, "fasta")
print(f"✅ Saved {len(records_out)} sequences to {output_file}")

