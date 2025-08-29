import random
from Bio import SeqIO
import os

# Input paths
genes_fasta = "/home/Users/rl152/rhea/hgtsim/selected_genes.fasta"
recipient_genomes_dir = "/home/Users/rl152/rhea/hgtsim/genomes_fasta"
output_distribution_file = "/home/Users/rl152/rhea/hgtsim/distribution_3.txt"

# Extract gene IDs from the multifasta file (remove the GCAxxx part)
gene_ids = []
for record in SeqIO.parse(genes_fasta, "fasta"):
    gene_name = record.id
    gene_ids.append(gene_name)

# Extract recipient genome base names (without .fna extension)
recipient_genomes = []
for filename in os.listdir(recipient_genomes_dir):
    if filename.endswith(".fna"):
        genome_base = os.path.splitext(filename)[0]  # Removes .fna
        recipient_genomes.append(genome_base)

# Create a distribution file
with open(output_distribution_file, "w") as f:
    for genome in recipient_genomes:
        num_genes = 10
        selected_genes = random.sample(gene_ids, num_genes)
        f.write(f"{genome},{','.join(selected_genes)}\n")

print(f"Distribution file has been created at: {output_distribution_file}")

