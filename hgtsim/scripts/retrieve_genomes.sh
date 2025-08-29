#!/bin/sh
# Make an output folder
mkdir -p genomes_fasta

# Download each genome as a zip
while read acc; do
    datasets download genome accession $acc --include genome --filename ${acc}.zip
    unzip -o ${acc}.zip -d ${acc}_unzipped
    mv ${acc}_unzipped/ncbi_dataset/data/*/*.fna genomes_fasta/
    rm -r ${acc}.zip ${acc}_unzipped
done < 6_species_genomes.txt

