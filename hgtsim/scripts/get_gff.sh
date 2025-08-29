#!/bin/sh
mkdir -p genomes_gff

while read acc; do
    # download gff
    datasets download genome accession $acc --include gff3 --filename ${acc}.zip
    
    # unzip to temp
    unzip -o ${acc}.zip -d ${acc}_unzipped
    
    # find the .gff file and rename it
    gff_file=$(find ${acc}_unzipped/ncbi_dataset/data/ -name "*.gff" | head -n 1)
    if [[ -f "$gff_file" ]]; then
        mv "$gff_file" "genomes_gff/${acc}.gff"
    else
        echo "⚠️ No GFF found for $acc"
    fi
    
    # cleanup
    rm -r ${acc}.zip ${acc}_unzipped
done < 6_species_genomes.txt

