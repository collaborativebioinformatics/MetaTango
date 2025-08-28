# MetaTango

MetaTango is a pipeline designed to evaluate structural variant calling in pipelines in microbiomes. Specifically, MetaTango compares reference-based structural variant tools (Sniffles) and graph-based structural variant tools (Rhea) to determine which method is better in different situations. 

![MetaTango Logo](https://github.com/collaborativebioinformatics/MetaTango/blob/main/img/metatango_logo_v1.png)

### Background
Studying bacterial genome dynamics is critical for understanding the mechanisms of bacterial growth, adaptation, and phenotypic expression. Structural variants, large genomic rearrangements of 50 basepairs or more. 

### Workflow
![MetaTango Workflow](https://github.com/collaborativebioinformatics/MetaTango/blob/main/img/MetaTango_Workflow_v1.png)
**Figure 1** Overview of MetaTango Workflow

### Methods

#### Simulated and Synthetic Data
In order to effectively evaluate the performance of graph-based and reference-based SV detection methods, we developed both simulated and synthetic datasets with known horizontal gene transfer events. For simulated data, we initially selected a set of 6 bacterial species (see Table 1) and simulated HGT events for each of them using hgtsim[CITE], resulting in a second set of post-HGT bacterial species. We then simulated Pacbio HiFi reads for longitudinal metagenomic samples containing each of the 6 species using pbsim3[CITE], with downstream timepoints having a larger percentage of reads derived from the genomes with simulated HGT events. The first sample contains reads derived solely from the original unmutated set, while each subsequent sample has an increasing proportion of reads derived from the mutated set (see Figure 2). 

| Bacterial Species                         | Percent abundance | HGT mutation rate |
|-------------------------------------------|-------------------|-------------------|
| *Escherichia coli* (ATCC 700926)          | 16.66%            | TBD               |
| *Staphylococcus epidermidis* (ATCC 12228) | 16.66%            | TBD               |
| *Streptococcus mutans* (ATCC 700610)      | 16.66%            | TBD               |
| *Porphyromonas gingivalis* (ATCC 33277)   | 16.66%            | TBD               |
| *Cereibacter sphaeroides* (ATCC 17029)    | 16.66%            | TBD               |
| *Neisseria meningitidis* (ATCC BAA-335)   | 16.66%            | TBD               |
**Table 1.** Bacterial species contained within simulated samples and their HGT mutation rate.

![MetaTango Simulation Overview](https://github.com/collaborativebioinformatics/MetaTango/blob/main/img/Metagnomic_simulation.png)
**Figure 2** Overview of simulation process

#### Initial taxonomic classification and reference selection
In order to perform reference-based SV detection and be able to anchor graph-based methods, it is necessary to initially identify a set of references. To do this, we ran Sylph[CITE] on our reads to get an initial true positive set of species. Then, for each species, we pulled the reference from NCBI to get a final set of reference genomes to map to.  

#### Reference-based SV detection
For reference-based SV detection, we used a common framework. First, we mapped all reads using Minimap2 to the set of reference genomes identified with Sylph. Then, we used Sniffles2 to call SVs for each reference and returned a VCF for each of them. 

#### Graph-based SV detection
Raw longitudinal samples were input into Rhea[CITE] for the graph based SV detection. Rhea first creates a co-assembly graph of all longitudinal samples before re-mapping each sample individually back to the graph. The software then detects logfold changes for certain nodes contained in structures that correspond to different forms of structural variants (insertions, deletions, tandem repeats, etc). Rhea does not return a VCF, but rather a list of structural variants and the edges responsible for those variants. Therefore, for each structural variants detected, we backtracked into the assembly graph to get the specific sequence of the structural variant. Because an assembly graph does not provide any positional information relative to a reference, we then mapped the sequences against the reference genomes using Minimap2[CITE] in order to gain start and end positions for each structural variant. If a SV did not map to any reference, we noted this as well.  


### Running MetaTango


### Results
