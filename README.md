# MetaTango

MetaTango is a benchmarking pipeline for evaluating **structural variant (SV) calling** in microbiome datasets.  
It provides a direct comparison between **reference-based** methods (e.g., Sniffles2) and **graph-based** methods (e.g., Rhea).

---

![MetaTango Logo](https://github.com/collaborativebioinformatics/MetaTango/blob/main/img/metatango_logo_v1.png)


## Slides:
If you want to see our slides, check here:
[MetaTango Hackathon Slides](https://docs.google.com/presentation/d/1x4vynogMJxUEtn7epWnm7PmspWcWUX_RqV_KffVWkLY/edit?slide=id.g3791c9a4996_0_18#slide=id.g3791c9a4996_0_18)

## Background
Structural variants play a crucial role in microbial evolution and adaptation (https://pubmed.ncbi.nlm.nih.gov/21298028/). MetaTango simulates **time-series metagenomic communities** with controlled **Horizontal Gene Transfer (HGT)** events, enabling a direct point of comparison of reference-guided and graph-based metagenomic SV callers to show the pros and cons of each approach.

---


## Workflow
![MetaTango Workflow](https://github.com/collaborativebioinformatics/MetaTango/blob/main/img/MetaTango_Workflow_v1.png)
**Figure 1** Overview of MetaTango Workflow

## High-level overview of 3 days of hacking:
- [x] We generated simulated bacterial communities with controlled HGT events.  
- [x]  We simulated long-read sequencing (PacBio HiFi).  
- [x] We called SVs with:
   - [x] **Reference-based:** Minimap2 + Sniffles2  
   - [x] **Graph-based:** Rhea (assembly graph-based detection via MetaFlye and Co-Assembly graph analysis)  
- [x] We mapped graph-derived SVs back to recruited reference genomes or flag as novel.  
- [x]  We computed summary statistics and evaluate runtime performance
- [ ] To do: integrate everything into an end-to-end pipeline (nextflow)
- [ ] To do: package into conda or docker
- [ ] To do: run on cheese microbiome dataset (ripe with HGT)
---

## Simulated data overview

- **Species:** Six bacterial strains (e.g., *E. coli*, *S. epidermidis*, *S. mutans*, *P. gingivalis*, *C. sphaeroides*, *N. meningitidis*).  
- **Abundance:** Equal abundance (~16.6% each).  
- **HGT Simulation:** Introduced with `hgtsim`.  
- **Sequencing:** PacBio HiFi reads simulated with `pbsim3`.  
- **Longitudinal Data:** Gradual increase in HGT-derived reads across timepoints.  

## Methods

### Simulated and Synthetic Data
In order to effectively evaluate the performance of graph-based and reference-based SV detection methods, we developed both simulated and synthetic datasets with known horizontal gene transfer events. For simulated data, we initially selected a set of 6 bacterial species (see Table 1) and simulated HGT events for each of them using hgtsim[https://doi.org/10.7717/peerj.4015], resulting in a second set of post-HGT bacterial species. We then simulated Pacbio HiFi reads for longitudinal metagenomic samples containing each of the 6 species using pbsim3[https://doi.org/10.1093/nargab/lqac092], with downstream timepoints having a larger percentage of reads derived from the genomes with simulated HGT events. The first sample contains reads derived solely from the original unmutated set, while each subsequent sample has an increasing proportion of reads derived from the mutated set (see Figure 2). 

| Bacterial Species                         | Percent abundance | HGT mutation rate |
|-------------------------------------------|-------------------|-------------------|
| *Escherichia coli* (ATCC 700926)          | 16.66%            | 10%               |
| *Staphylococcus epidermidis* (ATCC 12228) | 16.66%            | 10%               |
| *Streptococcus mutans* (ATCC 700610)      | 16.66%            | 10%               |
| *Porphyromonas gingivalis* (ATCC 33277)   | 16.66%            | 10%               |
| *Cereibacter sphaeroides* (ATCC 17029)    | 16.66%            | 10%               |
| *Neisseria meningitidis* (PartJ-Nmeningitidis-RM8376)   | 16.66%            | 10%               |
**Table 1.** Bacterial species contained within simulated samples and their HGT mutation rate.

![MetaTango Simulation Overview](https://github.com/collaborativebioinformatics/MetaTango/blob/main/img/Metagnomic_simulation.png)
**Figure 2** Overview of simulation process

### Initial taxonomic classification and reference selection
In order to perform reference-based SV detection and be able to anchor graph-based methods, it is necessary to initially identify a set of references. To do this, we ran Sylph (Shaw & Yu, 2025) on our reads to get an initial true positive set of species. Then, for each species, we pulled the reference from NCBI to get a final set of reference genomes to map to.  

### Reference-based SV detection
For reference-based SV detection, we used a common framework. First, we mapped all reads using Minimap2 (Li 2018; Li 2021) to the set of reference genomes identified with Sylph. Then, we used Sniffles2 (Smolka et al. 2024) to call SVs for each reference and returned a VCF for each of them. 

### Graph-based SV detection
Raw longitudinal samples were input into Rhea (Curry et al., 2024) for the graph based SV detection. Rhea first creates a co-assembly graph of all longitudinal samples before re-mapping each sample individually back to the graph. The software then detects logfold changes for certain nodes contained in structures that correspond to different forms of structural variants (insertions, deletions, tandem repeats, etc). Rhea does not return a VCF, but rather a list of structural variants and the edges responsible for those variants. Therefore, for each structural variants detected, we backtracked into the assembly graph to get the specific sequence of the structural variant. Because an assembly graph does not provide any positional information relative to a reference, we then mapped the sequences against the reference genomes using Minimap2 in order to gain start and end positions for each structural variant. If a SV did not map to any reference, we noted this as well. In addition, we added an HGT detection feature which calls HGT when there is at least one different reference hit among an SV and its neighboring sequences. 


## Running MetaTango
Running Rhea pipeline:
```bash
python make_vcf_and_detect_hgts.py t0.fq t1.fq --refs_folder path/to/reference/genomes
```
If you would like to change Rhea settings, you can use the flag `--rhea_flags` before adding Rhea flags. 

## Example Results

**Table 1 – HGT detection Metrics**

| Method        | Precision | Recall | F1-score |
|---------------|-----------|--------|----------|
| Sniffles2     | 0.92      | 0.85   | 0.89     |
| Rhea          | 0.91      | 0.93   | 0.92     |

### References
1. Curry, K. D., Yu, F. B., Vance, S. E., Segarra, S., Bhaya, D., Chikhi, R., Rocha, E. P. C., & Treangen, T. J. (2024). Reference-free structural variant detection in microbiomes via long-read co-assembly graphs. Bioinformatics (Oxford, England), 40(Suppl 1), i58–i67.
   
2. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics (Oxford, England), 34(18), 3094–3100.

3. Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics (Oxford, England), 37(23), 4572–4574.

4. Shaw, J., & Yu, Y. W. (2025). Rapid species-level metagenome profiling and containment estimation with sylph. Nature Biotechnology, 43(8), 1348–1359.

5. Smolka, M., Paulin, L. F., Grochowski, C. M., Horner, D. W., Mahmoud, M., Behera, S., Kalef-Ezra, E., Gandhi, M., Hong, K., Pehlivan, D., Scholz, S. W., Carvalho, C. M. B., Proukakis, C., & Sedlazeck, F. J. (2024). Publisher Correction: Detection of mosaic and population-level structural variants with Sniffles2. Nature Biotechnology, 42(10), 1616.

6. Yukiteru Ono, Michiaki Hamada, Kiyoshi Asai, PBSIM3: a simulator for all types of PacBio and ONT long reads, NAR Genomics and Bioinformatics, 4(4), lqac092
