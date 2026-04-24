# Pan-gene Set Analysis of Candida auris

**Author:** Faiza Khan  
**Supervisor:** Prof. Andy Jones  
**Programme:** MSc Bioinformatics, University of Liverpool  
**Academic Year:** 2025/2026

## About This Project

This repository contains the code and analysis for my MSc research project on the pan-gene set of *Candida auris*, a multidrug-resistant fungal pathogen responsible for hospital outbreaks worldwide.

The project applies the GET_PANGENES pipeline to six *C. auris* genomes to construct a pan-gene presence/absence matrix, classify genes by conservation (cloud, shell, softcore, core) and analyse gene expression differences between fluconazole-sensitive and resistant strains.

## Repository Structure

```
cauris-pangene/
├── README.md
├── requirements.txt
├── pangene_occupancy_analysis.py   Main analysis script
└── data/                           Data folder (add your files here)
```

## Data Sources

1. GET_PANGENES pan-gene matrix (`cauris_get_pangenes_t0_run.xlsx`) produced by Luc Elliott on the University of Liverpool HPC cluster using the `-t 0` parameter recommended by Bruno Contreras-Moreira.
2. Reference gene summary from VEuPathDB (`Cauris_ref_genes_Summary_project.tsv`) containing B8441 annotations, InterPro domains, protein lengths and expression data.

## How to Run

1. Install dependencies:
   ```
   pip install pandas numpy matplotlib scipy openpyxl
   ```

2. Place the two data files in a `data/` folder beside the script.

3. Run the script:
   ```
   python pangene_occupancy_analysis.py --matrix data/cauris_get_pangenes_t0_run.xlsx --summary data/Cauris_ref_genes_Summary_project.tsv --outdir figures/
   ```

4. All plots will be saved into the `figures/` folder.

## Outputs

The script produces a combined six-panel figure and individual plots showing:
- Protein length by occupancy
- Percentage of pan-genes with InterPro domains
- Paralog count
- Gene expression (log2 TPM+1)
- Molecular weight
- Number of pan-gene clusters per occupancy class

All plots use Kruskal-Wallis tests for statistical comparison across occupancy classes.

## Key Findings

- 7,929 pan-gene clusters identified across 6 genomes
- Classic U-shaped pan-genome distribution: 31% cloud (unique) + 59% core
- Cauris6684 genome has approximately 1,930 extra predicted genes compared to the mean of the other five genomes, suggesting over-annotation
- Core genes show higher InterPro annotation rates (88%) and longer proteins than cloud genes

## References

- Contreras-Moreira et al. (2026) Genome Research 36:226 to 238
- Satoh et al. (2009) Microbiology & Immunology 53(1):41 to 44
- Rhodes & Fisher (2019) Current Opinion in Microbiology 52:84 to 91
- Hoff et al. (2019) Bioinformatics 35(21):4433 to 4435

## Acknowledgements

- Prof. Andy Jones (supervisor)
- Bruno Contreras-Moreira (GET_PANGENES pipeline)
- Luc Elliott (HPC cluster run)
- MSc Bioinformatics teaching team, University of Liverpool
