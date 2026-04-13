# DNA Variant Analysis — Project 1
**Bioinformatics Portfolio | University of Cape Coast, Ghana**

## Author
Stephen Opoku Asante
University of Cape Coast, Cape Coast, Ghana

## Overview
End-to-end DNA variant analysis pipeline identifying
disease-associated genetic variants in a simulated
exome sequencing cohort of 20 samples.

## Study Design
- 20 samples: 10 patients vs 10 healthy controls
- 500 simulated variants across 22 chromosomes
- 263 variants retained after QC filtering

## Pipeline Steps
1. Data simulation — realistic exome variant data
2. Quality filtering — PASS vs LowQual
3. Statistical QC — t-test for quality bias
4. Variant annotation — genes, consequences, clinical significance
5. Disease association — Fishers exact test + Bonferroni correction
6. Visualization — 8 publication-quality plots

## Key Findings
| Gene | Odds Ratio | Adjusted P | Significant |
|------|-----------|------------|-------------|
| HBA1 | >12 | <0.05 | YES |
| BRCA2 | ~11 | <0.05 | YES |
| PCSK9 | ~9 | <0.05 | YES |

## Clinical Relevance
HBA1 enrichment is consistent with high prevalence of
haemoglobinopathies in West African populations (Ghana).
G6PD deficiency affects up to 25% of people in Ghana.

## Tools Used
- R 4.5.2
- tidyverse, ggplot2, dplyr
- vcfR, BiocManager
- Bioconductor: VariantAnnotation, GenomicRanges

## Files
- Project1_vcf_raw.csv — raw simulated variants
- Project1_vcf_pass.csv — quality filtered variants
- Project1_vcf_annotated.csv — fully annotated variants
- Project1_association_results.csv — statistical results
- plots/ — all 8 publication figures

