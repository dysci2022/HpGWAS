# HpGWAS

Analysis scripts for the manuscript:

**Identification of Helicobacter pylori Variation Associated with Gastric Cancer Risk in High-Incidence Regions**

This repository contains code used for the bacterial genome-wide association and downstream population-genomic analyses of *Helicobacter pylori* variants associated with gastric cancer risk, with emphasis on the leuS744V/SNP1625711 signal described in the manuscript.

## Contents

```text
scripts/
  rmvp_gwas.R                    rMVP GWAS workflow and visualization
  lambda_gc.R                    Genomic inflation factor calculation for pyseer output
  snippy_core_snp_workflow.sh     Snippy/core SNP calling and phylogeny notes
  finestructure_workflow.sh       fineSTRUCTURE preparation and staged run commands
  snpeff_annotation_workflow.sh   snpEff database build and SNP annotation commands
  prepare_snpgenie_inputs.py      Convert core VCF/GTF inputs for SNPGenie
  run_snpgenie_formal.sh          Plus/minus strand SNPGenie dN/dS workflow
docs/
  original_code_checksums.sha256  Checksums of the source files in the local code folder
```

## Study Context

The manuscript integrates public *H. pylori* genome data, clinical isolate validation, fecal qPCR validation, a prospective nested case-control cohort, and experimental follow-up. The scripts in this repository cover the computational components used for:

- core SNP calling and phylogenetic preparation with Snippy, PLINK, bcftools and IQ-TREE;
- bacterial GWAS using pyseer outputs and rMVP sensitivity analysis;
- genomic-control inflation assessment from GWAS p values;
- population-structure analysis using fineSTRUCTURE;
- candidate SNP annotation with snpEff;
- coding-sequence diversity analysis with SNPGenie.

## Required Software

Install the tools relevant to the workflow you plan to run:

- R with `rMVP`, `bigmemory`, `data.table`, `CMplot`, and `rgl`
- Snippy
- bcftools
- PLINK
- IQ-TREE 2
- fineSTRUCTURE
- GNU parallel
- snpEff
- SNPGenie
- Python 3
- Perl

## Input Files

The scripts expect user-provided genome and metadata inputs, including:

- assembled genomes (`*.fna`) for Snippy;
- a reference genome in FASTA/GenBank format;
- a core VCF from Snippy or another bacterial variant-calling workflow;
- phenotype and covariate files for GWAS;
- pyseer association output for lambda GC calculation;
- GTF/FASTA/VCF inputs for SNPGenie.

Raw sequencing data, clinical metadata, and unpublished subject-level data are not included in this code repository.

## Usage Notes

All scripts expose path variables near the top. Edit those variables before running the workflows on a new system.

Example:

```bash
Rscript scripts/lambda_gc.R SNPlmm_PCA_GRM_allSNPs.txt
bash scripts/snippy_core_snp_workflow.sh
bash scripts/finestructure_workflow.sh
bash scripts/snpeff_annotation_workflow.sh
bash scripts/run_snpgenie_formal.sh /path/to/dnds /path/to/SNPGenie-master
```

The original local scripts contained machine-specific paths. Those paths were converted to configurable variables where possible while preserving the analysis logic.

## Citation

If you use this code, please cite the associated manuscript once it is available.

