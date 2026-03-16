# Somatic Variant Interpretation — Pan-Cancer Liquid Biopsy

## Project Overview
A clinical genomics case study performing somatic mutation interpretation 
on TCGA lung adenocarcinoma (LUAD) data, mirroring pan-cancer liquid biopsy 
workflows used in oncology diagnostics.

## Clinical Context
Somatic variant interpretation is central to liquid biopsy oncology panels 
(e.g. comprehensive pan-cancer profiling). This project demonstrates:
- Somatic variant annotation using COSMIC, OncoKB, and ClinVar
- AMP/ASCO/CAP tiered classification of actionable mutations
- BAM-level verification using IGV
- Clinical report drafting for key oncogenic variants

## Tools & Databases
- samtools 1.23, bcftools
- ANNOVAR, VEP
- COSMIC, OncoKB, ClinVar, cBioPortal
- Python 3.10, R 4.3
- IGV (Integrative Genomics Viewer)

## Dataset
TCGA Lung Adenocarcinoma (LUAD) — somatic MAF file from cBioPortal

## Project Structure
data/raw/          # Raw MAF files from TCGA/cBioPortal
data/annotated/    # Annotated variant files
results/           # Interpretation reports and figures
scripts/           # Analysis scripts
notebooks/         # Step-by-step analysis notebooks
