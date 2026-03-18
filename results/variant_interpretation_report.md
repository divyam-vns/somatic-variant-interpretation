# Somatic Variant Interpretation Report
## TCGA Lung Adenocarcinoma (LUAD) — Pan-Cancer Liquid Biopsy Case Study

**Analyst:** Divya Mishra, Ph.D.  
**Date:** March 2026  
**Dataset:** TCGA-LUAD (n=513 patients, 196,952 PASS variants)  
**Guideline:** AMP/ASCO/CAP 2017 + ACMG/AMP framework  
**Tools:** CancerHotspots API, ClinVar, OncoKB (pending), GATK, IGV  

---

## Executive Summary

Somatic variant interpretation was performed on 513 TCGA lung adenocarcinoma 
(LUAD) patients using a clinical liquid biopsy-aligned workflow. A total of 
691 unique coding variants were identified across 17 clinically relevant genes. 
Using AMP/ASCO/CAP tiered classification:

| Tier | Count | Clinical Meaning |
|------|-------|-----------------|
| Tier I | 14 | Strong clinical significance — FDA-approved therapy exists |
| Tier II | 590 | Potential clinical significance — investigational or functional evidence |
| Tier IV | 87 | Unknown significance |

Key oncogenic drivers identified: KRAS (29.4%), TP53 (52.0%), EGFR (14.2%), 
KEAP1 (18.5%), STK11 (14.6%), BRAF (8.2%).

---

## Variant 1: KRAS p.G12C — Tier I

**Gene:** KRAS (Kirsten Rat Sarcoma Viral Proto-Oncogene)  
**Variant:** p.G12C (c.34G>T), Missense Mutation  
**Frequency in cohort:** 29.4% of LUAD patients (151/513)  
**Hotspot:** Yes — 2,175 tumors across cancer types (CancerHotspots)  
**VAF range:** 0.10–0.65 (median ~0.35)  

### Biological Mechanism
KRAS G12C substitutes glycine for cysteine at codon 12 in the P-loop of the 
GTPase domain. This locks KRAS in the GTP-bound active state, constitutively 
activating RAS/MAPK and PI3K/AKT signaling, driving uncontrolled proliferation.
The cysteine substitution uniquely creates a druggable covalent binding pocket.

### Clinical Significance
**FDA-approved therapies:**
- Sotorasib (Lumakras) — FDA approved May 2021 for KRAS G12C+ NSCLC 
  (CodeBreaK 100 trial, ORR 36%)
- Adagrasib (Krazati) — FDA approved December 2022 for KRAS G12C+ NSCLC
  (KRYSTAL-1 trial, ORR 43%)

**AMP/ASCO/CAP Classification: Tier I — Level A Evidence**  
Strong clinical significance with FDA-approved targeted therapy in NSCLC.

### Liquid Biopsy Considerations
At VAF <5%, KRAS G12C detection in cfDNA requires high-sensitivity assays 
(ddPCR or high-depth NGS). False negatives are possible in low tumor fraction 
samples — orthogonal confirmation recommended for VAFs <2%.

---

## Variant 2: EGFR p.L858R — Tier I

**Gene:** EGFR (Epidermal Growth Factor Receptor)  
**Variant:** p.L858R (c.2573T>G), Missense Mutation  
**Frequency in cohort:** 14.2% of LUAD patients (73/513)  
**Hotspot:** Yes — confirmed cancer hotspot  
**VAF range:** 0.15–0.70  

### Biological Mechanism
L858R substitution in EGFR kinase domain (exon 21) destabilizes the inactive 
conformation, locking the receptor in constitutively active state. This drives 
downstream RAS/MAPK and PI3K/AKT/mTOR signaling independently of ligand binding.

### Clinical Significance
**FDA-approved therapies (NSCLC):**
- Erlotinib, Gefitinib, Afatinib — 1st/2nd generation EGFR TKIs
- Osimertinib (Tagrisso) — 3rd generation, preferred first-line (FLAURA trial)
- Osimertinib also active against T790M resistance mutation

**AMP/ASCO/CAP Classification: Tier I — Level A Evidence**  
Standard of care biomarker for NSCLC — NCCN guidelines recommend EGFR testing 
for all advanced NSCLC patients.

### Co-mutation Context
EGFR L858R is typically mutually exclusive with KRAS mutations. Concurrent 
TP53 mutations are common and associated with reduced TKI response duration.

### Liquid Biopsy Considerations
EGFR L858R is one of the most validated cfDNA biomarkers. Sensitivity of 
liquid biopsy for EGFR mutations is 60-80% vs. tissue biopsy. Negative liquid 
biopsy does not rule out EGFR mutation — tissue confirmation recommended.

---

## Variant 3: BRAF p.V600E — Tier I

**Gene:** BRAF (B-Raf Proto-Oncogene)  
**Variant:** p.V600E (c.1799T>A), Missense Mutation  
**Frequency in cohort:** 8.2% of LUAD patients (42/513)  
**Hotspot:** Yes — 897 tumors (CancerHotspots)  

### Biological Mechanism
V600E substitution in the BRAF kinase activation loop constitutively activates 
BRAF kinase activity ~500-fold above baseline, bypassing RAS activation and 
driving ERK/MAPK signaling. Unlike other BRAF mutations, V600E signals as a 
monomer independent of RAS.

### Clinical Significance
**FDA-approved therapies (NSCLC):**
- Dabrafenib + Trametinib (BRAF + MEK inhibitor combination) — FDA approved 
  June 2017 for BRAF V600E+ metastatic NSCLC (ORR 64%, GEOMETRY trial)

**AMP/ASCO/CAP Classification: Tier I — Level A Evidence**

### Resistance Mechanisms
Primary resistance: concurrent RAS mutations, NF1 loss.  
Acquired resistance: KRAS amplification, MEK1/2 mutations, alternative MAPK 
pathway activation.

---

## Variant 4: EGFR p.T790M — Tier I

**Gene:** EGFR  
**Variant:** p.T790M (c.2369C>T), Missense Mutation  
**Clinical context:** Acquired resistance mutation  
**Hotspot:** Yes  

### Biological Mechanism
T790M ("gatekeeper" mutation) in EGFR exon 20 increases ATP binding affinity, 
sterically hindering 1st/2nd generation TKI binding. Accounts for ~60% of 
acquired resistance to erlotinib/gefitinib.

### Clinical Significance
**FDA-approved therapy:**
- Osimertinib (Tagrisso) — specifically designed to overcome T790M resistance
- Detectable in cfDNA at progression — liquid biopsy preferred over re-biopsy

**AMP/ASCO/CAP Classification: Tier I — Level A Evidence**  
T790M liquid biopsy testing at progression is NCCN-recommended standard of care.

### Liquid Biopsy — Key Clinical Application
This is where liquid biopsy shines: T790M can emerge in only a subset of 
metastatic lesions (spatial heterogeneity). A single tissue biopsy may miss 
T790M-positive clones that are detectable in cfDNA. Sensitivity of ddPCR-based 
cfDNA for T790M: ~70-80%.

---

## Variant 5: TP53 p.R175H — Tier II

**Gene:** TP53 (Tumor Protein P53)  
**Variant:** p.R175H (c.524G>A), Missense Mutation  
**Frequency in cohort:** TP53 mutated in 52% of LUAD patients  
**Hotspot:** Yes — one of the most common TP53 hotspot mutations  

### Biological Mechanism
R175H in the DNA-binding domain disrupts zinc coordination, abolishing 
sequence-specific DNA binding. This is a dominant-negative gain-of-function 
mutation — mutant p53 R175H not only loses tumor suppressor function but 
actively promotes oncogenesis by binding and inactivating wild-type p53, p63, 
and p73.

### Clinical Significance
No direct FDA-approved targeted therapy for TP53 mutations currently.  
**Clinical implications:**
- Associated with genomic instability and higher tumor mutational burden (TMB)
- High TMB (≥10 mut/Mb) → potential immunotherapy benefit (pembrolizumab)
- TP53 + STK11 co-mutation → poor immunotherapy response despite high TMB
- Emerging: APR-246 (eprenetapopt) targets mutant p53 — Phase III trials ongoing

**AMP/ASCO/CAP Classification: Tier II — Level B/C Evidence**  
Strong functional evidence but no direct targeted therapy approved.

---

## Summary Table — All Tier I Variants

| Gene | Variant | Patients | Hotspot | FDA Therapy | Tier |
|------|---------|----------|---------|-------------|------|
| KRAS | p.G12C | 151 (29.4%) | Yes | Sotorasib, Adagrasib | I |
| EGFR | p.L858R | 73 (14.2%) | Yes | Osimertinib | I |
| EGFR | p.T790M | detected | Yes | Osimertinib | I |
| EGFR | p.E746_A750del | detected | Yes | Osimertinib | I |
| BRAF | p.V600E | 42 (8.2%) | Yes | Dabrafenib+Trametinib | I |
| ERBB2 | p.G776delinsVC | 15 (2.9%) | Yes | Trastuzumab deruxtecan | I |
| MET | exon14 skip | 23 (4.5%) | — | Capmatinib, Tepotinib | I |

---

## Methodological Notes

### Variant Filtering Strategy
1. All variants → PASS filter (MUTECT/VARSCANS/MUSE consensus)
2. PASS variants → coding consequence filter
3. Coding variants → key cancer gene panel (17 genes)
4. Gene panel variants → hotspot + ClinVar + AMP/ASCO/CAP tiering

### Liquid Biopsy Relevance
VAF distribution analysis shows 1.7% of variants below 5% VAF — the typical 
detection threshold for liquid biopsy. This highlights the analytical challenge 
of cfDNA-based somatic variant detection and the importance of:
- High sequencing depth (>1000x for cfDNA)
- Molecular barcoding (UMIs) to suppress PCR errors
- Orthogonal confirmation for low-VAF variants

### Limitations
- ClinVar coverage of somatic variants is incomplete (germline-focused database)
- OncoKB therapeutic annotation pending API access
- BAM-level verification of low-VAF variants recommended (see IGV analysis, Day 4)

---

*Report generated as part of clinical genomics portfolio project.*  
*Data source: TCGA-LUAD via GDC (n=513 patients)*  
*Guidelines: AMP/ASCO/CAP 2017, ACMG/AMP 2015*
