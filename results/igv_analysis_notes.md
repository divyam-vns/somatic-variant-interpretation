# BAM-Level Variant Verification Notes
## IGV Analysis — TCGA LUAD Key Variants

**Tool:** IGV v2.16 (igv.org/app), genome hg19/GRCh37  
**Purpose:** Read-level verification of somatic variant calls  
**Analyst:** Divya Mishra, Ph.D.

---

## What to Verify in IGV for Each Somatic Variant

For every Tier I/II variant entering a clinical report, BAM-level review
checks for the following quality indicators:

### 1. Read Depth (t_depth)
- Minimum 20x for somatic calling confidence
- Liquid biopsy requires >500x (cfDNA low tumor fraction)
- Our TCGA cohort median depth: 78x

### 2. Variant Allele Frequency (VAF)
- Calculated as: t_alt_count / t_depth
- Tumor VAF >10% = reliable somatic call in tissue
- Liquid biopsy: VAF can be <1% — requires UMI error correction

### 3. Strand Bias
- Alt reads should appear on BOTH forward and reverse strands
- >95% alt reads on one strand = likely sequencing artifact
- Flag for manual review if strand bias score >60

### 4. Base Quality
- Alt bases should have Phred quality >20
- Low base quality at variant position = potential sequencing error

### 5. Mapping Quality
- Reads supporting variant should have MAPQ >20
- Low MAPQ = reads mapping to multiple locations (repeat regions)

---

## Variant-Specific IGV Findings

### KRAS p.G12C — chr12:25,398,284
**Expected IGV appearance:**
- Reference base: G (green in IGV)
- Alt base: T (red in IGV) — G>T transversion causing G12C
- Expected VAF in tumor: 35-45% (tissue biopsy)
- In liquid biopsy: may be as low as 0.5-2%
- Strand distribution: should be balanced (~50/50 F/R)
- Artifact flags: none expected at this position
- **Clinical note:** KRAS codon 12 is a low-complexity region;
  confirm no homopolymer-induced errors

### EGFR p.L858R — chr7:55,259,515
**Expected IGV appearance:**
- Reference base: T (red in IGV)
- Alt base: G (yellow in IGV) — T>G transversion causing L858R
- Expected VAF in tumor: 40-60%
- Exon 21 — well-mappable region, high confidence calls
- **Clinical note:** Check for concurrent T790M (chr7:55,249,071)
  which may co-occur or emerge at resistance

### EGFR p.T790M — chr7:55,249,071
**Expected IGV appearance:**
- Reference base: C
- Alt base: T — C>T transition causing T790M
- VAF in resistance setting: often subclonal (5-30%)
- **Critical liquid biopsy note:** T790M may be present in only
  a subset of tumor clones — liquid biopsy more sensitive than
  single-site tissue biopsy for detecting this resistance mutation
- Spatial heterogeneity: different metastases may have different
  T790M status — cfDNA captures all circulating clones

### BRAF p.V600E — chr7:140,453,136
**Expected IGV appearance:**
- Reference base: A
- Alt base: T — A>T transversion causing V600E
- Expected VAF: 40-55% (typically clonal in LUAD)
- Well-mappable region in BRAF exon 15
- **Clinical note:** Verify no adjacent V600K/V600D variants
  which have different drug sensitivity profiles

### Low VAF Artifact Example — EGFR p.L858V (VAF=0.04)
**This is a teaching case for artifact recognition:**
- VAF of 4% is at the noise threshold for standard NGS
- IGV flags to check:
  * Strand bias: if all 4% alt reads are on forward strand = artifact
  * Base quality: if alt bases have Q<20 = likely error
  * Read position: if alt only at read ends = damaged DNA artifact
  * In liquid biopsy context: 4% VAF could be real — needs ddPCR confirmation
- **Conclusion:** Would flag for orthogonal confirmation before reporting

---

## Liquid Biopsy-Specific BAM Analysis Considerations

Unlike tissue BAM files, cfDNA liquid biopsy BAMs require:

1. **Ultra-deep sequencing** (1000-5000x) — low tumor fraction means
   true variants are rare among background noise
2. **UMI-based error correction** — molecular barcodes identify
   PCR duplicates and suppress sequencing errors below VAF 1%
3. **Fragment length analysis** — tumor-derived cfDNA is shorter
   (120-150bp) vs normal cfDNA (160-180bp); IGV fragment view
   can help distinguish tumor vs normal signal
4. **Strand concordance** — true variants appear on both strands
   of the original double-stranded DNA molecule

---

## Summary

| Variant | Position | Ref | Alt | Expected VAF | Artifact Risk |
|---------|----------|-----|-----|-------------|---------------|
| KRAS G12C | chr12:25,398,284 | G | T | 35-45% | Low |
| EGFR L858R | chr7:55,259,515 | T | G | 40-60% | Low |
| EGFR T790M | chr7:55,249,071 | C | T | 5-30% | Medium (low VAF) |
| BRAF V600E | chr7:140,453,136 | A | T | 40-55% | Low |
| EGFR L858V | chr7:55,259,515 | T | C | 4% | High — artifact |

*Note: VAF expectations based on TCGA tissue data.*
*Liquid biopsy VAFs would be 10-100x lower depending on tumor fraction.*
