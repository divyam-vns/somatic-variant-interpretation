import pandas as pd
import requests
import time
import re
from tqdm import tqdm

print("Loading LUAD PASS variants...")
df = pd.read_csv("/workspace/data/raw/luad_pass.maf", sep="\t", low_memory=False)

coding_types = [
    "Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
    "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation"
]
df_coding = df[df["Variant_Classification"].isin(coding_types)].copy()
print(f"Coding variants: {len(df_coding):,}")

key_genes = [
    "EGFR", "KRAS", "TP53", "BRAF", "PIK3CA", "STK11",
    "KEAP1", "MET", "RET", "ALK", "ERBB2", "RB1",
    "CDKN2A", "NF1", "SMAD4", "ARID1A", "U2AF1"
]
df_key = df_coding[df_coding["Hugo_Symbol"].isin(key_genes)].copy()
print(f"Key gene variants: {len(df_key):,}")

df_unique = df_key[["Hugo_Symbol","HGVSp_Short","Variant_Classification",
                     "Chromosome","Start_Position",
                     "Reference_Allele","Tumor_Seq_Allele2"]].drop_duplicates().copy()
print(f"Unique variants to annotate: {len(df_unique):,}")

# ── CANCER HOTSPOTS ──────────────────────────────────────────────────────────
print("\nQuerying CancerHotspots API...")

# Load all hotspots per gene into a cache first (one call per gene)
hotspot_gene_cache = {}

def get_gene_hotspots(gene):
    if gene not in hotspot_gene_cache:
        try:
            url = f"https://www.cancerhotspots.org/api/hotspots/single?hugoSymbol={gene}"
            r = requests.get(url, timeout=10)
            hotspot_gene_cache[gene] = r.json() if r.status_code == 200 else []
        except:
            hotspot_gene_cache[gene] = []
    return hotspot_gene_cache[gene]

def extract_position(hgvsp):
    """Extract numeric position from protein change e.g. p.L858R -> 858"""
    if not isinstance(hgvsp, str):
        return None
    match = re.search(r'(\d+)', hgvsp.replace("p.", ""))
    return int(match.group(1)) if match else None

def check_hotspot(gene, hgvsp):
    hotspots = get_gene_hotspots(gene)
    if not hotspots:
        return "Not_Hotspot", 0
    pos = extract_position(hgvsp)
    if pos is None:
        return "Not_Hotspot", 0
    for h in hotspots:
        aa_pos = h.get("aminoAcidPosition", {})
        if aa_pos.get("start") == pos:
            return "Hotspot", h.get("tumorCount", 0)
    return "Not_Hotspot", 0

hotspot_results = []
for _, row in tqdm(df_unique.iterrows(), total=len(df_unique), desc="CancerHotspots"):
    status, count = check_hotspot(row["Hugo_Symbol"], str(row["HGVSp_Short"]))
    hotspot_results.append({"hotspot_status": status, "hotspot_tumor_count": count})
    time.sleep(0.05)

df_unique = df_unique.reset_index(drop=True)
df_unique = pd.concat([df_unique, pd.DataFrame(hotspot_results)], axis=1)

# ── CLINVAR ──────────────────────────────────────────────────────────────────
print("\nQuerying ClinVar API...")

def query_clinvar(gene, hgvsp):
    try:
        # Extract short protein change e.g. p.L858R -> L858R
        aa = str(hgvsp).replace("p.", "").strip()
        if not aa or aa == "nan":
            return "Not_found", ""

        url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=clinvar&term={gene}[gene]+AND+{aa}[variant+name]"
            f"&retmax=1&retmode=json"
        )
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return "Not_found", ""

        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return "Not_found", ""

        # Fetch clinical significance
        url2 = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=clinvar&id={ids[0]}&retmode=json"
        )
        r2 = requests.get(url2, timeout=10)
        result = r2.json().get("result", {}).get(ids[0], {})
        clin_sig = result.get("clinical_significance", {}).get("description", "Not_found")
        traits = result.get("trait_set", [])
        condition = traits[0].get("trait_name", "") if traits else ""
        return clin_sig, condition
    except:
        return "Not_found", ""

# Annotate all unique variants with ClinVar
clinvar_sig = []
clinvar_condition = []

for _, row in tqdm(df_unique.iterrows(), total=len(df_unique), desc="ClinVar"):
    sig, cond = query_clinvar(row["Hugo_Symbol"], row["HGVSp_Short"])
    clinvar_sig.append(sig)
    clinvar_condition.append(cond)
    time.sleep(0.34)

df_unique["clinvar_significance"] = clinvar_sig
df_unique["clinvar_condition"] = clinvar_condition

# ── AMP/ASCO/CAP TIER ASSIGNMENT ─────────────────────────────────────────────
print("\nAssigning AMP/ASCO/CAP tiers...")

# Known Tier I variants - strong clinical significance with FDA approved therapy
tier1_variants = {
    ("EGFR", "p.L858R"), ("EGFR", "p.T790M"), ("EGFR", "p.E746_A750del"),
    ("KRAS", "p.G12C"), ("BRAF", "p.V600E"),
    ("ERBB2", "p.G776delinsVC"), ("MET", "p.X1028_splice"),
    ("RET", "p.M918T"), ("ALK", "p.R1275Q"),
}

# Known Tier II variants - potential clinical significance
tier2_genes = {"PIK3CA", "STK11", "KEAP1", "NF1", "RB1", "CDKN2A"}

def assign_tier(row):
    key = (row["Hugo_Symbol"], str(row["HGVSp_Short"]))
    clin_sig = str(row["clinvar_significance"]).lower()
    hotspot = row["hotspot_status"] == "Hotspot"

    if key in tier1_variants:
        return "Tier_I"
    if "pathogenic" in clin_sig and hotspot:
        return "Tier_I"
    if hotspot or "pathogenic" in clin_sig or row["Hugo_Symbol"] in tier2_genes:
        return "Tier_II"
    if "uncertain" in clin_sig or "conflicting" in clin_sig:
        return "Tier_III"
    return "Tier_IV"

df_unique["amp_tier"] = df_unique.apply(assign_tier, axis=1)

# ── MERGE AND SAVE ────────────────────────────────────────────────────────────
print("\nMerging annotations...")
df_annotated = df_key.merge(
    df_unique[["Hugo_Symbol","HGVSp_Short","hotspot_status",
               "hotspot_tumor_count","clinvar_significance",
               "clinvar_condition","amp_tier"]],
    on=["Hugo_Symbol","HGVSp_Short"],
    how="left"
)

df_annotated.to_csv("/workspace/data/annotated/luad_annotated.tsv", sep="\t", index=False)
df_unique.to_csv("/workspace/data/annotated/unique_variants_annotated.tsv", sep="\t", index=False)

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print(f"\n{'='*55}")
print("ANNOTATION SUMMARY — TCGA LUAD KEY GENES")
print(f"{'='*55}")
print(f"Total coding variants in key genes : {len(df_key):,}")
print(f"Unique variants annotated          : {len(df_unique):,}")

print(f"\n--- Hotspot Status ---")
for status, count in df_unique["hotspot_status"].value_counts().items():
    print(f"  {status:<20} {count:>4}")

print(f"\n--- Top Hotspot Variants ---")
hotspots = df_unique[df_unique["hotspot_status"]=="Hotspot"].sort_values(
    "hotspot_tumor_count", ascending=False).head(15)
for _, row in hotspots.iterrows():
    print(f"  {row['Hugo_Symbol']:<8} {str(row['HGVSp_Short']):<15} "
          f"tumors: {row['hotspot_tumor_count']:>5,}  "
          f"tier: {row['amp_tier']}")

print(f"\n--- AMP/ASCO/CAP Tier Distribution ---")
for tier, count in df_unique["amp_tier"].value_counts().sort_index().items():
    pct = count / len(df_unique) * 100
    print(f"  {tier:<12} {count:>4}  ({pct:.1f}%)")

print(f"\n--- ClinVar Findings ---")
for sig, count in df_unique["clinvar_significance"].value_counts().head(8).items():
    print(f"  {sig:<40} {count:>4}")

print(f"\n--- Tier I & II Variants (Clinically Actionable) ---")
actionable = df_unique[df_unique["amp_tier"].isin(["Tier_I","Tier_II"])].sort_values("amp_tier")
for _, row in actionable.iterrows():
    print(f"  {row['amp_tier']}  {row['Hugo_Symbol']:<8} "
          f"{str(row['HGVSp_Short']):<18} "
          f"hotspot: {row['hotspot_status']:<12} "
          f"clinvar: {row['clinvar_significance']}")

print(f"\n{'='*55}")
print("ANNOTATION COMPLETE")
print(f"Saved: /workspace/data/annotated/luad_annotated.tsv")
print(f"Saved: /workspace/data/annotated/unique_variants_annotated.tsv")
print(f"{'='*55}")
