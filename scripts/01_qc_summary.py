import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("Loading LUAD MAF file...")

cols = [
    "Hugo_Symbol", "Chromosome", "Start_Position",
    "Variant_Classification", "Variant_Type",
    "Tumor_Sample_Barcode", "HGVSp_Short",
    "t_depth", "t_ref_count", "t_alt_count",
    "FILTER", "COSMIC"
]

df = pd.read_csv(
    "/workspace/data/raw/luad_only.maf",
    sep="\t",
    usecols=cols,
    low_memory=False
)

print(f"Loaded {len(df):,} variants")

# Calculate VAF
df["t_depth"] = pd.to_numeric(df["t_depth"], errors="coerce")
df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
df["VAF"] = df["t_alt_count"] / df["t_depth"]

# PASS filter
df_pass = df[df["FILTER"] == "PASS"]

print(f"\n{'='*50}")
print("TCGA LUAD — QC SUMMARY")
print(f"{'='*50}")
print(f"Total variants       : {len(df):,}")
print(f"PASS variants        : {len(df_pass):,}")
print(f"Filtered out         : {len(df) - len(df_pass):,}")
print(f"Unique patients      : {df_pass['Tumor_Sample_Barcode'].nunique():,}")
print(f"Unique genes         : {df_pass['Hugo_Symbol'].nunique():,}")
print(f"Median VAF           : {df_pass['VAF'].median():.3f}")
print(f"Median read depth    : {df_pass['t_depth'].median():.0f}")
print(f"Variants per patient : {len(df_pass)/df_pass['Tumor_Sample_Barcode'].nunique():.0f}")

print(f"\n--- Top 20 Most Mutated Genes (PASS only) ---")
top_genes = df_pass["Hugo_Symbol"].value_counts().head(20)
for gene, count in top_genes.items():
    bar = "█" * (count // 20)
    print(f"  {gene:<12} {count:>5,}  {bar}")

print(f"\n--- Variant Classification Breakdown ---")
for vclass, count in df_pass["Variant_Classification"].value_counts().items():
    pct = count / len(df_pass) * 100
    print(f"  {vclass:<35} {count:>6,}  ({pct:.1f}%)")

print(f"\n--- VAF Distribution ---")
bins = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
labels = ["<5%", "5-10%", "10-20%", "20-30%", "30-50%", ">50%"]
df_pass = df_pass.copy()
    df_pass["vaf_bin"] = pd.cut(df_pass["VAF"], bins=bins, labels=labels)
for label, count in df_pass["vaf_bin"].value_counts().sort_index().items():
    pct = count / len(df_pass) * 100
    print(f"  VAF {label:<10} {count:>6,}  ({pct:.1f}%)")

# Check key cancer genes
print(f"\n--- Key Oncogenes & Tumor Suppressors ---")
key_genes = ["EGFR","KRAS","TP53","BRAF","PIK3CA","STK11",
             "KEAP1","MET","RET","ALK","ERBB2","RB1","CDKN2A"]
for gene in key_genes:
    subset = df_pass[df_pass["Hugo_Symbol"] == gene]
    n_patients = subset["Tumor_Sample_Barcode"].nunique()
    if n_patients > 0:
        pct = n_patients / df_pass["Tumor_Sample_Barcode"].nunique() * 100
        top_mut = subset["HGVSp_Short"].value_counts().index[0] if len(subset) > 0 else "N/A"
        print(f"  {gene:<10} {n_patients:>4} patients ({pct:.1f}%)  top mut: {top_mut}")

# Save PASS-only file
df_pass.to_csv("/workspace/data/raw/luad_pass.maf", sep="\t", index=False)
print(f"\nSaved PASS variants: /workspace/data/raw/luad_pass.maf")

# Plot 1 - top mutated genes
fig, ax = plt.subplots(figsize=(10, 7))
top_genes.plot(kind="barh", ax=ax, color="#185FA5")
ax.set_xlabel("Number of Mutations", fontsize=12)
ax.set_title("Top 20 Most Mutated Genes — TCGA LUAD (n=585)", fontsize=13)
ax.invert_yaxis()
plt.tight_layout()
plt.savefig("/workspace/results/top_mutated_genes.png", dpi=150)

# Plot 2 - VAF distribution
fig2, ax2 = plt.subplots(figsize=(10, 5))
df_pass["VAF"].dropna().hist(bins=50, ax=ax2, color="#185FA5", edgecolor="white")
ax2.set_xlabel("Variant Allele Frequency (VAF)", fontsize=12)
ax2.set_ylabel("Number of Variants", fontsize=12)
ax2.set_title("VAF Distribution — TCGA LUAD PASS Variants", fontsize=13)
ax2.axvline(x=0.05, color="red", linestyle="--", label="VAF=0.05 (liquid biopsy threshold)")
ax2.legend()
plt.tight_layout()
plt.savefig("/workspace/results/vaf_distribution.png", dpi=150)

print("Saved plots to /workspace/results/")
print(f"\n{'='*50}")
print("QC COMPLETE")
print(f"{'='*50}")
