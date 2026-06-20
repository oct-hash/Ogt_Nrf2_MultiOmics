import csv, json
from collections import Counter
from scipy.stats import fisher_exact

# Load OGT genes
with open(r"D:\Cardiac_Ogt_MultiOmics\validation_data\OGT_high_genes.json") as f:
    ogt_high = set(json.load(f))
print(f"OGT high-stringency interactors: {len(ogt_high)}")

# Load DEGs
deg_data = {}
with open(r"D:\miri_paper\results\tables\GSE193997_DESeq2_Sham_vs_IR6h.csv") as f:
    for row in csv.DictReader(f):
        try:
            if float(row["padj"]) < 0.05:
                g = row["symbol"].strip()
                deg_data[g] = {
                    "log2FC": float(row["log2FoldChange"]),
                    "padj": float(row["padj"]),
                }
        except (ValueError, KeyError):
            pass

deg_set = set(deg_data)
print(f"DEGs (padj < 0.05): {len(deg_set)}")

# Overlap
overlap = deg_set & ogt_high
print(f"\n{'='*50}")
print(f"OGT x DEG Overlap: {len(overlap)} genes")
print(f"{'='*50}")

# Fisher exact test
n_deg = len(deg_set)
n_ogt = len(ogt_high)
n_over = len(overlap)
n_bg = 20000

a = n_over
b = n_deg - a
c = n_ogt - a
d = n_bg - n_deg - c

or_val, p_val = fisher_exact([[a, b], [c, max(d, 1)]], alternative="greater")
print(f"\nFisher exact test (background = {n_bg} genes):")
print(f"  {a}/{n_deg} DEGs are OGT interactors ({100*a/n_deg:.2f}%)")
print(f"  {n_ogt}/{n_bg} background OGT interactors ({100*n_ogt/n_bg:.2f}%)")
print(f"  OR = {or_val:.2f}, P = {p_val:.8f}")
if p_val < 0.05:
    print(f"  >>> SIGNIFICANT enrichment")
else:
    power = 100 * a / n_over if n_over > 0 else 0
    print(f"  >>> NOT significant (P = {p_val:.4f})")

# Direction
up = sum(1 for g in overlap if deg_data[g]["log2FC"] > 0)
down = len(overlap) - up
print(f"\nOverlap genes direction: {up} UP, {down} DOWN")

# Key Nrf2/Ferroptosis genes
print(f"\n--- Key Nrf2/Ferroptosis Genes ---")
print(f"{'Gene':10s} {'In OGT-PIN':12s} {'DEG status':30s}")
print("-" * 55)
key = [
    "OGT", "OGA", "NFE2L2", "KEAP1", "BACH1",
    "HMOX1", "SLC7A11", "GPX4", "FTH1", "TFRC",
    "NQO1", "GCLC", "GCLM", "SRXN1", "SQSTM1",
    "RELA", "NFKB1", "STAT3", "HIF1A",
]
for g in key:
    in_ogt = "YES" if g in ogt_high else "-"
    if g in deg_data:
        d = deg_data[g]
        deg_str = f"DEG log2FC={d['log2FC']:+.2f} padj={d['padj']:.1e}"
    else:
        deg_str = "not differentially expressed"
    print(f"  {g:10s} {in_ogt:12s} {deg_str:30s}")

# Top 20 overlap genes by |log2FC|
print(f"\n--- Top 20 Overlap Genes (by |log2FC|) ---")
top_overlap = sorted(overlap, key=lambda g: abs(deg_data[g]["log2FC"]), reverse=True)[:20]
for g in top_overlap:
    d = deg_data[g]
    direction = "UP" if d["log2FC"] > 0 else "DOWN"
    print(f"  {g:12s} log2FC={d['log2FC']:+7.3f}  padj={d['padj']:.1e}  {direction}")

print(f"\n=== DONE ===")
