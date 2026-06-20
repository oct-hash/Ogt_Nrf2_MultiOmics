"""Step ③: OGT enrichment with mouse-human mapping (case-insensitive)"""
import csv, json
from scipy.stats import fisher_exact

# Load OGT human genes
with open(r"D:\Cardiac_Ogt_MultiOmics\validation_data\OGT_high_genes.json") as f:
    ogt_high = set(json.load(f))
# Also make uppercase set for case-insensitive matching
ogt_upper = {g.upper() for g in ogt_high}

print(f"OGT high-stringency human interactors: {len(ogt_high)}")

# Load mouse DEGs
deg_data = {}
deg_mouse = {}
with open(r"D:\miri_paper\results\tables\GSE193997_DESeq2_Sham_vs_IR6h.csv") as f:
    for row in csv.DictReader(f):
        try:
            if float(row["padj"]) < 0.05:
                mouse_sym = row["symbol"].strip()
                human_guess = mouse_sym.upper()  # simplest ortholog mapping
                deg_data[mouse_sym] = {
                    "log2FC": float(row["log2FoldChange"]),
                    "padj": float(row["padj"]),
                    "human_guess": human_guess,
                }
                deg_mouse[mouse_sym] = human_guess
        except (ValueError, KeyError):
            pass

deg_set = set(deg_data)
deg_upper = {d["human_guess"] for d in deg_data.values()}
print(f"Mouse DEGs: {len(deg_set)}")
print(f"Unique human-guess symbols: {len(deg_upper)}")

# Intersection via case-insensitive matching
overlap_mouse = set()
overlap_genes = []
for mouse_g, d in deg_data.items():
    human_g = d["human_guess"]
    if human_g in ogt_upper:
        overlap_mouse.add(mouse_g)
        overlap_genes.append((human_g, mouse_g, d))

print(f"\n{'='*50}")
print(f"Overlap (case-insensitive): {len(overlap_mouse)} genes")
print(f"{'='*50}")

# Fisher
n_deg = len(deg_set)
n_ogt = len(ogt_high)
n_over = len(overlap_mouse)
n_bg = 20000

a = n_over; b = n_deg - a
c = n_ogt - a; d = n_bg - n_deg - c if n_bg - n_deg - c > 0 else 1

or_val, p_val = fisher_exact([[a, b], [c, d]], alternative="greater")
print(f"\nFisher exact test:")
print(f"  {a}/{n_deg} DEGs = OGT interactor orthologs ({100*a/n_deg:.2f}%)")
print(f"  {n_ogt}/{n_bg} background ({100*n_ogt/n_bg:.2f}%)")
print(f"  OR = {or_val:.3f}, P = {p_val:.6f}")
print(f"  {'>>> SIGNIFICANT' if p_val < 0.05 else '>>> NOT significant'}")

# Direction
up = sum(1 for _, _, d in overlap_genes if d["log2FC"] > 0)
down = len(overlap_genes) - up
print(f"\nDirection: {up} UP, {down} DOWN")

# Top overlap genes
print(f"\n--- Top 30 Overlap Genes (by |log2FC|) ---")
overlap_genes.sort(key=lambda x: abs(x[2]["log2FC"]), reverse=True)
for human_g, mouse_g, d in overlap_genes[:30]:
    direction = "UP" if d["log2FC"] > 0 else "DOWN"
    print(f"  {human_g:12s} ({mouse_g:12s}) log2FC={d['log2FC']:+7.3f}  {direction}")

# Key Nrf2/Ferroptosis genes
print(f"\n--- Key Nrf2/Ferroptosis Genes (case-insensitive match) ---")
key = {
    "OGT": "O-GlcNAc transferase",
    "OGA": "O-GlcNAcase (MGEA5)",
    "NFE2L2": "NRF2 master TF",
    "KEAP1": "NRF2 inhibitor",
    "BACH1": "NRF2 competitor",
    "HMOX1": "Ferroptosis / heme deg.",
    "SLC7A11": "Ferroptosis / cystine",
    "GPX4": "Ferroptosis / lipid ROS",
    "FTH1": "Iron storage",
    "TFRC": "Iron import",
    "NQO1": "NRF2 target",
    "GCLC": "Glutathione synthesis",
    "RELA": "NF-kB p65",
    "NFKB1": "NF-kB subunit",
    "STAT3": "STAT3 TF",
    "HIF1A": "HIF-1 alpha",
}
for human_g, desc in key.items():
    in_ogt = "YES" if human_g in ogt_high else "-"
    # Find mouse match
    mouse_match = None
    for mouse_g, d in deg_data.items():
        if d["human_guess"] == human_g:
            mouse_match = mouse_g
            fc = d["log2FC"]
            p = d["padj"]
            break
    if mouse_match:
        deg_str = f"DEG log2FC={fc:+.2f} padj={p:.1e}"
    else:
        deg_str = "-"
    print(f"  {human_g:10s} OGT: {in_ogt:4s}  {deg_str}  ({desc})")

print(f"\n=== DONE ===")
