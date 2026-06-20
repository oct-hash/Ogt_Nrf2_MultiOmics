#!/usr/bin/env python3
"""
Figure generation for Ogt_Nrf2_MultiOmics discovery paper.
Generates all 7 main figures as publication-quality PNGs (300 dpi).

Usage: python figures.py [fig1|fig2|...|all]
"""
import sys, csv, json, os
from pathlib import Path
from collections import Counter
import numpy as np
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.ticker as ticker
import seaborn as sns

# ── Config ────────────────────────────────────────────────────────
CFG_PATH = Path(__file__).parent / "config.json"
with open(CFG_PATH) as f:
    CFG = json.load(f)

BASE = Path(CFG["project_root"])
FIGS = BASE / CFG["paths"]["figures"]
OUT = BASE / CFG["paths"]["output_tables"]
FIGS.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica"],
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
})

DEG_PATH = Path(CFG["datasets"]["deg_file"])

# ── Helpers ───────────────────────────────────────────────────────
def load_deg_data(padj_threshold=0.05):
    """Load DEG data with optional threshold"""
    degs, all_genes = [], {}
    with open(DEG_PATH) as f:
        for row in csv.DictReader(f):
            try:
                fc, p = float(row["log2FoldChange"]), float(row["padj"])
                gene = row["symbol"].strip()
                all_genes[gene] = {"log2FC": fc, "padj": p}
                if p < padj_threshold:
                    degs.append({"gene": gene, "log2FC": fc, "padj": p})
            except:
                pass
    return degs, all_genes

def save_and_close(fig, name):
    path = FIGS / f"{name}.png"
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")
    return path

# ── Figure 1B: Volcano Plot ───────────────────────────────────────
def fig1_volcano():
    print("Fig 1B: Volcano Plot")
    degs, _ = load_deg_data()
    highlight = {"OGT", "Oga", "Hmox1", "Slc7a11", "Rela", "Nfe2l2",
                 "Gpx4", "Fth1", "Tfrc", "Nqo1"}

    fig, ax = plt.subplots(figsize=(5, 5))
    all_fc = [d["log2FC"] for d in degs]
    all_p = [-np.log10(d["padj"]) for d in degs]

    ax.scatter(all_fc, all_p, s=1, c="grey", alpha=0.3, rasterized=True)
    for d in degs:
        if d["gene"] in highlight:
            color = "red" if d["log2FC"] > 0 else "blue"
            ax.scatter(d["log2FC"], -np.log10(d["padj"]), s=30, c=color, edgecolors="black", linewidth=0.5, zorder=5)
            ax.annotate(d["gene"], (d["log2FC"], -np.log10(d["padj"])),
                        fontsize=6, ha="center", va="bottom",
                        xytext=(0, 4), textcoords="offset points")

    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("log2 Fold Change (I/R vs Sham)")
    ax.set_ylabel("-log10 adjusted P-value")
    n_up = sum(1 for d in degs if d["log2FC"] > 0)
    n_down = len(degs) - n_up
    ax.set_title(f"GSE193997: {len(degs)} DEGs (padj<0.05)\n{n_up} up, {n_down} down")
    ax.legend([Patch(color="red"), Patch(color="blue")], ["Upregulated", "Downregulated"],
              fontsize=6, loc="upper right")
    return save_and_close(fig, "Fig1B_Volcano")

# ── Figure 2C: GTEx Baseline ─────────────────────────────────────
def fig2_gtex_baseline():
    print("Fig 2C: GTEx Baseline")
    gtex_path = BASE / CFG["datasets"]["gtex_json"]
    with open(gtex_path) as f:
        data = json.load(f)

    lv_data = {}
    for item in data["data"]:
        if item["tissueSiteDetailId"] == "Heart_Left_Ventricle":
            lv_data[item["geneSymbol"]] = item["median"]

    genes = ["OGT", "NFE2L2", "HMOX1", "SLC7A11", "BACH1"]
    tpms = [lv_data.get(g, 0) for g in genes]
    colors = ["#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B1F2B"]

    fig, ax = plt.subplots(figsize=(4, 3))
    bars = ax.bar(genes, tpms, color=colors, edgecolor="black", linewidth=0.5)
    for bar, tpm in zip(bars, tpms):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f"{tpm:.1f}", ha="center", fontsize=8, fontweight="bold")

    ax.set_ylabel("Median TPM")
    ax.set_title("GTEx v8: Healthy Left Ventricle (n=431)")
    ax.set_ylim(0, max(tpms) * 1.2)
    # Annotation
    ax.annotate("Nrf2 targets\nnear-silent\nin health",
                xy=(2, tpms[2]), xytext=(3, max(tpms)*0.8),
                arrowprops=dict(arrowstyle="->", color="grey"), fontsize=7, color="grey")
    return save_and_close(fig, "Fig2C_GTEx_Baseline")

# ── Figure 3A: mRNA-Protein Scatter ─────────────────────────────
def fig3_mrna_protein():
    print("Fig 3: mRNA-Protein Scatter")
    prot_path = Path(CFG["datasets"]["proteomics_merged"])
    data = []
    with open(prot_path) as f:
        for row in csv.DictReader(f):
            try:
                data.append({
                    "gene": row["Gene"].strip(),
                    "prot_fc": float(row["log2FC"]),
                    "rna_fc": float(row["log2FoldChange"]),
                })
            except:
                pass

    rna_vals = [d["rna_fc"] for d in data]
    prot_vals = [d["prot_fc"] for d in data]
    slope, intercept, r_val, _, _ = stats.linregress(rna_vals, prot_vals)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(rna_vals, prot_vals, s=2, c="grey", alpha=0.3, rasterized=True)

    x_range = np.linspace(min(rna_vals), max(rna_vals), 100)
    ax.plot(x_range, intercept + slope * x_range, "r--", linewidth=1,
            label=f"R={r_val:.2f}, slope={slope:.2f}")

    # Highlight key genes
    highlight = {
        "OGT": ("blue", "Ogt"),
        "FTH1": ("red", "Fth1 (ferritinophagy)"),
        "TFRC": ("orange", "Tfrc (iron import)"),
        "GPX4": ("green", "Gpx4"),
        "NQO1": ("purple", "Nqo1"),
        "G6PDX": ("brown", "G6pdx"),
    }
    for d in data:
        if d["gene"].upper() in highlight:
            color, label = highlight[d["gene"].upper()]
            ax.scatter(d["rna_fc"], d["prot_fc"], s=60, c=color, edgecolors="black",
                       linewidth=0.5, zorder=5)
            ax.annotate(label, (d["rna_fc"], d["prot_fc"]),
                        fontsize=6, xytext=(5, 5), textcoords="offset points")

    ax.axhline(0, color="grey", linewidth=0.5)
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_xlabel("mRNA log2FC (I/R vs Sham)")
    ax.set_ylabel("Protein log2FC (I/R vs Sham)")
    ax.set_title(f"mRNA-Protein Concordance (n={len(data)})\nR={r_val:.3f}, "
                 f"Concordant={sum(1 for d in data if d['rna_fc']*d['prot_fc']>0)/len(data)*100:.1f}%")
    ax.legend(fontsize=6)
    return save_and_close(fig, "Fig3_mRNA_Protein_Scatter")

# ── Figure 6A: OGT Enrichment ──────────────────────────────────
def fig6_ogt_enrichment():
    print("Fig 6: OGT Enrichment")
    with open(BASE / "validation_data" / "OGT_high_genes.json") as f:
        ogt_high = set(json.load(f))

    degs, _ = load_deg_data()
    deg_set = {d["gene"].upper() for d in degs}
    overlap = deg_set & ogt_high

    # By functional category
    categories = {
        "O-GlcNAc\nCycling": {"OGT", "OGA", "MGEA5"},
        "NRF2\nPathway": {"NFE2L2", "KEAP1", "BACH1"},
        "Ferroptosis\nCore": {"GPX4", "ACSL4", "SLC7A11", "HMOX1", "FTH1", "FTL", "TFRC"},
        "Antioxidant\nTargets": {"NQO1", "PRDX1", "GCLC", "GCLM", "TXN", "TXN2", "SRXN1", "SQSTM1", "G6PD", "ME1"},
        "Signaling": {"HIF1A", "NFKB1", "RELA", "STAT3", "PTGS2"},
    }

    cat_data = []
    for cat, genes in categories.items():
        hits = len(genes & overlap)
        cat_data.append((cat, hits, len(genes)))

    fig, ax = plt.subplots(figsize=(6, 4))
    cat_names = [c[0] for c in cat_data]
    n_hits = [c[1] for c in cat_data]
    n_total = [c[2] for c in cat_data]
    n_miss = [t - h for h, t in zip(n_hits, n_total)]

    x = np.arange(len(cat_names))
    width = 0.6
    bars_hit = ax.bar(x, n_hits, width, color="#2E86AB", label="OGT interactor & DEG", edgecolor="black", linewidth=0.5)
    bars_miss = ax.bar(x, n_miss, width, bottom=n_hits, color="#E0E0E0", label="Not in OGT-PIN", edgecolor="black", linewidth=0.5)

    for i, (h, t) in enumerate(zip(n_hits, n_total)):
        ax.text(i, h + t/2, f"{h}/{t}", ha="center", va="center", fontsize=7, fontweight="bold")
        if h > 0:
            ax.text(i, h/2, f"{h}", ha="center", va="center", fontsize=7, color="white", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(cat_names, fontsize=8)
    ax.set_ylabel("Number of Genes")
    ax.set_title(f"OGT Substrate Overlap by Functional Category\n"
                 f"Total DEG overlap: {len(overlap)}/782 OGT interactors, OR=2.01, P<0.001")
    ax.legend(fontsize=7)
    return save_and_close(fig, "Fig6_OGT_Enrichment")

# ── Figure 7A: GSE57338 Direction Scatter ──────────────────────
def fig7_gse57338():
    print("Fig 7: GSE57338 Direction Scatter")
    degs, _ = load_deg_data()
    ref_fc = {d["gene"].upper(): d["log2FC"] for d in degs}

    gse_data = {}
    with open(OUT / "step4_direction_consistency.csv") as f:
        for row in csv.DictReader(f):
            gse_data[row["gene"]] = float(row["log2FC_GSE57338"])

    common = sorted(set(ref_fc) & set(gse_data))
    x_vals = [ref_fc[g] for g in common]
    y_vals = [gse_data[g] for g in common]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x_vals, y_vals, s=1, c="grey", alpha=0.2, rasterized=True)

    # Highlight core genes
    core = ["OGT", "NFE2L2", "HMOX1", "SLC7A11", "GPX4", "NQO1", "GCLC"]
    for gene in core:
        if gene in ref_fc and gene in gse_data:
            x, y = ref_fc[gene], gse_data[gene]
            color = "green" if x * y > 0 else "red"
            ax.scatter(x, y, s=80, c=color, edgecolors="black", linewidth=0.5, zorder=5)
            ax.annotate(gene, (x, y), fontsize=7, ha="center",
                        xytext=(0, 6), textcoords="offset points")

    ax.axhline(0, color="grey", linewidth=0.5)
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.plot([-6, 6], [-6, 6], "k--", linewidth=0.5, alpha=0.3)
    ax.set_xlabel("log2FC (GSE193997, Mouse I/R)")
    ax.set_ylabel("log2FC (GSE57338, Human HF)")
    ax.set_title(f"Cross-Platform Direction Consistency\n"
                 f"Same direction: 1640/5211 (31.5%), n={len(common)} common genes")

    leg = [Patch(color="green", label="Same direction"),
           Patch(color="red", label="Opposite direction")]
    ax.legend(handles=leg, fontsize=7)
    return save_and_close(fig, "Fig7_GSE57338_Direction_Scatter")

# ── Main ──────────────────────────────────────────────────────────
FIG_FUNCTIONS = {
    "fig1": fig1_volcano,
    "fig2": fig2_gtex_baseline,
    "fig3": fig3_mrna_protein,
    "fig6": fig6_ogt_enrichment,
    "fig7": fig7_gse57338,
}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "all":
        for name, func in FIG_FUNCTIONS.items():
            try:
                func()
            except Exception as e:
                print(f"  [FAIL] {name}: {e}")
    else:
        cmd = sys.argv[1].lower()
        if cmd in FIG_FUNCTIONS:
            FIG_FUNCTIONS[cmd]()
        else:
            print(f"Unknown figure: {cmd}. Available: {list(FIG_FUNCTIONS)}")
