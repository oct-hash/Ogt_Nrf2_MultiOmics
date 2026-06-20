#!/usr/bin/env python3
"""
Assemble individual panels into multi-panel Cell Reports figures.
Adds A/B/C/D labels and consistent styling.

Usage: python assemble_multipanel.py
Output: D:\Cardiac_Ogt_MultiOmics\figures\Fig{1-7}_Final.png
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle

BASE = Path(r"D:\Cardiac_Ogt_MultiOmics")
FIGS = BASE / "figures"
OUT = BASE / "figures"  # output to same directory

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica"],
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

def add_panel_label(ax, label, x=-0.02, y=1.02, fontsize=14):
    """Add bold panel label in top-left corner"""
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
            fontweight="bold", ha="left", va="bottom")

def save_final(fig, name):
    path = FIGS / name
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    print(f"  Saved: {path}")
    plt.close(fig)

# ================================================================
# Figure 1: Study Overview + Volcano
# ================================================================
def assemble_fig1():
    print("Assembling Figure 1...")
    img_a = mpimg.imread(str(FIGS / "Fig1A_StudyOverview.png"))
    img_b = mpimg.imread(str(FIGS / "Fig1B_Volcano.png"))

    fig = plt.figure(figsize=(16, 7))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.imshow(img_a)
    ax1.axis("off")
    add_panel_label(ax1, "A")

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(img_b)
    ax2.axis("off")
    add_panel_label(ax2, "B")

    plt.tight_layout(pad=0.5)
    save_final(fig, "Fig1_Final")

# ================================================================
# Figure 2: GSEA + Heatmap + GTEx Baseline
# ================================================================
def assemble_fig2():
    print("Assembling Figure 2...")
    img_a = mpimg.imread(str(FIGS / "Fig2A_GSEA_Nrf2ARE.png"))
    img_b = mpimg.imread(str(FIGS / "Fig2B_CoreGene_Heatmap.png"))
    img_c = mpimg.imread(str(FIGS / "Fig2C_GTEx_Baseline.png"))

    fig = plt.figure(figsize=(16, 12))
    # A: GSEA (top-left, wider)
    ax1 = fig.add_subplot(2, 3, (1, 3))  # spans top row, full width
    ax1.imshow(img_a)
    ax1.axis("off")
    add_panel_label(ax1, "A")

    # B: Heatmap (bottom-left, wider)
    ax2 = fig.add_subplot(2, 3, (4, 5))
    ax2.imshow(img_b)
    ax2.axis("off")
    add_panel_label(ax2, "B")

    # C: GTEx baseline (bottom-right)
    ax3 = fig.add_subplot(2, 3, 6)
    ax3.imshow(img_c)
    ax3.axis("off")
    add_panel_label(ax3, "C")

    plt.tight_layout(pad=0.3)
    save_final(fig, "Fig2_Final")

# ================================================================
# Figure 3: mRNA-Protein — already multi-panel, just copy with label
# ================================================================
def assemble_fig3():
    print("Assembling Figure 3 (already multi-panel)...")
    img = mpimg.imread(str(FIGS / "Fig3_mRNA_Protein_Final.png"))

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.imshow(img)
    ax.axis("off")
    # The existing figure already has sub-panel titles, add overall label
    ax.text(0, 1.01, "Figure 3", transform=ax.transAxes, fontsize=14,
            fontweight="bold", ha="left", va="bottom")

    save_final(fig, "Fig3_Final")

# ================================================================
# Figure 4: snRNA-seq UMAP + Feature + Violin — 4 panels
# ================================================================
def assemble_fig4():
    print("Assembling Figure 4...")
    img_a = mpimg.imread(str(FIGS / "Fig4A_UMAP_CellTypes.png"))
    img_b = mpimg.imread(str(FIGS / "Fig4B_UMAP_Ogt.png"))
    img_c = mpimg.imread(str(FIGS / "Fig4C_UMAP_Nrf2Targets.png"))
    img_d = mpimg.imread(str(FIGS / "Fig4D_Violin_KeyGenes.png"))

    fig = plt.figure(figsize=(16, 14))
    # A: UMAP cell types (top-left)
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.imshow(img_a)
    ax1.axis("off")
    add_panel_label(ax1, "A")

    # B: Ogt feature (top-middle)
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.imshow(img_b)
    ax2.axis("off")
    add_panel_label(ax2, "B")

    # C: Nrf2 targets (top-right)
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.imshow(img_c)
    ax3.axis("off")
    add_panel_label(ax3, "C")

    # D: Violin (bottom, full width)
    ax4 = fig.add_subplot(2, 3, (4, 6))
    ax4.imshow(img_d)
    ax4.axis("off")
    add_panel_label(ax4, "D")

    plt.tight_layout(pad=0.3)
    save_final(fig, "Fig4_Final")

# ================================================================
# Figure 5: CCC — already multi-panel
# ================================================================
def assemble_fig5():
    print("Assembling Figure 5 (already multi-panel)...")
    img = mpimg.imread(str(FIGS / "Fig5_CCC_Final.png"))
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.imshow(img)
    ax.axis("off")
    save_final(fig, "Fig5_Final")

# ================================================================
# Figure 6: OGT-PIN — already multi-panel
# ================================================================
def assemble_fig6():
    print("Assembling Figure 6 (already multi-panel)...")
    img = mpimg.imread(str(FIGS / "Fig6_OGT_RELA_Final.png"))
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.imshow(img)
    ax.axis("off")
    save_final(fig, "Fig6_Final")

# ================================================================
# Figure 7: Validation + Graphical Abstract — already multi-panel
# ================================================================
def assemble_fig7():
    print("Assembling Figure 7 (already multi-panel)...")
    img = mpimg.imread(str(FIGS / "Fig7_Validation_GraphicalAbstract.png"))
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.imshow(img)
    ax.axis("off")
    save_final(fig, "Fig7_Final")


if __name__ == "__main__":
    print("Assembling multi-panel figures...")
    print("="*60)
    assemble_fig1()
    assemble_fig2()
    assemble_fig3()
    assemble_fig4()
    assemble_fig5()
    assemble_fig6()
    assemble_fig7()
    print("="*60)
    print("ALL 7 MULTI-PANEL FIGURES ASSEMBLED")
    print(f"Output: {FIGS}/Fig{{1-7}}_Final.png")
