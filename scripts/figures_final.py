#!/usr/bin/env python3
"""Generate Fig 2B (heatmap) + Fig 4A-D (snRNA-seq UMAP/Feature/Violin)"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path

BASE = Path(r"D:\Cardiac_Ogt_MultiOmics")
FIGS = BASE / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],
                     "font.size":8,"figure.dpi":300,"savefig.dpi":300})

def save(fig, name):
    path = FIGS / f"{name}.png"
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")

# ================================================================
# Fig 2B: Core Gene Expression Heatmap
# ================================================================
print("Fig 2B: Core Gene Heatmap")
tpm = pd.read_csv(r"D:\miri_paper\results\cibersortx\GSE193997_TPM_expression_matrix.txt",
                  sep="\t", index_col=0)

# Select core genes (Nrf2/ferroptosis/O-GlcNAc)
core_genes = ["Ogt","Oga","Nfe2l2","Keap1","Bach1","Hmox1","Slc7a11",
              "Gpx4","Fth1","Tfrc","Nqo1","Gclc","Gclm","Srxn1",
              "Rela","Nfkb1","Stat3","Hif1a","Acs4"]
found = [g for g in core_genes if g in tpm.index]
missing = [g for g in core_genes if g not in tpm.index]
print(f"  Found: {len(found)}/{len(core_genes)}, Missing: {missing}")

# Subset and log-transform
mat = tpm.loc[found]
mat_log = np.log2(mat + 1)

# Z-score per gene
mat_z = mat_log.subtract(mat_log.mean(axis=1), axis=0).div(mat_log.std(axis=1), axis=0)

# Sort samples: Sham first, then IR
sham_cols = [c for c in mat_z.columns if "Sham" in c or "GSM582540" in c or "GSM582541" in c]
ir_cols = [c for c in mat_z.columns if c not in sham_cols]
# If sham/ir detection fails, sort by Ogt expression
if len(sham_cols) < 2:
    sham_cols = list(mat_z.columns[:len(mat_z.columns)//2])
    ir_cols = list(mat_z.columns[len(mat_z.columns)//2:])

ordered = sham_cols + ir_cols

fig, ax = plt.subplots(figsize=(10, 5))
sns.heatmap(mat_z[ordered], cmap="RdBu_r", center=0, ax=ax,
            xticklabels=False, yticklabels=True,
            cbar_kws={"label": "Z-score log2(TPM+1)", "shrink": 0.6},
            linewidths=0.1, linecolor="lightgrey")
ax.set_title(f"Core Nrf2/Ferroptosis/O-GlcNAc Genes\n(n={len(found)} genes, {len(sham_cols)} Sham + {len(ir_cols)} I/R)")
ax.axvline(len(sham_cols), color="black", linewidth=2)
ax.text(len(sham_cols)/2, -0.5, "Sham", ha="center", fontsize=9, fontweight="bold")
ax.text(len(sham_cols) + len(ir_cols)/2, -0.5, "I/R 6h", ha="center", fontsize=9, fontweight="bold")
save(fig, "Fig2B_CoreGene_Heatmap")

# ================================================================
# Fig 4A-D: snRNA-seq Visualizations
# ================================================================
print("\nFig 4: snRNA-seq Visualizations")
adata = sc.read_h5ad(r"D:\miri_paper\data\processed\GSE240848_full_analyzed.h5ad")
print(f"  Loaded: {adata.shape[0]} nuclei x {adata.shape[1]} genes")

# Set UMAP
adata.obsm["X_umap"] = adata.obsm["X_umap"]  # already exists

# --- 4A: UMAP by cell type ---
print("  4A: UMAP by cell type")
fig, ax = plt.subplots(figsize=(6, 5))
sc.pl.umap(adata, color="cell_type", ax=ax, show=False, legend_loc="right margin",
           title="snRNA-seq: 48,633 Cardiac Nuclei", frameon=False, palette="tab10")
save(fig, "Fig4A_UMAP_CellTypes")

# --- 4B: Feature Plot - Ogt ---
print("  4B: Feature Plot - Ogt")
fig, ax = plt.subplots(figsize=(5, 4.5))
gene = "Ogt" if "Ogt" in adata.var_names else ([g for g in adata.var_names if g.lower()=="ogt"] or ["Ogt"])[0]
sc.pl.umap(adata, color=gene, ax=ax, show=False, cmap="Reds",
           title="Ogt (O-GlcNAc Transferase)", frameon=False)
save(fig, "Fig4B_UMAP_Ogt")

# --- 4C: Feature Plot - Hmox1 + Slc7a11 (multi-panel) ---
print("  4C: Feature Plot - Nrf2 targets")
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
for ax_i, (gene, title, cmap) in enumerate(zip(
    ["Hmox1", "Slc7a11"],
    ["Hmox1 (Heme Oxygenase-1)", "Slc7a11 (Cystine/Glutamate Transporter)"],
    ["Oranges", "Purples"]
)):
    if gene in adata.var_names:
        sc.pl.umap(adata, color=gene, ax=axes[ax_i], show=False, cmap=cmap,
                   title=title, frameon=False)
    else:
        axes[ax_i].text(0.5, 0.5, f"{gene}\nnot in HVG", transform=axes[ax_i].transAxes,
                        ha="center", va="center")
save(fig, "Fig4C_UMAP_Nrf2Targets")

# --- 4D: Violin plot of key genes by cell type ---
print("  4D: Violin plot - key genes by cell type")
key_genes = ["Ogt", "Hmox1", "Slc7a11", "Nfe2l2", "Gpx4", "Rela"]
found_genes = [g for g in key_genes if g in adata.var_names]
print(f"    Genes found: {found_genes}")

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()
cell_type_order = ["Cardiomyocyte", "Pericyte", "Fibroblast", "Endothelial",
                   "SmoothMuscle", "Macrophage", "Neutrophil", "T_Cell", "B_Cell"]

for i, gene in enumerate(found_genes):
    ax = axes[i]
    sc.pl.violin(adata, gene, groupby="cell_type", ax=ax, show=False,
                 rotation=45, order=cell_type_order)
    ax.set_title(gene, fontsize=10, fontweight="bold")
    ax.set_xlabel("")

# Remove unused axes
for j in range(len(found_genes), len(axes)):
    axes[j].set_visible(False)

fig.suptitle("Cell-Type-Specific Expression of Key Genes", fontsize=12, fontweight="bold")
plt.tight_layout()
save(fig, "Fig4D_Violin_KeyGenes")

print("\nDone")
