"""
Step ②: mRNA-protein discordance by residual method
Uses PXD046631 matched RNA+protein abundance data
"""
import csv, os
from pathlib import Path
from statistics import mean, stdev
import numpy as np
from scipy import stats

BASE = Path(r"D:\Cardiac_Ogt_MultiOmics")
OUT = BASE / "output_tables"
OUT.mkdir(parents=True, exist_ok=True)

# ── 1. Load PXD046631 mRNA+protein merged data ─────────────────────
# This file has log2FC for both RNA and protein, plus mean abundance
print("=" * 60)
print("Step ②: mRNA-Protein Discordance (Residual Analysis)")
print("=" * 60)

data = []
with open(r"D:\miri_paper\results\proteomics\PXD046631_mRNA_protein_merged.csv") as f:
    reader = csv.DictReader(f)
    print(f"Columns: {reader.fieldnames}")
    for row in reader:
        try:
            gene = row["Gene"].strip()
            prot_log2fc = float(row["log2FC"])
            rna_log2fc = float(row["log2FoldChange"])
            mean_sham = float(row["mean_SHAM"])
            mean_ir = float(row["mean_MI_6h"])
            if mean_sham > 0 and mean_ir > 0:
                # protein abundance in log2 space
                prot_sham = np.log2(mean_sham)
                prot_ir = np.log2(mean_ir)
                data.append({
                    "gene": gene,
                    "prot_log2fc": prot_log2fc,
                    "rna_log2fc": rna_log2fc,
                    "prot_sham_log2": prot_sham,
                    "prot_ir_log2": prot_ir,
                    "mean_sham": mean_sham,
                    "mean_ir": mean_ir,
                })
        except (ValueError, KeyError, TypeError):
            pass

print(f"Genes with matched RNA+protein: {len(data)}")

if len(data) < 50:
    print("\n[WARNING] PXD046631 has limited matched genes. Trying alternative data...")
    # Check if there are larger proteomics datasets with matched RNA
    # PXD035154 proteinGroups.txt
    for alt in [
        r"D:\miri_paper\data\proteomics\PXD035154\fractionated_pool1_merged_results_20171214_1723_proteins.txt",
        r"D:\miri_paper\data\proteomics\PXD035154\fractionated_pool2_merged_results_20171217_1712_proteins.txt",
        r"D:\Cardiac_Ogt_MultiOmics\data_sources\02_原始数据\PXD035154\fractionated_pool1_merged_results_20171214_1723_proteins.txt",
    ]:
        if os.path.exists(alt):
            with open(alt) as fh:
                first = fh.readline().strip()
                print(f"  {os.path.basename(alt)}: header = {first[:100]}...")
            break

# ── 2. Residual analysis ──────────────────────────────────────────
if len(data) >= 50:
    rna_fc = [d["rna_log2fc"] for d in data]
    prot_fc = [d["prot_log2fc"] for d in data]

    # Linear regression: protein ~ RNA
    slope, intercept, r_val, p_val, std_err = stats.linregress(rna_fc, prot_fc)
    print(f"\nLinear regression: log2(protein_FC) = {intercept:.3f} + {slope:.3f} * log2(RNA_FC)")
    print(f"R = {r_val:.3f}, P = {p_val:.2e}")

    # Calculate residuals
    for d in data:
        predicted = intercept + slope * d["rna_log2fc"]
        d["residual"] = d["prot_log2fc"] - predicted

    # Studentize residuals
    residuals = [d["residual"] for d in data]
    res_mean = mean(residuals)
    res_std = stdev(residuals)
    for d in data:
        d["studentized_residual"] = (d["residual"] - res_mean) / res_std

    # Sort by |residual|
    data.sort(key=lambda d: abs(d["residual"]), reverse=True)

    # Top deviating genes
    print(f"\n--- Top 20 Genes with Largest |Residual| ---")
    print(f"{'Gene':12s} {'RNA_FC':>8s} {'Prot_FC':>8s} {'Residual':>8s}  {'Interpretation':>30s}")
    print("-" * 75)
    for d in data[:20]:
        if d["residual"] > 0:
            interp = "Protein >> RNA (transl. activation?)"
        else:
            interp = "Protein << RNA (degradation?)"
        print(f"{d['gene']:12s} {d['rna_log2fc']:+8.3f} {d['prot_log2fc']:+8.3f} {d['residual']:+8.3f}  {interp:30s}")

    # Concordance assessment
    up_both = sum(1 for d in data if d["rna_log2fc"] > 0 and d["prot_log2fc"] > 0)
    down_both = sum(1 for d in data if d["rna_log2fc"] < 0 and d["prot_log2fc"] < 0)
    discordant = sum(1 for d in data if d["rna_log2fc"] * d["prot_log2fc"] < 0)
    concordant = up_both + down_both
    total_with_dir = sum(1 for d in data if d["rna_log2fc"] != 0)

    print(f"\n--- Concordance Summary ---")
    print(f"Concordant (same direction): {concordant}/{len(data)} ({100*concordant/len(data):.1f}%)")
    print(f"  Both UP: {up_both}")
    print(f"  Both DOWN: {down_both}")
    print(f"Discordant (opposite direction): {discordant}/{len(data)} ({100*discordant/len(data):.1f}%)")

    # Genes of interest
    key_genes = ["OGT", "OGA", "NFE2L2", "HMOX1", "SLC7A11", "GPX4", "GCLC",
                 "RELA", "FTH1", "BACH1", "KEAP1", "TFRC", "NQO1"]
    print(f"\n--- Key Gene Residuals ---")
    print(f"{'Gene':10s} {'RNA_FC':>8s} {'Prot_FC':>8s} {'Residual':>8s} {'Status'}")
    print("-" * 50)
    for gene in key_genes:
        match = [d for d in data if d["gene"].upper() == gene.upper()]
        if match:
            d = match[0]
            status = "protein confirms RNA" if d["residual"]*d["rna_log2fc"] > 0 else "protein deviates"
            print(f"{gene:10s} {d['rna_log2fc']:+8.3f} {d['prot_log2fc']:+8.3f} {d['residual']:+8.3f} {status}")
        else:
            print(f"{gene:10s} {'not detected':>36s}")

    # Save results
    out_csv = OUT / "mRNA_protein_residuals.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene", "rna_log2fc", "prot_log2fc", "residual",
            "studentized_residual", "mean_sham", "mean_ir"
        ])
        writer.writeheader()
        for d in data:
            writer.writerow({k: d[k] for k in writer.fieldnames})
    print(f"\nSaved: {out_csv}")

else:
    print("\nInsufficient matched RNA-protein data for residual analysis.")
    print("Need: matched mRNA TPM + protein intensity at gene level across samples.")

print(f"\n=== DONE ===")
