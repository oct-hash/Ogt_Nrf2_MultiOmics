"""
Step ④ Part B: GSE57338 HF vs NF DEG + direction consistency
Corrected: parse heart failure yes/no from sample characteristics
"""
import csv, json, gzip
from pathlib import Path
from statistics import mean
import numpy as np
from scipy import stats

BASE = Path(r"D:\Cardiac_Ogt_MultiOmics")
OUT = BASE / "output_tables"
OUT.mkdir(parents=True, exist_ok=True)

gse_path = BASE / "validation_data" / "GSE57338_series_matrix.txt.gz"

print("Parsing GSE57338 series matrix...")

# Parse phenotype from sample characteristics
hf_mask = []
samples = []
with gzip.open(gse_path, "rt", errors="ignore") as f:
    lines = f.readlines()

i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.startswith("!Sample_title"):
        samples = [s.strip('"') for s in line.split("\t")[1:]]
    elif line.startswith("!Sample_characteristics_ch1") and "heart failure" in line:
        parts = line.split("\t")
        hf_mask = ["heart failure: yes" in p for p in parts[1:]]
    i += 1

hf_idx = [j for j, v in enumerate(hf_mask) if v]
nf_idx = [j for j, v in enumerate(hf_mask) if not v]
print(f"Samples: {len(samples)}, HF: {len(hf_idx)}, NF: {len(nf_idx)}")

# Parse expression matrix and probe-gene mapping
probe_to_gene = {}
expr_data = {}

i = 0
while i < len(lines):
    line = lines[i].strip()

    if line == "!platform_table_begin":
        i += 1
        i += 1  # skip header
        while i < len(lines) and not lines[i].strip().startswith("!platform_table_end"):
            cols = lines[i].strip().split("\t")
            if len(cols) >= 2:
                pid = cols[0].strip('"')
                gene = cols[1].strip('"')
                if gene and "///" not in gene:
                    probe_to_gene[pid] = gene
                elif gene:
                    probe_to_gene[pid] = gene.split("///")[0].strip()
            i += 1

    if line == "!series_matrix_table_begin":
        i += 1
        matrix_header = lines[i].strip().split("\t")
        i += 1
        while i < len(lines) and not lines[i].strip().startswith("!series_matrix_table_end"):
            cols = lines[i].strip().split("\t")
            pid = cols[0].strip('"')
            vals = []
            for v in cols[1:]:
                try:
                    vals.append(float(v.strip('"')))
                except:
                    vals.append(float('nan'))
            if len(vals) == len(samples):
                gene = probe_to_gene.get(pid, pid)
                if gene not in expr_data:
                    expr_data[gene] = vals
                else:
                    old_m = mean(expr_data[gene])
                    new_m = mean(vals)
                    if new_m > old_m:
                        expr_data[gene] = vals
            i += 1
    i += 1

print(f"Genes with expression: {len(expr_data)}")

# Run Welch t-test HF vs NF
print(f"\nRunning differential expression: {len(hf_idx)} HF vs {len(nf_idx)} NF...")
deg_results = []

for gene, values in expr_data.items():
    hf_vals = [values[i] for i in hf_idx]
    nf_vals = [values[i] for i in nf_idx]

    mean_hf = mean(hf_vals)
    mean_nf = mean(nf_vals)
    log2fc = mean_hf - mean_nf  # RMA data is in log2 space

    try:
        t_stat, p_val = stats.ttest_ind(hf_vals, nf_vals, equal_var=False)
    except:
        p_val = 1.0

    deg_results.append({
        "gene": gene,
        "log2FC": log2fc,
        "mean_HF": mean_hf,
        "mean_NF": mean_nf,
        "p_value": p_val,
    })

deg_results.sort(key=lambda x: x["p_value"])
n_sig = sum(1 for d in deg_results if d["p_value"] < 0.05)
print(f"Genes with P < 0.05: {n_sig}/{len(deg_results)}")

# Load reference DEGs (GSE193997, mouse -> human uppercased)
print("\nLoading reference DEGs (GSE193997, padj < 0.05)...")
ref_deg = {}
with open(r"D:\miri_paper\results\tables\GSE193997_DESeq2_Sham_vs_IR6h.csv") as f:
    for row in csv.DictReader(f):
        try:
            if float(row["padj"]) < 0.05:
                gene = row["symbol"].strip().upper()
                ref_deg[gene] = {
                    "log2FC": float(row["log2FoldChange"]),
                    "padj": float(row["padj"]),
                }
        except:
            pass
print(f"Reference DEGs: {len(ref_deg)}")

# Direction consistency
gse_dict = {d["gene"]: d for d in deg_results}
common = set(ref_deg.keys()) & set(gse_dict.keys())
print(f"Common genes: {len(common)}")

same_dir = 0
opp_dir = 0
no_dir = 0
for gene in common:
    ref_fc = ref_deg[gene]["log2FC"]
    gse_fc = gse_dict[gene]["log2FC"]
    if abs(gse_fc) < 0.01:
        no_dir += 1
    elif (ref_fc > 0 and gse_fc > 0) or (ref_fc < 0 and gse_fc < 0):
        same_dir += 1
    else:
        opp_dir += 1

total_with_dir = same_dir + opp_dir
pct = 100 * same_dir / total_with_dir if total_with_dir > 0 else 0
from scipy.stats import binomtest
bt = binomtest(same_dir, total_with_dir, p=0.5, alternative="greater")

print(f"\n{'='*50}")
print(f"DIRECTION CONSISTENCY")
print(f"{'='*50}")
print(f"Same direction:   {same_dir}/{total_with_dir} ({pct:.1f}%)")
print(f"Opposite direction: {opp_dir}/{total_with_dir}")
print(f"No change in GSE: {no_dir}")
print(f"Binomial test P = {bt.pvalue:.8f}")
if bt.pvalue < 0.001:
    print(f">>> HIGHLY SIGNIFICANT direction preservation <<<")

# Core genes
print(f"\n--- Core Gene Direction Check ---")
core = ["OGT", "OGA", "NFE2L2", "HMOX1", "SLC7A11", "GPX4",
        "GCLC", "RELA", "KEAP1", "BACH1", "FTH1", "TFRC",
        "SRXN1", "NQO1", "GCLM", "HIF1A", "STAT3"]
found = 0; matched = 0
for gene in core:
    ref = ref_deg.get(gene, {})
    gse = gse_dict.get(gene, {})
    ref_fc = ref.get("log2FC", None)
    gse_fc = gse.get("log2FC", None)
    if ref_fc and gse_fc:
        found += 1
        match = "OK" if ref_fc * gse_fc > 0 else "MISMATCH"
        if ref_fc * gse_fc > 0: matched += 1
        print(f"  {gene:10s}  Ref:{ref_fc:+6.2f}  GSE:{gse_fc:+6.2f}  P_gse={gse['p_value']:.3f}  {match}")
    elif ref_fc:
        print(f"  {gene:10s}  Ref:{ref_fc:+6.2f}  GSE: not found on array")
    elif gse_fc:
        print(f"  {gene:10s}  Ref: not DEG  GSE:{gse_fc:+6.2f}")

print(f"\nCore genes found: {found}/{len(core)}, direction match: {matched}/{found}")

# Save results
out_csv = OUT / "GSE57338_direction_consistency.csv"
with open(out_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene", "log2FC_GSE57338", "p_GSE57338",
                     "log2FC_GSE193997", "padj_GSE193997", "same_direction"])
    for gene in common:
        ref_fc = ref_deg[gene]["log2FC"]
        gse_fc = gse_dict[gene]["log2FC"]
        same = (ref_fc > 0 and gse_fc > 0) or (ref_fc < 0 and gse_fc < 0)
        writer.writerow([gene, gse_fc, gse_dict[gene]["p_value"],
                         ref_fc, ref_deg[gene]["padj"], same])
print(f"\nSaved: {out_csv}")

# GSEA ranking for external fgsea
gsea_data = []
for d in deg_results:
    if d["p_value"] > 0 and d["p_value"] < 1:
        signed_score = d["log2FC"] * (-np.log10(d["p_value"]))
        gsea_data.append((d["gene"], signed_score))
gsea_data.sort(key=lambda x: x[1], reverse=True)

out_rnk = OUT / "GSE57338_GSEA_ranking.rnk"
with open(out_rnk, "w") as f:
    for gene, score in gsea_data:
        f.write(f"{gene}\t{score}\n")
print(f"GSEA ranking saved: {out_rnk} ({len(gsea_data)} genes)")

print(f"\n=== DONE ===")
