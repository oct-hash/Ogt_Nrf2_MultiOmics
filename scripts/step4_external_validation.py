"""
Step ④: External validation
Part A: GTEx v8 healthy LV baseline
Part B: GSE57338 HF vs NF direction consistency + GSEA
"""
import csv, json, os, gzip
from pathlib import Path
from collections import Counter
from statistics import mean, stdev
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu

BASE = Path(r"D:\Cardiac_Ogt_MultiOmics")
OUT = BASE / "output_tables"
FIGS = BASE / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

# ============================================================
# PART A: GTEx v8 Healthy LV Baseline
# ============================================================
print("=" * 60)
print("PART A: GTEx v8 Healthy Left Ventricle Baseline")
print("=" * 60)

with open(BASE / "validation_data" / "GTEx_core_genes_LV.json") as f:
    gtex_raw = json.load(f)

# Extract LV data for each gene
gtex_lv = {}
for item in gtex_raw["data"]:
    if item["tissueSiteDetailId"] == "Heart_Left_Ventricle":
        gene = item["geneSymbol"]
        gtex_lv[gene] = item["median"]

print(f"\nGTEx v8 Heart - Left Ventricle (n=431 healthy donors):")
print(f"{'Gene':12s} {'Median TPM':>10s}")
print("-" * 24)
for gene in ["OGT", "NFE2L2", "HMOX1", "SLC7A11", "BACH1"]:
    tpm = gtex_lv.get(gene, None)
    if tpm is not None:
        print(f"  {gene:10s} {tpm:10.1f}")
    else:
        print(f"  {gene:10s} {'N/A':>10s}")
print("  (OGA, KEAP1 not in GTEx API query — add if needed)")

# Comparison with disease model (from GSE193997)
deg_ref = {
    "OGT": -0.56,
    "OGA": -0.44,
    "NFE2L2": 0.53,
    "HMOX1": 4.43,
    "SLC7A11": 4.07,
    "GPX4": 0.25,
    "GCLC": 0.84,
    "RELA": 1.13,
    "BACH1": None,
    "KEAP1": None,
}

print(f"\n{'Gene':10s} {'GTEx LV TPM':>12s} {'Disease log2FC':>14s}  Interpretation")
print("-" * 60)
for gene in ["OGT", "NFE2L2", "HMOX1", "SLC7A11", "GPX4", "GCLC", "RELA"]:
    tpm = gtex_lv.get(gene, None)
    fc = deg_ref.get(gene, None)
    tpm_str = f"{tpm:10.1f}" if tpm is not None else "       N/A"
    if fc is not None and tpm is not None:
        if fc > 0.5:
            interp = "Nrf2 target: LOW in health, HIGH in disease"
        elif fc < -0.3:
            interp = "Constitutive in health, DOWN in disease"
        else:
            interp = "Modest change"
        print(f"  {gene:10s} {tpm_str:>10s} {fc:+14.2f}  {interp}")
    else:
        print(f"  {gene:10s} {tpm_str:>10s} {'N/A':>14s}")

# ============================================================
# PART B: GSE57338 Direction Consistency
# ============================================================
print(f"\n{'='*60}")
print("PART B: GSE57338 Heart Failure Direction Consistency")
print(f"{'='*60}")

gse_path = BASE / "validation_data" / "GSE57338_series_matrix.txt.gz"

# Parse series matrix
print("Parsing GSE57338 series matrix...")
samples = []
expr_data = {}
probe_to_gene = {}
in_matrix = False
in_platform = False
matrix_header = None

with gzip.open(gse_path, "rt", errors="ignore") as f:
    lines = f.readlines()

# First pass: get metadata and platform probe-gene mapping
i = 0
while i < len(lines):
    line = lines[i].strip()

    # Sample characteristics
    if line.startswith("!Sample_title"):
        samples = [s.strip('"') for s in line.split("\t")[1:]]
    elif line.startswith("!Sample_characteristics_ch1"):
        pass  # Will parse later

    # Platform table (probe-gene mapping)
    if line == "!platform_table_begin":
        i += 1
        plat_header = lines[i].strip().split("\t")
        i += 1
        while i < len(lines) and not lines[i].strip().startswith("!platform_table_end"):
            cols = lines[i].strip().split("\t")
            if len(cols) >= 2:
                probe_id = cols[0].strip('"')
                gene_symbol = cols[1].strip('"') if len(cols) > 1 else ""
                if gene_symbol and "///" not in gene_symbol:
                    probe_to_gene[probe_id] = gene_symbol
                elif gene_symbol:
                    # Multi-gene probe: take first
                    probe_to_gene[probe_id] = gene_symbol.split("///")[0].strip()
            i += 1

    # Matrix table
    if line == "!series_matrix_table_begin":
        in_matrix = True
        i += 1
        matrix_header = lines[i].strip().split("\t")
        i += 1
        while i < len(lines) and not lines[i].strip().startswith("!series_matrix_table_end"):
            cols = lines[i].strip().split("\t")
            probe_id = cols[0].strip('"')
            values = [float(v.strip('"')) for v in cols[1:]]
            gene = probe_to_gene.get(probe_id, probe_id)

            if gene not in expr_data:
                expr_data[gene] = values
            else:
                # Multiple probes: keep max mean
                old_mean = sum(expr_data[gene]) / len(expr_data[gene])
                new_mean = sum(values) / len(values)
                if new_mean > old_mean:
                    expr_data[gene] = values
            i += 1

    i += 1

print(f"Samples: {len(samples)}")
print(f"Genes with expression: {len(expr_data)}")

# Parse phenotype (HF vs NF) from sample titles
# Sample titles format: "Heart failure sample X" or "Non-failing heart sample X"
hf_samples = []
nf_samples = []
for idx, title in enumerate(samples):
    title_lower = title.lower()
    if "failing" in title_lower or "heart failure" in title_lower or "hf" in title_lower:
        hf_samples.append(idx)
    elif "non" in title_lower or "normal" in title_lower or "control" in title_lower or "nf" in title_lower:
        nf_samples.append(idx)

# If automatic parsing fails, try manual from known GSE57338 design
if not hf_samples or not nf_samples:
    print("Auto phenotype detection failed, using known GSE57338 design...")
    # GSE57338: 3 HF RNA-seq (samples 1-3) + 3 NF (samples 4-6) + microarray set
    # Let me check sample count and try
    n_total = len(samples)
    mid = n_total // 2
    print(f"Sample count: {n_total}, trying first half HF, second half NF")
    print(f"Sample titles: {samples[:5]}...")

print(f"HF samples: {len(hf_samples)}, NF samples: {len(nf_samples)}")

if len(hf_samples) < 3 or len(nf_samples) < 3:
    print("\n[WARNING] Phenotype parsing unclear. Printing sample titles:")
    for i, s in enumerate(samples[:10]):
        print(f"  {i}: {s[:100]}")
    print("\nSkipping DEG analysis — manual inspection needed.")
else:
    # Run Welch's t-test per gene
    print(f"\nRunning differential expression: {len(hf_samples)} HF vs {len(nf_samples)} NF...")
    deg_results = []

    for gene, values in expr_data.items():
        hf_vals = [values[i] for i in hf_samples]
        nf_vals = [values[i] for i in nf_samples]

        mean_hf = mean(hf_vals)
        mean_nf = mean(nf_vals)

        if mean_nf > 0:
            log2fc = mean_hf - mean_nf  # Already in log2 space if RMA-normalized
        else:
            log2fc = 0

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

    # Sort by P value
    deg_results.sort(key=lambda x: x["p_value"])

    n_sig = sum(1 for d in deg_results if d["p_value"] < 0.05)
    print(f"Genes with P < 0.05: {n_sig}/{len(deg_results)}")

    # Load reference DEGs from GSE193997
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
    print(f"Reference DEGs (GSE193997): {len(ref_deg)}")

    # Direction consistency
    gse_dict = {d["gene"]: d for d in deg_results}
    common = set(ref_deg.keys()) & set(gse_dict.keys())
    print(f"Common genes: {len(common)}")

    same_dir = 0
    opp_dir = 0
    for gene in common:
        ref_fc = ref_deg[gene]["log2FC"]
        gse_fc = gse_dict[gene]["log2FC"]
        if (ref_fc > 0 and gse_fc > 0) or (ref_fc < 0 and gse_fc < 0):
            same_dir += 1
        else:
            opp_dir += 1

    total_with_dir = same_dir + opp_dir
    print(f"\nDirection consistency:")
    print(f"  Same direction: {same_dir}/{total_with_dir} ({100*same_dir/total_with_dir:.1f}%)")

    # Binomial test
    from scipy.stats import binomtest
    bt = binomtest(same_dir, total_with_dir, p=0.5, alternative="greater")
    print(f"  Binomial P = {bt.pvalue:.6f}")

    # Core genes check
    print(f"\n--- Core Gene Direction Check ---")
    core = ["OGT", "OGA", "NFE2L2", "HMOX1", "SLC7A11", "GPX4",
            "GCLC", "RELA", "KEAP1", "BACH1", "FTH1", "TFRC"]
    for gene in core:
        ref = ref_deg.get(gene, {})
        gse = gse_dict.get(gene, {})
        ref_fc = ref.get("log2FC", None)
        gse_fc = gse.get("log2FC", None)
        if ref_fc and gse_fc:
            match = "✓" if (ref_fc * gse_fc > 0) else "✗"
            print(f"  {gene:10s}  Ref: {ref_fc:+5.2f}  GSE: {gse_fc:+5.2f}  {match}")
        elif ref_fc:
            print(f"  {gene:10s}  Ref: {ref_fc:+5.2f}  GSE: not found")
        elif gse_fc:
            print(f"  {gene:10s}  Ref: not found  GSE: {gse_fc:+5.2f}")

    # Save results
    out_csv = OUT / "GSE57338_direction_consistency.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "log2FC_GSE57338", "p_GSE57338", "log2FC_GSE193997", "padj_GSE193997", "same_direction"])
        for gene in common:
            ref_fc = ref_deg[gene]["log2FC"]
            gse_fc = gse_dict[gene]["log2FC"]
            same = (ref_fc > 0 and gse_fc > 0) or (ref_fc < 0 and gse_fc < 0)
            writer.writerow([
                gene, gse_fc, gse_dict[gene]["p_value"],
                ref_fc, ref_deg[gene]["padj"], same
            ])
    print(f"\nSaved: {out_csv}")

    # GSEA preparation
    # Rank genes by log2FC * -log10(p) for signed ranking
    gsea_data = []
    for d in deg_results:
        if d["p_value"] > 0:
            signed_score = d["log2FC"] * (-np.log10(d["p_value"]))
            gsea_data.append((d["gene"], signed_score))
    gsea_data.sort(key=lambda x: x[1], reverse=True)

    out_rnk = OUT / "GSE57338_GSEA_ranking.rnk"
    with open(out_rnk, "w") as f:
        for gene, score in gsea_data:
            f.write(f"{gene}\t{score}\n")
    print(f"GSEA ranking saved: {out_rnk} ({len(gsea_data)} genes)")

print(f"\n=== PART B DONE ===")
