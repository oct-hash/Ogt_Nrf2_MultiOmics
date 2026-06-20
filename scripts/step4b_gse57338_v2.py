"""
Step ④ Part B v2: GSE57338 with proper GPL11532 probe-to-gene mapping
"""
import csv, gzip
from pathlib import Path
from statistics import mean
import numpy as np
from scipy import stats
from scipy.stats import binomtest

BASE = Path(r"D:\Cardiac_Ogt_MultiOmics")
OUT = BASE / "output_tables"
OUT.mkdir(parents=True, exist_ok=True)

# ── 1. Parse GPL11532 probe-to-gene ───────────────────────────────
print("Parsing GPL11532 annotation...")
probe_to_gene = {}
with gzip.open(BASE / "validation_data" / "GPL11532.annot.gz", "rt", errors="ignore") as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if len(cols) >= 3:
            probe_id = cols[0].strip()
            gene_symbol = cols[2].strip()
            if gene_symbol and gene_symbol != "---" and "///" not in gene_symbol:
                probe_to_gene[probe_id] = gene_symbol
            elif gene_symbol and "///" in gene_symbol:
                probe_to_gene[probe_id] = gene_symbol.split("///")[0].strip()
print(f"Probes mapped to genes: {len(probe_to_gene)}")

# ── 2. Parse GSE57338 matrix ──────────────────────────────────────
print("Parsing GSE57338 expression matrix...")
with gzip.open(BASE / "validation_data" / "GSE57338_series_matrix.txt.gz", "rt", errors="ignore") as f:
    lines = f.readlines()

# Phenotype
hf_mask = []
for line in lines:
    ls = line.strip()
    if ls.startswith("!Sample_characteristics_ch1") and "heart failure" in ls:
        parts = ls.split("\t")
        hf_mask = ["heart failure: yes" in p for p in parts[1:]]
        break

hf_idx = [j for j, v in enumerate(hf_mask) if v]
nf_idx = [j for j, v in enumerate(hf_mask) if not v]
print(f"HF: {len(hf_idx)}, NF: {len(nf_idx)}")

# Matrix
expr_data = {}
i = 0
while i < len(lines):
    ls = lines[i].strip()
    if ls == "!series_matrix_table_begin":
        i += 1
        i += 1  # skip header
        while i < len(lines) and not lines[i].strip().startswith("!series_matrix_table_end"):
            cols = lines[i].strip().split("\t")
            pid = cols[0].strip('"')
            gene = probe_to_gene.get(pid, pid)
            vals = []
            for v in cols[1:]:
                try:
                    vals.append(float(v.strip('"')))
                except:
                    vals.append(float("nan"))
            if len(vals) == len(hf_mask):
                if gene not in expr_data:
                    expr_data[gene] = vals
                else:
                    if mean(vals) > mean(expr_data[gene]):
                        expr_data[gene] = vals
            i += 1
    i += 1

print(f"Unique genes: {len(expr_data)}")
print(f"Sample genes: {list(expr_data.keys())[:10]}")

# ── 3. DEG analysis ───────────────────────────────────────────────
print(f"\nDifferential expression: {len(hf_idx)} HF vs {len(nf_idx)} NF...")
deg_results = []
for gene, vals in expr_data.items():
    hf_v = [vals[i] for i in hf_idx]
    nf_v = [vals[i] for i in nf_idx]
    log2fc = mean(hf_v) - mean(nf_v)
    try:
        t_stat, p_val = stats.ttest_ind(hf_v, nf_v, equal_var=False)
    except:
        p_val = 1.0
    deg_results.append({"gene": gene, "log2FC": log2fc,
                         "mean_HF": mean(hf_v), "mean_NF": mean(nf_v),
                         "p_value": p_val})

deg_results.sort(key=lambda x: x["p_value"])
n_sig = sum(1 for d in deg_results if d["p_value"] < 0.05)
print(f"P < 0.05: {n_sig}/{len(deg_results)}")

# ── 4. Direction consistency with GSE193997 ─────────────────────────
print("\nLoading reference DEGs (GSE193997, uppercased -> human)...")
ref_deg = {}
with open(r"D:\miri_paper\results\tables\GSE193997_DESeq2_Sham_vs_IR6h.csv") as f:
    for row in csv.DictReader(f):
        try:
            if float(row["padj"]) < 0.05:
                gene = row["symbol"].strip().upper()
                ref_deg[gene] = {"log2FC": float(row["log2FoldChange"]),
                                  "padj": float(row["padj"])}
        except: pass
print(f"Reference DEGs: {len(ref_deg)}")

gse_dict = {d["gene"]: d for d in deg_results}
common = set(ref_deg.keys()) & set(gse_dict.keys())
print(f"Common genes (intersection): {len(common)}")

if len(common) == 0:
    print("\nERROR: Still no overlap. GSE57338 genes might use different symbol format.")
    print("GSE57338 sample genes:", sorted(list(gse_dict.keys()))[:20])
    print("REF sample genes:", sorted(list(ref_deg.keys()))[:20])
else:
    same = opp = nodir = 0
    for gene in common:
        rf = ref_deg[gene]["log2FC"]
        gf = gse_dict[gene]["log2FC"]
        if abs(gf) < 0.01: nodir += 1
        elif rf * gf > 0: same += 1
        else: opp += 1

    twd = same + opp
    pct = 100 * same / twd if twd > 0 else 0
    bt = binomtest(same, twd, p=0.5, alternative="greater")

    print(f"\n{'='*50}")
    print(f"DIRECTION CONSISTENCY")
    print(f"{'='*50}")
    print(f"Same: {same}/{twd} ({pct:.1f}%), Opp: {opp}, NoDir: {nodir}")
    print(f"Binomial P = {bt.pvalue:.8f}")
    print(f"{'>>> HIGHLY SIGNIFICANT' if bt.pvalue < 0.001 else ''}")

    # Core genes
    print(f"\n--- Core Gene Direction ---")
    core = ["OGT", "NFE2L2", "HMOX1", "SLC7A11", "GPX4", "GCLC", "RELA",
            "KEAP1", "BACH1", "FTH1", "TFRC", "SRXN1", "NQO1", "GCLM",
            "HIF1A", "STAT3", "OGA"]
    found = 0; matched = 0
    for gene in core:
        ref = ref_deg.get(gene, {})
        gse = gse_dict.get(gene, {})
        rf = ref.get("log2FC")
        gf = gse.get("log2FC")
        if rf and gf:
            found += 1
            ok = "OK" if rf * gf > 0 else "X"
            if rf * gf > 0: matched += 1
            print(f"  {gene:10s} Ref:{rf:+6.2f} GSE:{gf:+6.2f} P={gse['p_value']:.3f} {ok}")
        elif rf:
            print(f"  {gene:10s} Ref:{rf:+6.2f} GSE: not on array")
        elif gf:
            print(f"  {gene:10s} Ref: not DEG  GSE:{gf:+6.2f}")
    print(f"Found: {found}/{len(core)}, Match: {matched}/{found}")

    # Save
    with open(OUT / "GSE57338_direction_consistency.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "log2FC_GSE57338", "p_GSE57338",
                    "log2FC_GSE193997", "padj_GSE193997", "same_direction"])
        for gene in common:
            rf = ref_deg[gene]["log2FC"]
            gf = gse_dict[gene]["log2FC"]
            same_b = (rf > 0 and gf > 0) or (rf < 0 and gf < 0)
            w.writerow([gene, gf, gse_dict[gene]["p_value"],
                        rf, ref_deg[gene]["padj"], same_b])
    print(f"\nSaved: GSE57338_direction_consistency.csv")

    # GSEA ranking
    rnk = []
    for d in deg_results:
        if d["p_value"] > 0 and d["p_value"] < 1:
            rnk.append((d["gene"], d["log2FC"] * (-np.log10(d["p_value"]))))
    rnk.sort(key=lambda x: x[1], reverse=True)
    with open(OUT / "GSE57338_GSEA_ranking.rnk", "w") as f:
        for gene, score in rnk:
            f.write(f"{gene}\t{score}\n")
    print(f"GSEA ranking: {len(rnk)} genes")

print(f"\n=== DONE ===")
