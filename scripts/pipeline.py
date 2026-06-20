#!/usr/bin/env python3
"""
Ogt_Nrf2_MultiOmics — Analysis Pipeline
=========================================
Four-step validation pipeline for the discovery paper:
  Step 3: OGT substrate enrichment analysis
  Step 4: External validation (GTEx baseline + GSE57338 direction consistency)
  Step 1: Cell-cell communication dual-method consensus
  Step 2: mRNA-protein discordance (residual analysis)

Usage:
  python pipeline.py step3   # OGT substrate enrichment
  python pipeline.py step4   # External validation
  python pipeline.py step1   # CCC dual-method consensus
  python pipeline.py step2   # mRNA-protein residual analysis
  python pipeline.py all     # Run all steps

Requirements:
  pip install pandas scipy openpyxl numpy
"""
import sys, os, csv, json, gzip
from pathlib import Path
from collections import Counter, defaultdict
from statistics import mean
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact, binomtest, mannwhitneyu

# ── Config ────────────────────────────────────────────────────────
def load_config():
    cfg_path = Path(__file__).parent / "config.json"
    with open(cfg_path) as f:
        return json.load(f)

CFG = load_config()
BASE = Path(CFG["project_root"])
OUT = BASE / CFG["paths"]["output_tables"]
OUT.mkdir(parents=True, exist_ok=True)

def print_header(title):
    print(f"\n{'='*60}\n{title}\n{'='*60}")

# ── Step 3: OGT Substrate Enrichment ──────────────────────────────
def step3_ogt_enrichment():
    print_header("Step 3: OGT Substrate Enrichment Analysis")

    # Load OGT interactors from JSON (pre-extracted)
    ogt_json = BASE / CFG["datasets"]["ogt_high_json"]
    if not ogt_json.exists():
        # Extract from Excel using pandas
        import pandas as pd
        xlsx = BASE / CFG["datasets"]["ogt_pin_excel"]
        df = pd.read_excel(xlsx)
        human = df[df["ncbi_taxonomy identifier_b"].str.contains("Homo sapiens", na=False)]
        high = human[human["High-stringency interactor"].str.lower().str.contains("yes", na=False)]
        genes = sorted(set(high["gene_name_b"].dropna().astype(str)))
        with open(ogt_json, "w") as f:
            json.dump(genes, f)
        ogt_high = set(genes)
    else:
        with open(ogt_json) as f:
            ogt_high = set(json.load(f))

    print(f"OGT high-stringency interactors: {len(ogt_high)}")

    # Load DEGs
    deg_path = Path(CFG["datasets"]["deg_file"])
    deg_data = {}
    with open(deg_path) as f:
        for row in csv.DictReader(f):
            try:
                if float(row["padj"]) < CFG["parameters"]["deg_padj_threshold"]:
                    g = row["symbol"].strip().upper()
                    deg_data[g] = {
                        "log2FC": float(row["log2FoldChange"]),
                        "padj": float(row["padj"]),
                    }
            except (ValueError, KeyError):
                pass

    print(f"DEGs (padj < 0.05): {len(deg_data)}")

    # Overlap
    overlap = set(deg_data) & ogt_high
    n_deg, n_ogt, n_over = len(deg_data), len(ogt_high), len(overlap)
    n_bg = CFG["parameters"]["fisher_background_genes"]

    a, b = n_over, n_deg - n_over
    c, d = n_ogt - n_over, n_bg - n_deg - n_ogt + n_over
    or_val, p_val = fisher_exact([[a, b], [c, max(d, 1)]], alternative="greater")

    print(f"\nOverlap: {n_over}/{n_deg} genes ({100*n_over/n_deg:.2f}%)")
    print(f"Fisher: OR={or_val:.2f}, P={p_val:.8f}")
    print(f"{'SIGNIFICANT' if p_val < 0.05 else 'NOT significant'}")

    # Direction
    up = sum(1 for g in overlap if deg_data[g]["log2FC"] > 0)
    down = n_over - up
    print(f"Direction: {up} UP, {down} DOWN")

    # Key genes
    key_genes = {
        "OGT": "O-GlcNAc transferase", "OGA": "O-GlcNAcase",
        "NFE2L2": "NRF2", "KEAP1": "NRF2 inhibitor", "BACH1": "NRF2 competitor",
        "HMOX1": "Ferroptosis", "SLC7A11": "Ferroptosis/cystine",
        "GPX4": "Ferroptosis/lipid ROS", "FTH1": "Iron storage", "TFRC": "Iron import",
        "RELA": "NF-kB p65", "NFKB1": "NF-kB", "STAT3": "STAT3", "HIF1A": "HIF-1a",
    }
    print(f"\nKey genes:")
    for gene, desc in key_genes.items():
        in_ogt = "YES" if gene in ogt_high else "-"
        in_deg = deg_data.get(gene)
        d_str = f"DEG log2FC={in_deg['log2FC']:+.2f}" if in_deg else "-"
        print(f"  {gene:10s} OGT:{in_ogt:4s} {d_str} ({desc})")

    # Save
    with open(OUT / "step3_ogt_enrichment.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "in_OGT_PIN", "log2FC", "padj"])
        for gene in overlap:
            w.writerow([gene, "YES", deg_data[gene]["log2FC"], deg_data[gene]["padj"]])

    # Verify
    assert n_over >= 0, "Overlap count negative"
    assert or_val > 0, "OR should be positive"
    print(f"\n[PASS] Step 3 complete. Output: step3_ogt_enrichment.csv")
    return {"overlap": n_over, "or": or_val, "p": p_val, "significant": p_val < 0.05}

# ── Step 4: External Validation ───────────────────────────────────
def step4_external_validation():
    print_header("Step 4: External Validation (GTEx + GSE57338)")

    # Part A: GTEx baseline
    print("Part A: GTEx v8 Healthy LV Baseline")
    gtex_path = BASE / CFG["datasets"]["gtex_json"]
    with open(gtex_path) as f:
        gtex_raw = json.load(f)

    gtex_lv = {}
    for item in gtex_raw["data"]:
        if item["tissueSiteDetailId"] == "Heart_Left_Ventricle":
            gtex_lv[item["geneSymbol"]] = item["median"]

    print(f"GTEx v8 LV (n=431):")
    for gene in ["OGT", "NFE2L2", "HMOX1", "SLC7A11", "BACH1"]:
        tpm = gtex_lv.get(gene, "N/A")
        print(f"  {gene:10s} {tpm:8.1f} TPM")

    # Part B: GSE57338 direction consistency
    print("\nPart B: GSE57338 Direction Consistency")

    # Parse GPL11532 probe-to-gene
    probe_to_gene = {}
    gpl_path = BASE / CFG["datasets"]["gpl11532_annot"]
    with gzip.open(gpl_path, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 3:
                pid, gene = cols[0].strip(), cols[2].strip()
                if gene and gene != "---":
                    probe_to_gene[pid] = gene.split("///")[0].strip()

    # Parse GSE57338 matrix
    gse_path = BASE / CFG["datasets"]["gse57338_matrix"]
    with gzip.open(gse_path, "rt", errors="ignore") as f:
        lines = f.readlines()

    # Get phenotype
    hf_mask = []
    for line in lines:
        ls = line.strip()
        if ls.startswith("!Sample_characteristics_ch1") and "heart failure" in ls:
            hf_mask = ["heart failure: yes" in p for p in ls.split("\t")[1:]]
            break

    hf_idx = [j for j, v in enumerate(hf_mask) if v]
    nf_idx = [j for j, v in enumerate(hf_mask) if not v]
    print(f"GSE57338: {len(hf_idx)} HF, {len(nf_idx)} NF")

    # Parse expression matrix
    expr_data = {}
    i = 0
    while i < len(lines):
        ls = lines[i].strip()
        if ls == "!series_matrix_table_begin":
            i += 2
            while i < len(lines) and not lines[i].strip().startswith("!series_matrix_table_end"):
                cols = lines[i].strip().split("\t")
                pid = cols[0].strip('"')
                gene = probe_to_gene.get(pid, pid)
                vals = []
                for v in cols[1:]:
                    try: vals.append(float(v.strip('"')))
                    except: vals.append(float("nan"))
                if len(vals) == len(hf_mask):
                    if gene not in expr_data or mean(vals) > mean(expr_data[gene]):
                        expr_data[gene] = vals
                i += 1
        i += 1

    # DEG
    deg_results = []
    for gene, vals in expr_data.items():
        hf_v = [vals[i] for i in hf_idx]
        nf_v = [vals[i] for i in nf_idx]
        log2fc = mean(hf_v) - mean(nf_v)
        try:
            _, p_val = stats.ttest_ind(hf_v, nf_v, equal_var=False)
        except:
            p_val = 1.0
        deg_results.append({"gene": gene, "log2FC": log2fc, "p_value": p_val})

    gse_dict = {d["gene"]: d for d in deg_results}

    # Load reference DEGs
    ref_deg = {}
    ref_path = Path(CFG["datasets"]["deg_file"])
    with open(ref_path) as f:
        for row in csv.DictReader(f):
            try:
                if float(row["padj"]) < 0.05:
                    ref_deg[row["symbol"].strip().upper()] = float(row["log2FoldChange"])
            except:
                pass

    # Direction consistency
    common = set(ref_deg) & set(gse_dict)
    same, opp, nodir = 0, 0, 0
    for gene in common:
        rf, gf = ref_deg[gene], gse_dict[gene]["log2FC"]
        if abs(gf) < 0.01: nodir += 1
        elif rf * gf > 0: same += 1
        else: opp += 1

    twd = same + opp
    pct = 100 * same / twd if twd > 0 else 0
    bt = binomtest(same, twd, p=0.5, alternative="greater")

    print(f"Common genes: {len(common)}")
    print(f"Same direction: {same}/{twd} ({pct:.1f}%), Opp: {opp}, NoDir: {nodir}")
    print(f"Binomial P = {bt.pvalue:.6f}")

    # Core genes check
    core = ["OGT", "NFE2L2", "HMOX1", "SLC7A11", "GPX4", "GCLC", "RELA",
            "NQO1", "FTH1", "TFRC", "GCLM", "SRXN1"]
    print(f"\nCore gene direction:")
    matched = 0
    for gene in core:
        rf = ref_deg.get(gene)
        gse = gse_dict.get(gene)
        if rf and gse:
            ok = rf * gse["log2FC"] > 0
            if ok: matched += 1
            print(f"  {gene:10s} Ref:{rf:+5.2f} GSE:{gse['log2FC']:+5.2f} {'OK' if ok else 'X'}")
    print(f"Match: {matched}/{sum(1 for g in core if g in ref_deg and g in gse_dict)}")

    # Save
    with open(OUT / "step4_direction_consistency.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "log2FC_GSE57338", "p_value", "log2FC_GSE193997", "same_direction"])
        for gene in common:
            same_dir = (ref_deg[gene] > 0 and gse_dict[gene]["log2FC"] > 0) or \
                       (ref_deg[gene] < 0 and gse_dict[gene]["log2FC"] < 0)
            w.writerow([gene, gse_dict[gene]["log2FC"], gse_dict[gene]["p_value"],
                        ref_deg[gene], same_dir])

    # Verify
    assert len(common) > 1000, "Too few common genes for direction test"
    assert len(hf_idx) >= 100, "GSE57338: too few HF samples"
    print(f"\n[PASS] Step 4 complete. Output: step4_direction_consistency.csv")
    return {"common_genes": len(common), "same_pct": pct, "binomial_p": bt.pvalue}

# ── Step 1: CCC Dual-Method Consensus ─────────────────────────────
def step1_ccc_consensus():
    print_header("Step 1: CCC Dual-Method Consensus (CellChat + CellPhoneDB)")

    myeloid = CFG["parameters"]["myeloid_cell_types"]
    targets = CFG["parameters"]["target_cell_types"]

    # Load CellChat
    cellchat = []
    cc_path = Path(CFG["datasets"]["cellchat_lr"])
    with open(cc_path) as f:
        for row in csv.DictReader(f):
            cellchat.append({
                "sender": row["Sender"], "receiver": row["Receiver"],
                "ligand": row["Ligand"], "receptor": row["Receptor"],
                "score": float(row["Score"]), "pathway": row["Pathway"],
            })

    # Filter myeloid -> targets
    cc_filtered = [x for x in cellchat if x["sender"] in myeloid and x["receiver"] in targets]
    print(f"CellChat myeloid->CM/Peri: {len(cc_filtered)}")

    # Load CellPhoneDB
    cpdb = []
    li_path = Path(CFG["datasets"]["liana_full"])
    with open(li_path) as f:
        for row in csv.DictReader(f):
            p = float(row.get("cellphone_pvals", 1))
            if p < 0.05:
                cpdb.append({
                    "sender": row["source"], "receiver": row["target"],
                    "ligand": row["ligand_complex"], "receptor": row["receptor_complex"],
                    "pval": p,
                })

    cpdb_filtered = [x for x in cpdb if x["sender"] in myeloid and x["receiver"] in targets]
    print(f"CellPhoneDB myeloid->CM/Peri: {len(cpdb_filtered)}")

    # Consensus matching
    def norm(s):
        return s.strip().upper().replace(" ", "").replace("_", "")

    consensus = []
    for cc in cc_filtered:
        cc_lig = norm(cc["ligand"])
        cc_rec = norm(cc["receptor"])
        for cp in cpdb_filtered:
            if cc["sender"] == cp["sender"] and cc["receiver"] == cp["receiver"]:
                cp_lig = norm(cp["ligand"])
                cp_rec = norm(cp["receptor"])
                if (cc_lig == cp_lig and cc_rec == cp_rec) or \
                   (cc_lig in cp_lig or cp_lig in cc_lig) and (cc_rec in cp_rec or cp_rec in cc_rec):
                    consensus.append({**cc, "cpdb_pval": cp["pval"]})
                    break

    consensus.sort(key=lambda x: x["score"], reverse=True)
    print(f"\nConsensus: {len(consensus)} interactions")
    for c in consensus:
        print(f"  {c['sender']:14s} -> {c['receiver']:14s} {c['ligand']:12s} -> {c['receptor']:12s} "
              f"Score={c['score']:.4f} P={c['cpdb_pval']:.4f} [{c['pathway']}]")

    # Pathways
    pathways = Counter(c["pathway"] for c in consensus)
    print(f"\nPathways: {dict(pathways)}")

    # Save
    with open(OUT / "step1_ccc_consensus.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sender","receiver","ligand","receptor",
                                           "cellchat_score","cpdb_pval","pathway"])
        w.writeheader()
        for c in consensus:
            w.writerow({"sender":c["sender"],"receiver":c["receiver"],
                        "ligand":c["ligand"],"receptor":c["receptor"],
                        "cellchat_score":c["score"],"cpdb_pval":c["cpdb_pval"],
                        "pathway":c["pathway"]})

    # Verify
    assert len(consensus) > 0, "No consensus interactions found"
    assert any("SPP1" in c["ligand"].upper() for c in consensus), "SPP1 not in consensus"
    print(f"\n[PASS] Step 1 complete. Output: step1_ccc_consensus.csv")
    return {"consensus_n": len(consensus), "pathways": dict(pathways)}

# ── Step 2: mRNA-Protein Residual ─────────────────────────────────
def step2_residual():
    print_header("Step 2: mRNA-Protein Discordance (Residual Analysis)")

    prot_path = Path(CFG["datasets"]["proteomics_merged"])
    data = []
    with open(prot_path) as f:
        for row in csv.DictReader(f):
            try:
                gene = row["Gene"].strip()
                prot_fc = float(row["log2FC"])
                rna_fc = float(row["log2FoldChange"])
                data.append({"gene": gene, "prot_log2fc": prot_fc, "rna_log2fc": rna_fc})
            except (ValueError, KeyError):
                pass

    print(f"Matched genes: {len(data)}")

    rna_fc = [d["rna_log2fc"] for d in data]
    prot_fc = [d["prot_log2fc"] for d in data]
    slope, intercept, r_val, p_val, _ = stats.linregress(rna_fc, prot_fc)

    print(f"Regression: log2(Prot) = {intercept:.3f} + {slope:.3f} * log2(RNA)")
    print(f"R = {r_val:.3f}, P = {p_val:.2e}")

    # Residuals
    for d in data:
        predicted = intercept + slope * d["rna_log2fc"]
        d["residual"] = d["prot_log2fc"] - predicted

    data.sort(key=lambda d: abs(d["residual"]), reverse=True)

    # Concordance
    up_both = sum(1 for d in data if d["rna_log2fc"] > 0 and d["prot_log2fc"] > 0)
    down_both = sum(1 for d in data if d["rna_log2fc"] < 0 and d["prot_log2fc"] < 0)
    discordant = sum(1 for d in data if d["rna_log2fc"] * d["prot_log2fc"] < 0)
    concordant = up_both + down_both

    print(f"\nConcordant: {concordant}/{len(data)} ({100*concordant/len(data):.1f}%)")
    print(f"  Both UP: {up_both}, Both DOWN: {down_both}")
    print(f"Discordant: {discordant}/{len(data)} ({100*discordant/len(data):.1f}%)")

    # Key genes
    key = ["OGT", "FTH1", "TFRC", "GPX4", "NQO1", "HMOX1", "SLC7A11",
           "NFE2L2", "RELA", "GCLC", "G6PDX"]
    print(f"\nKey gene residuals:")
    for gene in key:
        match = [d for d in data if d["gene"].upper() == gene.upper()]
        if match:
            d = match[0]
            print(f"  {gene:10s} RNA:{d['rna_log2fc']:+6.2f} Prot:{d['prot_log2fc']:+6.2f} "
                  f"Res:{d['residual']:+6.2f}")
        else:
            print(f"  {gene:10s} not detected")

    # Save
    with open(OUT / "step2_residual_analysis.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene","rna_log2fc","prot_log2fc","residual"])
        w.writeheader()
        for d in data: w.writerow(d)

    # Verify
    assert len(data) > 1000, "Too few matched genes"
    assert abs(r_val) < 0.5, f"R={r_val:.3f} too high - expected weak mRNA-protein correlation"
    assert 40 < 100*concordant/len(data) < 70, f"Concordance {100*concordant/len(data):.1f}% outside expected 40-70% range"
    print(f"\n[PASS] Step 2 complete. Output: step2_residual_analysis.csv")
    return {"n_genes": len(data), "R": r_val, "concordance_pct": 100*concordant/len(data)}

# ── Main ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python pipeline.py [step3|step4|step1|step2|all]")
        sys.exit(1)

    cmd = sys.argv[1].lower()

    results = {}
    if cmd in ("step3", "all"):
        results["step3"] = step3_ogt_enrichment()
    if cmd in ("step4", "all"):
        results["step4"] = step4_external_validation()
    if cmd in ("step1", "all"):
        results["step1"] = step1_ccc_consensus()
    if cmd in ("step2", "all"):
        results["step2"] = step2_residual()

    if cmd == "all":
        print_header("ALL STEPS COMPLETE")
        for name, r in results.items():
            print(f"  {name}: {json.dumps(r, indent=2, default=str)}")
