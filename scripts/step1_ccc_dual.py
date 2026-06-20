"""
Step ①: CCC dual-method consensus (CellChat + CellPhoneDB)
Myeloid -> Cardiomyocyte/Pericyte paracrine interactions
"""
import csv
from pathlib import Path
from collections import defaultdict

BASE = Path(r"D:\Cardiac_Ogt_MultiOmics")
OUT = BASE / "output_tables"
OUT.mkdir(parents=True, exist_ok=True)

# ── 1. Load CellChat results ──────────────────────────────────────
print("=" * 60)
print("Step ①: Dual-Method CCC Consensus")
print("=" * 60)

cellchat_lr = []
with open(r"D:\miri_paper\results\tables\SuppTable_CellChat_LR_interactions.csv") as f:
    for row in csv.DictReader(f):
        sender = row["Sender"]
        receiver = row["Receiver"]
        ligand = row["Ligand"]
        receptor = row["Receptor"]
        score = float(row["Score"]) if row["Score"] else 0
        pathway = row["Pathway"]
        cellchat_lr.append({
            "sender": sender, "receiver": receiver,
            "ligand": ligand, "receptor": receptor,
            "score": score, "pathway": pathway,
        })

print(f"CellChat interactions: {len(cellchat_lr)}")
cellchat_set = set()
for x in cellchat_lr:
    key = (x["sender"], x["receiver"], x["ligand"], x["receptor"])
    cellchat_set.add(key)

# ── 2. Load CellPhoneDB results from LIANA ─────────────────────────
cpdb_lr = []
with open(r"D:\miri_paper\results\tables\LIANA_liana_res_full.csv") as f:
    for row in csv.DictReader(f):
        cpdb_p = float(row["cellphone_pvals"]) if row["cellphone_pvals"] else 1
        if cpdb_p < 0.05:  # significant in CellPhoneDB
            cpdb_lr.append({
                "sender": row["source"],
                "receiver": row["target"],
                "ligand": row["ligand_complex"],
                "receptor": row["receptor_complex"],
                "cpdb_pval": cpdb_p,
            })

print(f"CellPhoneDB significant interactions: {len(cpdb_lr)}")

cpdb_set = set()
for x in cpdb_lr:
    key = (x["sender"], x["receiver"], x["ligand"], x["receptor"])
    cpdb_set.add(key)

# ── 3. Filter myeloid->cardiomyocyte/pericyte ──────────────────────
myeloid_types = ["Macrophage", "Neutrophil"]
target_types = ["Cardiomyocyte", "Pericyte"]

# CellChat myeloid -> CM/Pericyte
cc_mye_target = []
for x in cellchat_lr:
    if x["sender"] in myeloid_types and x["receiver"] in target_types:
        cc_mye_target.append(x)

print(f"\nCellChat myeloid->CM/Pericyte: {len(cc_mye_target)}")

# CellPhoneDB myeloid -> CM/Pericyte
cpdb_mye_target = []
for x in cpdb_lr:
    if x["sender"] in myeloid_types and x["receiver"] in target_types:
        cpdb_mye_target.append(x)

print(f"CellPhoneDB myeloid->CM/Pericyte: {len(cpdb_mye_target)}")

# ── 4. Consensus: found in BOTH methods ────────────────────────────
# Match by (sender, receiver, ligand, receptor) - need fuzzy ligand/receptor matching
# because CellChat and CellPhoneDB use slightly different naming

def normalize(name):
    """Normalize gene/protein name for cross-method comparison"""
    return name.strip().upper().replace(" ", "").replace("_", "")

# Build CellChat index by sender-receiver
cc_index = defaultdict(list)
for x in cc_mye_target:
    key = (x["sender"], x["receiver"])
    cc_index[key].append(x)

# Build CellPhoneDB index
cpdb_index = defaultdict(list)
for x in cpdb_mye_target:
    key = (x["sender"], x["receiver"])
    cpdb_index[key].append(x)

# Find consensus: ligand/receptor must match (normalized) in both methods
consensus = []
for sr_key in set(list(cc_index.keys()) + list(cpdb_index.keys())):
    cc_entries = cc_index.get(sr_key, [])
    cpdb_entries = cpdb_index.get(sr_key, [])

    for cc in cc_entries:
        cc_lig = normalize(cc["ligand"])
        cc_rec = normalize(cc["receptor"])
        for cp in cpdb_entries:
            cp_lig = normalize(cp["ligand"])
            cp_rec = normalize(cp["receptor"])
            # Match: same ligand AND same receptor
            if cc_lig == cp_lig and cc_rec == cp_rec:
                consensus.append({
                    "sender": sr_key[0],
                    "receiver": sr_key[1],
                    "ligand": cc["ligand"],
                    "receptor": cc["receptor"],
                    "cellchat_score": cc["score"],
                    "cpdb_pval": cp["cpdb_pval"],
                    "pathway": cc.get("pathway", ""),
                })
            # Also match if ligand is substring of other (complex naming)
            elif (cc_lig in cp_lig or cp_lig in cc_lig) and (cc_rec in cp_rec or cp_rec in cc_rec):
                consensus.append({
                    "sender": sr_key[0],
                    "receiver": sr_key[1],
                    "ligand": cc["ligand"],
                    "receptor": cc["receptor"],
                    "ligand_cpdb": cp["ligand"],
                    "receptor_cpdb": cp["receptor"],
                    "cellchat_score": cc["score"],
                    "cpdb_pval": cp["cpdb_pval"],
                    "pathway": cc.get("pathway", ""),
                })

print(f"\n{'='*50}")
print(f"CONSENSUS (both methods): {len(consensus)} interactions")
print(f"{'='*50}")

# Sort by CellChat score
consensus.sort(key=lambda x: x["cellchat_score"], reverse=True)

# Show all consensus interactions
print(f"\n{'Sender':14s} {'Receiver':14s} {'Ligand':16s} {'Receptor':16s} {'CC Score':>8s} {'CPDB P':>8s}  Pathway")
print("-" * 100)
for c in consensus:
    print(f"{c['sender']:14s} {c['receiver']:14s} {c['ligand']:16s} {c['receptor']:16s} {c['cellchat_score']:8.4f} {c['cpdb_pval']:8.4f}  {c['pathway']}")

# ── 5. Nrf2/Ferroptosis relevance ─────────────────────────────────
ferr_ligands = {"HMOX1", "SLC7A11", "GPX4", "FTH1", "TFRC", "HIF1A"}
nrf2_ligands = {"NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "SRXN1", "TXN", "PRDX1"}
inflam_ligands = {"TNF", "IL1B", "IL6", "IL1A", "HMGB1", "S100A8", "S100A9", "TGFB1", "CCL2"}

print(f"\n--- Nrf2/Ferroptosis/Inflammatory Ligands in Consensus ---")
for c in consensus:
    lig = c["ligand"].upper()
    tags = []
    if any(g.upper() in lig for g in ferr_ligands): tags.append("FERROPTOSIS")
    if any(g.upper() in lig for g in nrf2_ligands): tags.append("NRF2")
    if any(g.upper() in lig for g in inflam_ligands): tags.append("INFLAMMATORY")
    if tags:
        print(f"  {c['sender']:14s} -> {c['receiver']:14s} {c['ligand']:16s} -> {c['receptor']:16s}  [{', '.join(tags)}]")

# ── 6. Pathway summary ────────────────────────────────────────────
print(f"\n--- Pathway Distribution in Consensus ---")
from collections import Counter
pathways = Counter(c["pathway"] for c in consensus)
for p, n in pathways.most_common():
    ligands = [c["ligand"] for c in consensus if c["pathway"] == p]
    print(f"  {p}: {n} interactions ({', '.join(ligands[:5])}{'...' if len(ligands)>5 else ''})")

# ── 7. Save ───────────────────────────────────────────────────────
out_csv = OUT / "CCC_consensus_CellChat_CellPhoneDB.csv"
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "sender", "receiver", "ligand", "receptor", "cellchat_score",
        "cpdb_pval", "pathway", "nrf2_relevant", "ferroptosis_relevant"
    ])
    writer.writeheader()
    for c in consensus:
        c["nrf2_relevant"] = any(g.upper() in c["ligand"].upper() for g in nrf2_ligands)
        c["ferroptosis_relevant"] = any(g.upper() in c["ligand"].upper() for g in ferr_ligands)
        writer.writerow(c)
print(f"\nSaved: {out_csv}")

# Summary stats
print(f"\n=== SUMMARY ===")
print(f"CellChat myeloid->CM/Peri: {len(cc_mye_target)}")
print(f"CellPhoneDB myeloid->CM/Peri: {len(cpdb_mye_target)}")
print(f"Consensus: {len(consensus)}")
print(f"Unique pathways: {len(pathways)}")
