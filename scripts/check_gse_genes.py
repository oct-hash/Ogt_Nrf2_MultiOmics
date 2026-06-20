import gzip
lines = []
with gzip.open(r"D:\Cardiac_Ogt_MultiOmics\validation_data\GSE57338_series_matrix.txt.gz", "rt", errors="ignore") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    ls = line.strip()
    if ls == "!platform_table_begin":
        # Print header + first 15 data rows
        hdr = lines[i+1].strip().split("\t")
        print(f"Platform header ({len(hdr)} cols): {hdr[:4]}")
        for j in range(i+2, min(i+17, len(lines))):
            cols = lines[j].strip().split("\t")
            if cols[0].strip('"').startswith("!"):
                break
            gene = cols[1].strip('"') if len(cols) > 1 else "NA"
            print(f"  {cols[0][:25]:25s} -> {gene[:30]}")
        break
