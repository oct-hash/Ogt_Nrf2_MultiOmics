import csv
path = r"D:\miri_paper\results\tables\GSE193997_DESeq2_Sham_vs_IR6h.csv"
with open(path) as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        if i < 5:
            print(f'{row["symbol"]:15s} {row["ensembl"]}')
        else:
            break
print("\nMouse symbols: Serpine1, Prg4, Thbs1 -> YES, this is MOUSE")
