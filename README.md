# Ogt_Nrf2_MultiOmics

**Multi-omics dissection of the OGT-Nrf2-ferroptosis regulatory axis in acute myocardial ischemia-reperfusion injury**

## Overview

This repository contains the complete computational analysis pipeline for a multi-omics study investigating the O-GlcNAc transferase (OGT)-Nrf2-ferroptosis regulatory axis in acute myocardial ischemia-reperfusion (I/R) injury. The study integrates:

- **Bulk RNA-seq**: GSE193997 (discovery, n=24), GSE307385 (validation, n=15)
- **Single-nucleus RNA-seq**: GSE240848 (48,633 nuclei, 9 cell types)
- **Proteomics**: PXD046631 (LFQ), PXD040135 (LFQ), PXD035154 (TMT)
- **External validation**: GTEx v8 (healthy human LV, n=431), GSE57338 (chronic HF, 177 HF vs 136 NF)

## Key Findings

1. Ogt is significantly downregulated (log2FC = −0.56, padj = 3.09×10⁻⁷) at 6h post-I/R
2. Coordinated Nrf2-ARE and ferroptosis pathway activation (NES = 1.77, FDR = 0.001)
3. mRNA-protein discordance (R = 0.253) reveals predominant post-transcriptional regulation
4. Single-nucleus analysis uncovers cell-type dichotomy: Ogt in cardiomyocytes/pericytes, Nrf2 effectors in myeloid cells
5. SPP1-CD44 identified as the dominant myeloid-to-parenchymal paracrine signal
6. RELA/NF-kB is the sole transcription factor linking OGT substrate network to Nrf2 activation
7. Disease-stage specificity confirmed: Nrf2-ferroptosis signature is acute I/R-specific

## Repository Structure

```
├── scripts/                 # All analysis scripts
│   ├── figures.py           # Main figure generation
│   ├── figures_final.py     # Final multi-panel assembly
│   ├── generate_supplementary_figures.py  # Figures S1-S12
│   ├── create_supplementary_tables.py     # Tables S1-S15 (Excel)
│   ├── step1_ccc_dual.py    # Cell-cell communication (CellChat + CellPhoneDB)
│   ├── step2_residual.py    # mRNA-protein residual analysis
│   ├── step3_ogt_enrichment.py  # OGT-PIN enrichment analysis
│   ├── step4_external_validation.py  # GSE57338 validation
│   ├── pipeline.py          # Master analysis pipeline
│   └── ...
├── output_tables/           # Analysis outputs (CSV)
│   ├── mRNA_protein_residuals.csv
│   ├── CCC_consensus_CellChat_CellPhoneDB.csv
│   ├── step3_ogt_enrichment.csv
│   ├── GSE57338_direction_consistency.csv
│   └── ...
├── figures/                 # Generated figures (300 dpi PNG)
├── manuscript/              # Manuscript text files
│   ├── _results.txt
│   ├── _methods.txt
│   ├── _intro_disc.txt
│   ├── _references.txt
│   ├── _figure_legends.txt
│   └── _author_info.txt
└── Supplementary_Tables_S1_S15.xlsx  # Combined supplementary tables
```

## Data Availability

| Dataset | Accession | Type | Species |
|---------|-----------|------|---------|
| Discovery bulk RNA-seq | GSE193997 | RNA-seq | Mouse |
| Validation bulk RNA-seq | GSE307385 | RNA-seq | Mouse |
| snRNA-seq | GSE240848 | 10x Genomics | Rat |
| Chronic HF validation | GSE57338 | Microarray | Human |
| Proteomics LFQ | PXD046631 | Mass Spec | Mouse |
| Proteomics LFQ | PXD040135 | Mass Spec | Mouse |
| Proteomics TMT | PXD035154 | Mass Spec | Mouse |
| Healthy baseline | GTEx v8 | RNA-seq | Human |

## Requirements

- Python 3.10+
- R 4.2+
- Key Python packages: pandas, numpy, matplotlib, seaborn, scanpy, openpyxl
- Key R packages: DESeq2, fgsea, WGCNA, Seurat, CellChat

## Author

Guoxin Zhang  
Quanzhou First Hospital Affiliated to Fujian Medical University  
Quanzhou, Fujian 362000, China  
Email: 9201941391@fjmu.edu.cn

## License

This project is licensed under the MIT License.
