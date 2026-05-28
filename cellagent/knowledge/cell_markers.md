# Cell Type Marker Gene Reference

## Metadata

**Short_Description**: Comprehensive reference of canonical marker genes for major cell types across human and mouse tissues, organized by tissue type.
**Authors**: CellAgent Knowledge Base
**Version**: 1.0
**Last_Updated**: 2025-01

---

## Overview

This document provides curated marker gene lists for cell type annotation in single-cell RNA-seq analysis. Markers are organized by tissue type and validated across multiple studies.

## Human PBMC (Peripheral Blood Mononuclear Cells)

| Cell Type | Primary Markers | Secondary Markers |
|-----------|----------------|-------------------|
| CD4+ T cells | CD3D, CD3E, CD4 | IL7R, CCR7 (naive), FOXP3 (Treg) |
| CD8+ T cells | CD3D, CD3E, CD8A, CD8B | GZMB, PRF1 (cytotoxic) |
| Naive T cells | CCR7, LEF1, TCF7, SELL | CD27, CD28 |
| Memory T cells | S100A4, IL7R, GZMK | CD44 |
| Regulatory T cells | FOXP3, IL2RA, CTLA4 | TIGIT, IKZF2 |
| NK cells | NKG7, GNLY, KLRD1 | NCAM1 (CD56), KLRF1, GZMB |
| B cells | CD19, MS4A1 (CD20), CD79A | CD79B, PAX5 |
| Plasma cells | MZB1, SDC1 (CD138), JCHAIN | XBP1, IGHG1 |
| CD14+ Monocytes | CD14, LYZ, S100A9 | S100A8, VCAN |
| CD16+ Monocytes | FCGR3A (CD16), LST1, MS4A7 | LILRB2 |
| Dendritic cells (cDC) | FCER1A, CD1C | CLEC10A, HLA-DQA1 |
| Plasmacytoid DC | LILRA4, IRF7, CLEC4C | TCF4, IL3RA |
| Platelets | PPBP, PF4, GP9 | ITGA2B |
| Erythrocytes | HBB, HBA1, HBA2 | ALAS2 |

## Human Brain

| Cell Type | Primary Markers | Secondary Markers |
|-----------|----------------|-------------------|
| Excitatory neurons | SLC17A7, SATB2, CAMK2A | NRGN, SYT1 |
| Inhibitory neurons | GAD1, GAD2, SLC32A1 | PVALB, SST, VIP |
| Astrocytes | GFAP, AQP4, SLC1A2 | ALDH1L1, S100B |
| Oligodendrocytes | MBP, MOG, PLP1 | OLIG1, OLIG2 (OPC) |
| Microglia | CX3CR1, P2RY12, TMEM119 | AIF1, CSF1R |
| Endothelial cells | CLDN5, FLT1, PECAM1 | VWF |
| OPC | PDGFRA, CSPG4, OLIG2 | SOX10 |

## Human Lung

| Cell Type | Primary Markers | Secondary Markers |
|-----------|----------------|-------------------|
| AT1 cells | AGER, PDPN, AQP5 | HOPX |
| AT2 cells | SFTPC, SFTPA1, SFTPB | ABCA3, NKX2-1 |
| Ciliated cells | FOXJ1, PIFO, TPPP3 | DNAH5 |
| Club cells | SCGB1A1, SCGB3A1 | CYP2F1 |
| Basal cells | KRT5, KRT17, TP63 | NGFR |
| Fibroblasts | COL1A1, DCN, LUM | VIM, COL3A1 |
| Endothelial cells | PECAM1, VWF, CDH5 | ERG |
| Smooth muscle | ACTA2, MYH11, TAGLN | DES |

## Human Liver

| Cell Type | Primary Markers | Secondary Markers |
|-----------|----------------|-------------------|
| Hepatocytes | ALB, APOA1, HP | CYP3A4, TF |
| Cholangiocytes | KRT19, KRT7, EPCAM | SOX9 |
| Stellate cells | RGS5, ACTA2, DES | LRAT |
| Kupffer cells | CD68, MARCO, CD163 | CLEC4F |
| LSEC | CLEC4G, FCGR2B, STAB2 | LYVE1 |

## Mouse Brain

| Cell Type | Primary Markers | Secondary Markers |
|-----------|----------------|-------------------|
| Excitatory neurons | Slc17a7, Satb2 | Neurod6 |
| Inhibitory neurons | Gad1, Gad2 | Pvalb, Sst, Vip |
| Astrocytes | Gfap, Aqp4, Slc1a3 | Aldh1l1 |
| Oligodendrocytes | Mbp, Mog, Plp1 | Mag |
| Microglia | Cx3cr1, P2ry12, Tmem119 | Aif1 |

## Mouse Immune

| Cell Type | Primary Markers | Secondary Markers |
|-----------|----------------|-------------------|
| T cells | Cd3d, Cd3e | Cd4, Cd8a |
| B cells | Cd19, Ms4a1, Cd79a | Pax5 |
| NK cells | Nkg7, Klrb1c | Ncr1 |
| Macrophages | Adgre1 (F4/80), Cd68 | Csf1r |
| Neutrophils | S100a8, S100a9, Ly6g | Cxcr2 |
| Dendritic cells | Itgax (Cd11c), H2-Ab1 | Flt3 |

## Annotation Best Practices

1. **Always check gene presence**: Verify marker genes exist in your dataset before scoring
2. **Use multiple markers**: Never rely on a single gene; use 3-5 markers per cell type
3. **Consider tissue context**: The same gene may mark different cell types in different tissues
4. **Validate with dot plots**: Visualize marker expression across clusters
5. **Cross-reference databases**: Use CellMarker, PanglaoDB, or CellTypist for additional markers
6. **Check for batch effects**: Ensure marker expression is consistent across batches
7. **Hierarchical annotation**: Start with broad categories, then refine subtypes
