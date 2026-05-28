# Multi-Omics Integration Guide

## Metadata

**Short_Description**: Guide for integrating multiple single-cell omics modalities (RNA, ATAC, protein) and multi-sample batch correction strategies.
**Authors**: CellAgent Knowledge Base
**Version**: 1.0
**Last_Updated**: 2025-01

---

## Overview

Modern single-cell experiments increasingly measure multiple modalities simultaneously (e.g., CITE-seq for RNA+protein, 10X Multiome for RNA+ATAC). This guide covers strategies for integrating these data types.

## Multi-Sample Integration

### When to Integrate
- Multiple samples from different patients/conditions
- Technical replicates with batch effects
- Data from different experiments or labs

### Method Selection Guide

| Method | Speed | Batch Effect Strength | Preserves Biology | Scalability |
|--------|-------|----------------------|-------------------|-------------|
| Harmony | Fast | Moderate | Good | Excellent |
| BBKNN | Fast | Mild-Moderate | Very Good | Good |
| scVI | Slow | Strong | Good | Good |
| Scanorama | Moderate | Moderate | Good | Good |
| Combat | Fast | Mild | Good | Excellent |

### Recommended Workflow
1. Process each sample independently through QC
2. Merge samples: `adata = ad.concat([adata1, adata2], label='sample')`
3. Run HVG selection with batch_key
4. Run PCA on merged data
5. Apply integration method
6. Evaluate integration quality

## CITE-seq (RNA + Protein)

### Processing Protein Data
1. CLR normalization (centered log-ratio) for ADT counts
2. DSB normalization (preferred, uses empty droplets as background)
3. Store in `adata.obsm['protein']`

### Joint Analysis
- Use WNN (Weighted Nearest Neighbors) from Seurat v4
- Or concatenate normalized RNA and protein features
- Protein data is especially useful for surface marker-based cell typing

## 10X Multiome (RNA + ATAC)

### ATAC Processing
1. Peak calling and count matrix generation
2. TF-IDF normalization
3. LSI (Latent Semantic Indexing) for dimensionality reduction

### Joint Analysis
- Use MOFA+ for multi-omics factor analysis
- Or use ArchR/Signac for integrated analysis
- Gene activity scores from ATAC can validate RNA expression

## Evaluation Metrics

### Batch Mixing
- **Batch silhouette score**: Should be close to 0 (well-mixed)
- **kBET**: Measures local batch mixing
- **Batch entropy**: Higher = better mixing

### Biological Conservation
- **Cell type silhouette score**: Should be high (well-separated)
- **ARI/NMI**: Compare clustering with known cell types
- **Marker gene expression**: Should be consistent after integration
