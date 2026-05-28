# Single-Cell RNA-seq Analysis Workflow Guide

## Metadata

**Short_Description**: Step-by-step guide for standard scRNA-seq analysis workflows, from raw data to biological interpretation, with decision points and parameter recommendations.
**Authors**: CellAgent Knowledge Base
**Version**: 1.0
**Last_Updated**: 2025-01

---

## Overview

This document provides a comprehensive workflow guide for single-cell RNA-seq data analysis, covering the standard pipeline from data loading to biological interpretation.

## Standard Analysis Pipeline

### Phase 1: Data Loading and QC

**Steps:**
1. Load data (h5ad, 10X MTX, or other formats)
2. Calculate QC metrics (genes/cell, UMI/cell, %MT, %ribo)
3. Visualize QC distributions
4. Filter low-quality cells and genes
5. Detect and remove doublets (optional but recommended)

**Decision Points:**
- Fixed vs. adaptive (MAD-based) thresholds
- Whether to remove doublets or just flag them
- Whether ambient RNA correction is needed

### Phase 2: Normalization and Feature Selection

**Steps:**
1. Normalize total counts per cell (target_sum=1e4)
2. Log-transform (log1p)
3. Select highly variable genes (HVGs)
4. Save raw counts for later use

**Key Parameters:**
- n_top_genes: 2000-5000 (2000 is standard)
- HVG method: seurat_v3 (recommended), seurat, cell_ranger

### Phase 3: Dimensionality Reduction

**Steps:**
1. Scale data (optional, max_value=10)
2. Run PCA (50 components)
3. Examine elbow plot to determine n_pcs
4. Compute neighbor graph

**Key Parameters:**
- n_pcs: 20-50 (use elbow plot to decide)
- n_neighbors: 10-30 (15 is standard)

### Phase 4: Clustering and Visualization

**Steps:**
1. Run Leiden clustering (recommended over Louvain)
2. Compute UMAP embedding
3. Visualize clusters on UMAP
4. Evaluate clustering quality

**Key Parameters:**
- resolution: 0.1-2.0 (0.5 is a good starting point)
  - Lower resolution = fewer, larger clusters
  - Higher resolution = more, smaller clusters
- UMAP min_dist: 0.1-0.5

### Phase 5: Cell Type Annotation

**Approaches (in order of recommendation):**
1. **Marker-based**: Use known marker genes with dot plots
2. **LLM-assisted**: Use AI to interpret marker genes
3. **Reference-based**: Transfer labels from annotated reference
4. **Automated**: Use CellTypist or scType

**Workflow:**
1. Find marker genes per cluster (rank_genes_groups)
2. Compare markers with known cell type signatures
3. Assign cell type labels
4. Validate with visualization

### Phase 6: Differential Expression

**Steps:**
1. Define comparison groups
2. Run DE test (Wilcoxon rank-sum recommended)
3. Filter by log2FC and adjusted p-value
4. Generate volcano plots
5. Perform pathway enrichment analysis

### Phase 7: Advanced Analysis (Optional)

**Trajectory Analysis:**
- PAGA for trajectory topology
- Diffusion pseudotime for ordering
- RNA velocity for directionality

**Integration:**
- Harmony for fast batch correction
- BBKNN for neighbor-based integration
- scVI for deep learning-based integration

**Cell Communication:**
- CellChat or LIANA for ligand-receptor analysis

## Common Pitfalls

1. **Over-filtering**: Removing too many cells/genes loses biological signal
2. **Wrong resolution**: Too high creates artificial subclusters, too low merges distinct types
3. **Ignoring batch effects**: Multi-sample data needs integration
4. **Single-marker annotation**: Always use multiple markers
5. **Circular analysis**: Don't use the same data for discovery and validation
