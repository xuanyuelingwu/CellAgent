# Troubleshooting Guide for scRNA-seq Analysis

## Metadata

**Short_Description**: Solutions for common problems encountered during single-cell RNA-seq analysis, including clustering issues, annotation failures, and technical artifacts.
**Authors**: CellAgent Knowledge Base
**Version**: 1.0
**Last_Updated**: 2025-01

---

## Overview

This guide covers common issues and their solutions in single-cell RNA-seq analysis.

## Clustering Issues

### Problem: Too many or too few clusters
**Symptoms**: Clusters don't correspond to known biology
**Solutions**:
- Adjust Leiden resolution parameter (lower for fewer clusters, higher for more)
- Try multiple resolutions and compare with known markers
- Check if batch effects are driving clustering

### Problem: Clusters driven by technical artifacts
**Symptoms**: Clusters correlate with QC metrics (n_genes, total_counts, pct_mt)
**Solutions**:
- Improve QC filtering
- Regress out confounding variables: `sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])`
- Use more stringent filtering thresholds

### Problem: One large cluster with everything mixed
**Symptoms**: Most cells in one cluster, poor separation
**Solutions**:
- Increase resolution parameter
- Check normalization (was log1p applied?)
- Ensure HVG selection was performed
- Try different number of PCs

## Annotation Issues

### Problem: No clear marker genes for a cluster
**Symptoms**: Top DEGs are housekeeping genes or unknown
**Solutions**:
- Check if the cluster is a doublet population
- Look at broader marker categories
- Consider the cluster may be a transitional state
- Use LLM-assisted annotation for interpretation

### Problem: Same cell type in multiple clusters
**Symptoms**: Multiple clusters express the same markers
**Solutions**:
- Lower clustering resolution
- The subclusters may represent valid subtypes
- Check for batch effects splitting cell types
- Merge clusters with similar marker profiles

## Technical Issues

### Problem: Memory errors with large datasets
**Solutions**:
- Use sparse matrices (CSR format)
- Subsample for initial exploration
- Use backed mode: `sc.read_h5ad(file, backed='r')`
- Reduce n_top_genes for HVG selection

### Problem: Slow computation
**Solutions**:
- Use approximate methods for neighbors (default in scanpy)
- Reduce n_pcs for neighbor computation
- Subsample for parameter optimization
- Use GPU-accelerated tools (rapids-singlecell)

### Problem: Gene names not matching
**Symptoms**: Known markers not found in data
**Solutions**:
- Check if gene names are symbols vs. Ensembl IDs
- Convert between naming conventions
- Check organism (human vs. mouse capitalization)
- Look for aliases: `adata.var[adata.var_names.str.contains('CD3')]`

## Integration Issues

### Problem: Over-correction removes biological signal
**Symptoms**: Known distinct cell types merge after integration
**Solutions**:
- Use less aggressive integration (lower Harmony sigma)
- Try BBKNN instead of Harmony
- Evaluate with cell type silhouette scores
- Compare integrated vs. unintegrated clustering

### Problem: Batch effects persist after integration
**Symptoms**: Clusters still separate by batch
**Solutions**:
- Try more aggressive integration parameters
- Use scVI for complex batch effects
- Check if biological differences exist between batches
- Ensure consistent preprocessing across batches

## Visualization Issues

### Problem: UMAP looks fragmented
**Solutions**:
- Increase n_neighbors (try 30-50)
- Increase min_dist (try 0.3-0.5)
- Use PAGA-initialized UMAP
- Check if data is properly normalized

### Problem: UMAP looks like a blob
**Solutions**:
- Decrease min_dist (try 0.1)
- Increase n_neighbors
- Check if HVGs were selected
- Verify PCA was computed on HVGs only
