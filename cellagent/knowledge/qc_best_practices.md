# Quality Control Best Practices for scRNA-seq

## Metadata

**Short_Description**: Comprehensive guide for quality control of single-cell RNA-seq data, including filtering thresholds, doublet detection, and ambient RNA removal.
**Authors**: CellAgent Knowledge Base
**Version**: 1.0
**Last_Updated**: 2025-01

---

## Overview

Quality control (QC) is the critical first step in any single-cell RNA-seq analysis. Poor QC leads to unreliable downstream results. This document covers best practices for assessing and filtering single-cell data quality.

## Key QC Metrics

### Per-Cell Metrics
1. **Number of genes detected (n_genes_by_counts)**: Cells with very few genes may be empty droplets or dead cells. Cells with extremely high gene counts may be doublets.
2. **Total UMI counts (total_counts)**: Similar to gene count, extreme values indicate quality issues.
3. **Mitochondrial gene percentage (pct_counts_mt)**: High mitochondrial content (>20%) indicates stressed or dying cells where cytoplasmic mRNA has leaked out.
4. **Ribosomal gene percentage (pct_counts_ribo)**: Extremely high ribosomal content may indicate specific cell states.

### Per-Gene Metrics
1. **Number of cells expressing the gene**: Genes detected in very few cells provide little information and add noise.

## Recommended Filtering Thresholds

| Metric | Typical Range | Notes |
|--------|--------------|-------|
| min_genes | 200-500 | Lower for droplet-based, higher for plate-based |
| max_genes | 2500-8000 | Depends on cell type complexity |
| min_cells | 3-10 | Remove very rare genes |
| max_pct_mt | 5-20% | Tissue-dependent; brain tissue tolerates higher |
| max_counts | Dataset-specific | Use MAD-based outlier detection |

## Tissue-Specific Considerations

- **PBMC**: max_pct_mt = 10-15%, typically 200-5000 genes/cell
- **Brain**: max_pct_mt = 15-25%, neurons have higher MT content
- **Tumor**: Higher heterogeneity, use adaptive thresholds
- **Embryonic**: Rapidly dividing cells may have different profiles

## Doublet Detection

Doublets occur when two cells are captured in one droplet (~4-8% for 10X Chromium).

**Recommended approach**:
1. Run Scrublet or DoubletFinder before filtering
2. Expected doublet rate: ~0.8% per 1000 cells captured
3. For 10,000 cells: expect ~8% doublets
4. Remove predicted doublets before downstream analysis

## Ambient RNA Removal

Cell-free RNA from lysed cells contaminates all droplets.

**Tools**: SoupX, CellBender, DecontX
**When to use**: When you see unexpected expression of cell-type-specific genes across all clusters

## Adaptive Thresholds (MAD-based)

Instead of fixed thresholds, use Median Absolute Deviation (MAD):
```python
import numpy as np

def mad_outlier(data, nmads=5):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    lower = median - nmads * mad
    upper = median + nmads * mad
    return lower, upper
```

Use 5 MADs for a lenient filter, 3 MADs for stringent.

## QC Workflow

1. Calculate QC metrics (before any filtering)
2. Visualize distributions (violin plots, scatter plots)
3. Identify outliers using adaptive or fixed thresholds
4. Detect doublets
5. Remove ambient RNA (if needed)
6. Filter cells and genes
7. Re-visualize to confirm filtering quality
