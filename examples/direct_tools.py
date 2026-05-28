"""Direct Tool Usage Example.

Demonstrates how to use CellAgent's tools directly without the
agent loop, for users who want more control over the analysis.
"""

import os
import sys
import numpy as np
import anndata as ad
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cellagent.tools.preprocessing import (
    calculate_qc_metrics,
    filter_cells_and_genes,
    normalize_and_log_transform,
    select_highly_variable_genes,
    run_pca,
)
from cellagent.tools.clustering import (
    compute_neighbors,
    run_umap,
    run_leiden_clustering,
    evaluate_clustering_quality,
)
from cellagent.tools.annotation import (
    find_marker_genes,
    annotate_with_markers,
    annotate_with_llm,
)
from cellagent.tools.differential import (
    differential_expression,
    volcano_plot,
)
from cellagent.tools.visualization import (
    plot_umap_colored,
    plot_cell_composition,
    generate_analysis_summary,
)


def main():
    """Run a step-by-step analysis using tools directly."""
    output_dir = "./output/direct_tools"
    os.makedirs(output_dir, exist_ok=True)

    # ============================================================
    # Step 1: Load Data
    # ============================================================
    print("Step 1: Loading PBMC3k dataset...")
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # ============================================================
    # Step 2: Quality Control
    # ============================================================
    print("\nStep 2: Quality Control...")
    qc_report = calculate_qc_metrics(adata, output_dir=output_dir)
    print(qc_report)

    # ============================================================
    # Step 3: Filtering
    # ============================================================
    print("\nStep 3: Filtering...")
    filter_report = filter_cells_and_genes(
        adata, min_genes=200, max_genes=2500, max_pct_mt=5
    )
    print(filter_report)

    # ============================================================
    # Step 4: Normalization
    # ============================================================
    print("\nStep 4: Normalization...")
    norm_report = normalize_and_log_transform(adata)
    print(norm_report)

    # ============================================================
    # Step 5: Feature Selection
    # ============================================================
    print("\nStep 5: HVG Selection...")
    hvg_report = select_highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    print(hvg_report)

    # Subset to HVGs for PCA
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    # ============================================================
    # Step 6: PCA
    # ============================================================
    print("\nStep 6: PCA...")
    pca_report = run_pca(adata_hvg, n_comps=50, output_dir=output_dir)
    print(pca_report)

    # ============================================================
    # Step 7: Neighbors and Clustering
    # ============================================================
    print("\nStep 7: Clustering...")
    neighbors_report = compute_neighbors(adata_hvg, n_neighbors=10, n_pcs=40)
    print(neighbors_report)

    cluster_report = run_leiden_clustering(adata_hvg, resolution=0.5)
    print(cluster_report)

    quality_report = evaluate_clustering_quality(adata_hvg)
    print(quality_report)

    # ============================================================
    # Step 8: UMAP
    # ============================================================
    print("\nStep 8: UMAP...")
    umap_report = run_umap(adata_hvg, min_dist=0.3, output_dir=output_dir)
    print(umap_report)

    # ============================================================
    # Step 9: Marker Genes
    # ============================================================
    print("\nStep 9: Finding Marker Genes...")
    marker_report = find_marker_genes(adata_hvg, output_dir=output_dir)
    print(marker_report)

    # ============================================================
    # Step 10: Cell Type Annotation
    # ============================================================
    print("\nStep 10: Cell Type Annotation...")

    # Known PBMC markers
    pbmc_markers = {
        "CD4+ T cells": ["CD3D", "CD3E", "IL7R"],
        "CD8+ T cells": ["CD3D", "CD3E", "CD8A", "CD8B"],
        "NK cells": ["NKG7", "GNLY", "KLRD1"],
        "B cells": ["CD79A", "MS4A1", "CD79B"],
        "CD14+ Monocytes": ["CD14", "LYZ", "S100A9"],
        "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7", "LST1"],
        "Dendritic cells": ["FCER1A", "CD1C"],
        "Platelets": ["PPBP", "PF4"],
    }

    annotation_report = annotate_with_markers(
        adata_hvg, pbmc_markers, output_dir=output_dir
    )
    print(annotation_report)

    # ============================================================
    # Step 11: Visualization
    # ============================================================
    print("\nStep 11: Visualization...")
    plot_report = plot_umap_colored(
        adata_hvg,
        color=["leiden", "cell_type_marker"],
        output_dir=output_dir,
        filename="umap_annotated.png",
    )
    print(plot_report)

    comp_report = plot_cell_composition(
        adata_hvg, "cell_type_marker", output_dir=output_dir
    )
    print(comp_report)

    # ============================================================
    # Step 12: Summary Report
    # ============================================================
    print("\nStep 12: Generating Summary...")
    summary = generate_analysis_summary(
        adata_hvg, output_dir=output_dir, title="PBMC3k Analysis Report"
    )
    print(summary)

    print(f"\nAll outputs saved to: {output_dir}/")
    print(f"Files: {os.listdir(output_dir)}")


if __name__ == "__main__":
    main()
