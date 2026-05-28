"""Visualization tools for single-cell data.

Provides publication-quality plotting functions for various aspects
of single-cell analysis results.
"""

TOOL_DESCRIPTIONS = [
    {
        "name": "plot_umap_colored",
        "description": (
            "Generate UMAP plots colored by one or more variables. "
            "Supports coloring by cluster labels, cell types, gene expression, "
            "QC metrics, or any other variable in adata.obs."
        ),
        "module": "cellagent.tools.visualization",
        "category": "visualization",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with UMAP computed"},
            {"name": "color", "type": "list",
             "description": "List of variables to color by (obs columns or gene names)"},
        ],
        "optional_parameters": [
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save plots"},
            {"name": "filename", "type": "str", "default": "umap_colored.png",
             "description": "Output filename"},
            {"name": "ncols", "type": "int", "default": 3,
             "description": "Number of columns in multi-panel plot"},
        ],
        "returns": "str",
        "example": 'result = plot_umap_colored(adata, color=["leiden", "cell_type", "CD3D"])',
    },
    {
        "name": "plot_dotplot",
        "description": (
            "Generate a dot plot showing expression of marker genes across clusters. "
            "Dot size represents fraction of cells expressing the gene, "
            "color represents mean expression level."
        ),
        "module": "cellagent.tools.visualization",
        "category": "visualization",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "marker_genes", "type": "dict",
             "description": "Dictionary mapping categories to gene lists"},
        ],
        "optional_parameters": [
            {"name": "groupby", "type": "str", "default": "leiden",
             "description": "Grouping variable"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save plot"},
        ],
        "returns": "str",
        "example": 'result = plot_dotplot(adata, {"T cells": ["CD3D"], "B cells": ["CD19"]})',
    },
    {
        "name": "plot_heatmap",
        "description": (
            "Generate a heatmap of gene expression across clusters or cell types. "
            "Shows expression patterns of selected genes across groups."
        ),
        "module": "cellagent.tools.visualization",
        "category": "visualization",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "genes", "type": "list", "description": "List of genes to plot"},
        ],
        "optional_parameters": [
            {"name": "groupby", "type": "str", "default": "leiden",
             "description": "Grouping variable"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save plot"},
        ],
        "returns": "str",
        "example": 'result = plot_heatmap(adata, ["CD3D", "CD19", "CD14"], groupby="cell_type")',
    },
    {
        "name": "plot_stacked_violin",
        "description": (
            "Generate stacked violin plots for gene expression across groups. "
            "Compact visualization showing expression distributions."
        ),
        "module": "cellagent.tools.visualization",
        "category": "visualization",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "genes", "type": "list", "description": "List of genes to plot"},
        ],
        "optional_parameters": [
            {"name": "groupby", "type": "str", "default": "leiden",
             "description": "Grouping variable"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save plot"},
        ],
        "returns": "str",
        "example": 'result = plot_stacked_violin(adata, ["CD3D", "CD19"], groupby="leiden")',
    },
    {
        "name": "plot_cell_composition",
        "description": (
            "Generate cell type composition bar plots showing the proportion "
            "of each cell type across samples, conditions, or batches."
        ),
        "module": "cellagent.tools.visualization",
        "category": "visualization",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "cell_type_key", "type": "str",
             "description": "Key in adata.obs for cell type labels"},
        ],
        "optional_parameters": [
            {"name": "condition_key", "type": "str", "default": None,
             "description": "Key for condition/sample grouping"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save plot"},
        ],
        "returns": "str",
        "example": 'result = plot_cell_composition(adata, "cell_type", condition_key="sample")',
    },
    {
        "name": "generate_analysis_summary",
        "description": (
            "Generate a comprehensive summary report of the single-cell analysis, "
            "including dataset statistics, QC metrics, clustering results, "
            "cell type annotations, and key findings. Outputs an HTML report."
        ),
        "module": "cellagent.tools.visualization",
        "category": "visualization",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with analysis results"},
        ],
        "optional_parameters": [
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save report"},
            {"name": "title", "type": "str", "default": "Single-Cell Analysis Report",
             "description": "Report title"},
        ],
        "returns": "str",
        "example": 'result = generate_analysis_summary(adata, title="PBMC Analysis")',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def plot_umap_colored(adata, color: list, output_dir: str = "./output",
                      filename: str = "umap_colored.png", ncols: int = 3) -> str:
    """Generate colored UMAP plots."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if "X_umap" not in adata.obsm:
        return "Error: UMAP embedding not found. Run run_umap before plot_umap_colored."

    # Filter to valid color keys
    valid_colors = []
    for c in color:
        if c in adata.obs.columns or c in adata.var_names:
            valid_colors.append(c)

    if not valid_colors:
        return f"Error: None of the specified variables found: {color}"

    sc.pl.umap(adata, color=valid_colors, ncols=ncols, show=False)
    path = os.path.join(output_dir, filename)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    report = (
        f"UMAP Plot Generated:\n"
        f"  Variables: {valid_colors}\n"
        f"  Saved to: {path}\n"
    )
    return report


def plot_dotplot(adata, marker_genes: dict, groupby: str = "leiden",
                 output_dir: str = "./output") -> str:
    """Generate dot plot."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if groupby not in adata.obs:
        return f"Error: Group key '{groupby}' not found in adata.obs."

    # Filter to valid genes
    valid_markers = {}
    for ct, genes in marker_genes.items():
        valid = [g for g in genes if g in adata.var_names]
        if valid:
            valid_markers[ct] = valid

    if not valid_markers:
        return "Error: No valid marker genes found in the dataset."

    sc.pl.dotplot(adata, valid_markers, groupby=groupby, show=False)
    path = os.path.join(output_dir, "dotplot.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    report = f"Dot Plot Generated:\n  Groupby: {groupby}\n  Saved to: {path}\n"
    return report


def plot_heatmap(adata, genes: list, groupby: str = "leiden",
                 output_dir: str = "./output") -> str:
    """Generate heatmap."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if groupby not in adata.obs:
        return f"Error: Group key '{groupby}' not found in adata.obs."

    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        return "Error: No valid genes found."

    sc.pl.heatmap(adata, valid_genes, groupby=groupby, show=False)
    path = os.path.join(output_dir, "heatmap.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    report = f"Heatmap Generated:\n  Genes: {len(valid_genes)}\n  Saved to: {path}\n"
    return report


def plot_stacked_violin(adata, genes: list, groupby: str = "leiden",
                        output_dir: str = "./output") -> str:
    """Generate stacked violin plot."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if groupby not in adata.obs:
        return f"Error: Group key '{groupby}' not found in adata.obs."

    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        return "Error: No valid genes found."

    sc.pl.stacked_violin(adata, valid_genes, groupby=groupby, show=False)
    path = os.path.join(output_dir, "stacked_violin.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    report = f"Stacked Violin Plot Generated:\n  Genes: {valid_genes}\n  Saved to: {path}\n"
    return report


def plot_cell_composition(adata, cell_type_key: str, condition_key: str = None,
                          output_dir: str = "./output") -> str:
    """Plot cell type composition."""
    import pandas as pd
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if cell_type_key not in adata.obs:
        return f"Error: Cell type key '{cell_type_key}' not found in adata.obs."
    if condition_key and condition_key not in adata.obs:
        return f"Error: Condition key '{condition_key}' not found in adata.obs."

    if condition_key and condition_key in adata.obs.columns:
        ct = pd.crosstab(adata.obs[condition_key], adata.obs[cell_type_key], normalize="index")
        ax = ct.plot(kind="bar", stacked=True, figsize=(12, 6))
        ax.set_ylabel("Proportion")
        ax.set_title("Cell Type Composition by Condition")
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        counts = adata.obs[cell_type_key].value_counts()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        counts.plot(kind="bar", ax=ax1)
        ax1.set_title("Cell Type Counts")
        ax1.set_ylabel("Number of Cells")
        counts.plot(kind="pie", autopct="%1.1f%%", ax=ax2)
        ax2.set_title("Cell Type Proportions")

    plt.tight_layout()
    path = os.path.join(output_dir, "cell_composition.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    report = f"Cell Composition Plot Generated:\n"
    report += f"  Cell type key: {cell_type_key}\n"
    if condition_key:
        report += f"  Condition key: {condition_key}\n"
    report += f"  Saved to: {path}\n"
    return report


def generate_analysis_summary(adata, output_dir: str = "./output",
                               title: str = "Single-Cell Analysis Report") -> str:
    """Generate comprehensive analysis summary."""
    import os

    os.makedirs(output_dir, exist_ok=True)

    lines = [f"# {title}\n"]
    lines.append("## Dataset Overview\n")
    lines.append(f"- **Cells**: {adata.n_obs}")
    lines.append(f"- **Genes**: {adata.n_vars}")
    lines.append(f"- **Obs columns**: {', '.join(adata.obs.columns[:20])}")
    lines.append(f"- **Embeddings**: {', '.join(adata.obsm.keys())}")

    if "pct_counts_mt" in adata.obs.columns:
        lines.append("\n## QC Metrics\n")
        lines.append(f"- Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
        lines.append(f"- Median UMI/cell: {adata.obs['total_counts'].median():.0f}")
        lines.append(f"- Median %MT: {adata.obs['pct_counts_mt'].median():.2f}%")

    for key in ["leiden", "louvain"]:
        if key in adata.obs.columns:
            lines.append(f"\n## Clustering ({key})\n")
            lines.append(f"- Number of clusters: {adata.obs[key].nunique()}")
            for cluster in sorted(adata.obs[key].unique(), key=lambda x: int(x) if str(x).isdigit() else x):
                n = (adata.obs[key] == cluster).sum()
                lines.append(f"- Cluster {cluster}: {n} cells ({n/adata.n_obs*100:.1f}%)")

    for key in ["cell_type_marker", "cell_type_llm"]:
        if key in adata.obs.columns:
            lines.append(f"\n## Cell Type Annotation ({key})\n")
            for ct in adata.obs[key].value_counts().index:
                n = (adata.obs[key] == ct).sum()
                lines.append(f"- {ct}: {n} cells ({n/adata.n_obs*100:.1f}%)")

    if "dpt_pseudotime" in adata.obs.columns:
        lines.append("\n## Trajectory Analysis\n")
        lines.append(f"- Pseudotime range: [{adata.obs['dpt_pseudotime'].min():.4f}, "
                     f"{adata.obs['dpt_pseudotime'].max():.4f}]")

    report_text = "\n".join(lines)
    report_path = os.path.join(output_dir, "analysis_report.md")
    with open(report_path, "w") as f:
        f.write(report_text)

    return f"Analysis summary saved to: {report_path}\n\n{report_text}"


TOOL_FUNCTIONS = {
    "plot_umap_colored": plot_umap_colored,
    "plot_dotplot": plot_dotplot,
    "plot_heatmap": plot_heatmap,
    "plot_stacked_violin": plot_stacked_violin,
    "plot_cell_composition": plot_cell_composition,
    "generate_analysis_summary": generate_analysis_summary,
}
