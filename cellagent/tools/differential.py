"""Differential expression and gene set enrichment analysis tools.

Implements differential gene expression testing between conditions/clusters
and pathway enrichment analysis for biological interpretation.
"""

TOOL_DESCRIPTIONS = [
    {
        "name": "differential_expression",
        "description": (
            "Perform differential expression analysis between two groups of cells. "
            "Supports comparison between clusters, conditions, or any categorical "
            "variable in adata.obs. Uses Wilcoxon rank-sum test by default."
        ),
        "module": "cellagent.tools.differential",
        "category": "differential_expression",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "groupby", "type": "str",
             "description": "Key in adata.obs for grouping variable"},
        ],
        "optional_parameters": [
            {"name": "groups", "type": "list", "default": None,
             "description": "Specific groups to compare (default: all vs rest)"},
            {"name": "reference", "type": "str", "default": "rest",
             "description": "Reference group for comparison"},
            {"name": "method", "type": "str", "default": "wilcoxon",
             "description": "Test method: 'wilcoxon', 't-test', 't-test_overestim_var'"},
            {"name": "n_genes", "type": "int", "default": 50,
             "description": "Number of top DEGs to report"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save DEG results"},
        ],
        "returns": "str",
        "example": 'result = differential_expression(adata, groupby="condition", groups=["disease"], reference="control")',
    },
    {
        "name": "volcano_plot",
        "description": (
            "Generate a volcano plot for differential expression results. "
            "Plots log2 fold change vs -log10 adjusted p-value, highlighting "
            "significantly up- and down-regulated genes."
        ),
        "module": "cellagent.tools.differential",
        "category": "differential_expression",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with DEG results"},
        ],
        "optional_parameters": [
            {"name": "group", "type": "str", "default": None,
             "description": "Specific group to plot"},
            {"name": "log2fc_threshold", "type": "float", "default": 1.0,
             "description": "Log2 fold change threshold for significance"},
            {"name": "pval_threshold", "type": "float", "default": 0.05,
             "description": "Adjusted p-value threshold"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save plot"},
        ],
        "returns": "str",
        "example": 'result = volcano_plot(adata, group="disease")',
    },
    {
        "name": "gene_set_enrichment",
        "description": (
            "Perform gene set enrichment analysis on differentially expressed genes. "
            "Supports Gene Ontology (GO), KEGG pathways, Reactome, and other databases "
            "via the GSEApy/Enrichr interface. Returns top enriched pathways."
        ),
        "module": "cellagent.tools.differential",
        "category": "differential_expression",
        "required_parameters": [
            {"name": "gene_list", "type": "list",
             "description": "List of gene symbols for enrichment analysis"},
        ],
        "optional_parameters": [
            {"name": "organism", "type": "str", "default": "human",
             "description": "Organism: 'human' or 'mouse'"},
            {"name": "gene_sets", "type": "str", "default": "GO_Biological_Process_2023",
             "description": "Gene set database to use"},
            {"name": "top_k", "type": "int", "default": 20,
             "description": "Number of top pathways to return"},
            {"name": "background", "type": "list", "default": None,
             "description": "Background gene list"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save enrichment plots"},
        ],
        "returns": "str",
        "example": 'result = gene_set_enrichment(["BRCA1", "TP53", "EGFR"], gene_sets="KEGG_2021_Human")',
    },
    {
        "name": "gsea_prerank",
        "description": (
            "Run Gene Set Enrichment Analysis (GSEA) using pre-ranked gene list. "
            "Genes are ranked by their differential expression statistic. "
            "Identifies pathways enriched at the top or bottom of the ranked list."
        ),
        "module": "cellagent.tools.differential",
        "category": "differential_expression",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with DEG results"},
            {"name": "group", "type": "str", "description": "Group name from DEG analysis"},
        ],
        "optional_parameters": [
            {"name": "gene_sets", "type": "str", "default": "GO_Biological_Process_2023",
             "description": "Gene set database"},
            {"name": "n_top", "type": "int", "default": 20,
             "description": "Number of top pathways to report"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save GSEA plots"},
        ],
        "returns": "str",
        "example": 'result = gsea_prerank(adata, group="0", gene_sets="KEGG_2021_Human")',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def differential_expression(
    adata, groupby: str, groups=None, reference: str = "rest",
    method: str = "wilcoxon", n_genes: int = 50, output_dir: str = "./output"
) -> str:
    """Perform differential expression analysis."""
    import scanpy as sc
    import pandas as pd
    import os

    os.makedirs(output_dir, exist_ok=True)

    if groupby not in adata.obs:
        return f"Error: Group key '{groupby}' not found in adata.obs."

    sc.tl.rank_genes_groups(
        adata, groupby, groups=groups, reference=reference,
        method=method, n_genes=n_genes
    )

    result = adata.uns["rank_genes_groups"]
    group_names = result["names"].dtype.names

    report = f"Differential Expression Analysis:\n"
    report += f"  Method: {method}\n"
    report += f"  Groupby: {groupby}\n"
    report += f"  Reference: {reference}\n\n"

    all_degs = []
    for group in group_names:
        genes = list(result["names"][group][:n_genes])
        scores = list(result["scores"][group][:n_genes])
        pvals = list(result["pvals_adj"][group][:n_genes])
        logfc = list(result["logfoldchanges"][group][:n_genes])

        n_sig = sum(1 for p in pvals if p < 0.05)
        report += f"  Group {group}: {n_sig} significant DEGs (padj<0.05)\n"
        report += f"    Top 5 upregulated:\n"
        for i in range(min(5, len(genes))):
            report += f"      {genes[i]}: logFC={logfc[i]:.3f}, padj={pvals[i]:.2e}\n"

        # Save full results
        for i in range(len(genes)):
            all_degs.append({
                "group": group, "gene": genes[i],
                "score": scores[i], "logfoldchange": logfc[i],
                "pval_adj": pvals[i],
            })

    # Save to CSV
    df = pd.DataFrame(all_degs)
    csv_path = os.path.join(output_dir, "deg_results.csv")
    df.to_csv(csv_path, index=False)
    report += f"\n  Full results saved to: {csv_path}\n"

    return report


def volcano_plot(adata, group=None, log2fc_threshold: float = 1.0,
                 pval_threshold: float = 0.05, output_dir: str = "./output") -> str:
    """Generate volcano plot."""
    import numpy as np
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if "rank_genes_groups" not in adata.uns:
        return "Error: Differential expression results not found. Run differential_expression first."

    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    if group is None:
        group = groups[0]
    if group not in groups:
        return f"Error: Group '{group}' not found in differential expression results."

    genes = np.array(result["names"][group])
    logfc = np.array(result["logfoldchanges"][group])
    pvals = np.array(result["pvals_adj"][group])

    # Replace 0 p-values
    pvals[pvals == 0] = 1e-300
    neg_log10_p = -np.log10(pvals)

    # Classify genes
    sig_up = (logfc > log2fc_threshold) & (pvals < pval_threshold)
    sig_down = (logfc < -log2fc_threshold) & (pvals < pval_threshold)
    not_sig = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(logfc[not_sig], neg_log10_p[not_sig], c="gray", alpha=0.3, s=5, label="NS")
    ax.scatter(logfc[sig_up], neg_log10_p[sig_up], c="red", alpha=0.6, s=10, label=f"Up ({sig_up.sum()})")
    ax.scatter(logfc[sig_down], neg_log10_p[sig_down], c="blue", alpha=0.6, s=10, label=f"Down ({sig_down.sum()})")

    # Label top genes
    top_idx = np.argsort(neg_log10_p)[-10:]
    for idx in top_idx:
        ax.annotate(genes[idx], (logfc[idx], neg_log10_p[idx]), fontsize=7, alpha=0.8)

    ax.axhline(-np.log10(pval_threshold), ls="--", c="gray", alpha=0.5)
    ax.axvline(log2fc_threshold, ls="--", c="gray", alpha=0.5)
    ax.axvline(-log2fc_threshold, ls="--", c="gray", alpha=0.5)
    ax.set_xlabel("Log2 Fold Change")
    ax.set_ylabel("-Log10 Adjusted P-value")
    ax.set_title(f"Volcano Plot - Group {group}")
    ax.legend()

    plot_path = os.path.join(output_dir, f"volcano_{group}.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()

    report = (
        f"Volcano Plot Generated:\n"
        f"  Group: {group}\n"
        f"  Upregulated (logFC>{log2fc_threshold}, padj<{pval_threshold}): {sig_up.sum()}\n"
        f"  Downregulated (logFC<-{log2fc_threshold}, padj<{pval_threshold}): {sig_down.sum()}\n"
        f"  Plot saved to: {plot_path}\n"
    )
    return report


def gene_set_enrichment(
    gene_list: list, organism: str = "human",
    gene_sets: str = "GO_Biological_Process_2023",
    top_k: int = 20, background: list = None,
    output_dir: str = "./output"
) -> str:
    """Perform gene set enrichment analysis using Enrichr."""
    import gseapy as gp
    import os

    os.makedirs(output_dir, exist_ok=True)

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=os.path.join(output_dir, "enrichr"),
            background=background,
            no_plot=False,
        )

        results = enr.results.head(top_k)

        report = f"Gene Set Enrichment Analysis:\n"
        report += f"  Database: {gene_sets}\n"
        report += f"  Input genes: {len(gene_list)}\n"
        report += f"  Significant pathways (padj<0.05): {(results['Adjusted P-value'] < 0.05).sum()}\n\n"

        for _, row in results.iterrows():
            report += (
                f"  {row['Term']}\n"
                f"    P-value: {row['Adjusted P-value']:.2e}, "
                f"Overlap: {row['Overlap']}, "
                f"Combined Score: {row['Combined Score']:.2f}\n"
            )

        report += f"\n  Results saved to: {output_dir}/enrichr/\n"

    except Exception as e:
        report = f"Enrichment analysis failed: {str(e)}\n"
        report += "Tip: Check internet connection and gene names.\n"

    return report


def gsea_prerank(
    adata, group: str, gene_sets: str = "GO_Biological_Process_2023",
    n_top: int = 20, output_dir: str = "./output"
) -> str:
    """Run GSEA with pre-ranked gene list."""
    import gseapy as gp
    import pandas as pd
    import numpy as np
    import os

    os.makedirs(output_dir, exist_ok=True)

    if "rank_genes_groups" not in adata.uns:
        return "Error: Differential expression results not found. Run differential_expression first."

    result = adata.uns["rank_genes_groups"]
    if group not in result["names"].dtype.names:
        return f"Error: Group '{group}' not found in differential expression results."

    genes = np.array(result["names"][group])
    scores = np.array(result["scores"][group])

    # Create ranked gene list
    rnk = pd.DataFrame({"gene": genes, "score": scores})
    rnk = rnk.sort_values("score", ascending=False)
    rnk = rnk.drop_duplicates(subset="gene")

    try:
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            outdir=os.path.join(output_dir, "gsea_prerank"),
            min_size=5,
            max_size=500,
            permutation_num=1000,
            no_plot=False,
        )

        results = pre_res.res2d.head(n_top)
        report = f"GSEA Pre-ranked Results:\n"
        report += f"  Group: {group}\n"
        report += f"  Gene sets: {gene_sets}\n\n"

        for _, row in results.iterrows():
            report += (
                f"  {row['Term']}\n"
                f"    NES: {row['NES']:.3f}, "
                f"FDR: {row['FDR q-val']:.3e}\n"
            )

        report += f"\n  Results saved to: {output_dir}/gsea_prerank/\n"

    except Exception as e:
        report = f"GSEA pre-rank analysis failed: {str(e)}\n"

    return report


TOOL_FUNCTIONS = {
    "differential_expression": differential_expression,
    "volcano_plot": volcano_plot,
    "gene_set_enrichment": gene_set_enrichment,
    "gsea_prerank": gsea_prerank,
}
