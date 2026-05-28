"""Cell type annotation tools for single-cell RNA-seq data.

Implements marker-based, automated, and LLM-assisted cell type annotation
following best practices from the single-cell community.
"""

TOOL_DESCRIPTIONS = [
    {
        "name": "find_marker_genes",
        "description": (
            "Identify marker genes for each cluster using statistical tests. "
            "Computes differentially expressed genes between each cluster and the rest "
            "using Wilcoxon rank-sum test (recommended) or t-test. "
            "Returns top marker genes per cluster with statistics."
        ),
        "module": "cellagent.tools.annotation",
        "category": "annotation",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with clustering"},
        ],
        "optional_parameters": [
            {"name": "groupby", "type": "str", "default": "leiden",
             "description": "Key in adata.obs for cluster labels"},
            {"name": "method", "type": "str", "default": "wilcoxon",
             "description": "Statistical test: 'wilcoxon', 't-test', 't-test_overestim_var'"},
            {"name": "n_genes", "type": "int", "default": 25,
             "description": "Number of top marker genes to report per cluster"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save marker gene plots"},
        ],
        "returns": "str",
        "example": 'result = find_marker_genes(adata, groupby="leiden", n_genes=25)',
    },
    {
        "name": "annotate_with_markers",
        "description": (
            "Annotate cell types based on known marker gene expression. "
            "Computes dot plots and scores for provided marker gene sets, "
            "then assigns cell type labels based on highest scoring markers. "
            "Supports custom marker dictionaries for any tissue type."
        ),
        "module": "cellagent.tools.annotation",
        "category": "annotation",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with clustering"},
            {"name": "marker_dict", "type": "dict",
             "description": "Dictionary mapping cell types to marker gene lists, e.g., {'T cells': ['CD3D', 'CD3E']}"},
        ],
        "optional_parameters": [
            {"name": "groupby", "type": "str", "default": "leiden",
             "description": "Cluster key for annotation"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save annotation plots"},
        ],
        "returns": "str",
        "example": (
            'markers = {"T cells": ["CD3D", "CD3E"], "B cells": ["CD19", "MS4A1"]}\n'
            'result = annotate_with_markers(adata, markers, groupby="leiden")'
        ),
    },
    {
        "name": "annotate_with_llm",
        "description": (
            "Use an LLM to annotate cell types based on differentially expressed genes. "
            "Extracts top marker genes per cluster, sends them to the LLM with tissue "
            "context, and returns predicted cell type labels with confidence scores. "
            "This leverages the LLM's knowledge of cell type markers across tissues."
        ),
        "module": "cellagent.tools.annotation",
        "category": "annotation",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with clustering and DEGs"},
            {"name": "tissue_info", "type": "str",
             "description": "Tissue/sample information, e.g., 'human PBMC', 'mouse brain cortex'"},
        ],
        "optional_parameters": [
            {"name": "groupby", "type": "str", "default": "leiden",
             "description": "Cluster key"},
            {"name": "n_genes", "type": "int", "default": 15,
             "description": "Number of top marker genes to use per cluster"},
        ],
        "returns": "str",
        "example": 'result = annotate_with_llm(adata, tissue_info="human PBMC")',
    },
    {
        "name": "score_cell_types",
        "description": (
            "Score cells for specific cell type signatures using gene set scoring. "
            "Computes a score for each cell based on the expression of a set of marker genes. "
            "Uses Scanpy's score_genes method (based on Tirosh et al. 2016)."
        ),
        "module": "cellagent.tools.annotation",
        "category": "annotation",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "gene_sets", "type": "dict",
             "description": "Dictionary mapping cell type names to gene lists"},
        ],
        "optional_parameters": [
            {"name": "ctrl_size", "type": "int", "default": 50,
             "description": "Number of control genes for scoring"},
        ],
        "returns": "str",
        "example": (
            'gene_sets = {"T_cell_score": ["CD3D", "CD3E", "CD2"]}\n'
            'result = score_cell_types(adata, gene_sets)'
        ),
    },
    {
        "name": "compare_annotations",
        "description": (
            "Compare two sets of cell type annotations to assess agreement. "
            "Computes confusion matrix, adjusted Rand index, and normalized mutual information. "
            "Useful for validating automated annotations against manual labels."
        ),
        "module": "cellagent.tools.annotation",
        "category": "annotation",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "key1", "type": "str", "description": "First annotation key in adata.obs"},
            {"name": "key2", "type": "str", "description": "Second annotation key in adata.obs"},
        ],
        "optional_parameters": [],
        "returns": "str",
        "example": 'result = compare_annotations(adata, "manual_annotation", "llm_annotation")',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def find_marker_genes(adata, groupby: str = "leiden", method: str = "wilcoxon",
                      n_genes: int = 25, output_dir: str = "./output") -> str:
    """Find marker genes for each cluster."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if groupby not in adata.obs:
        return f"Error: Group key '{groupby}' not found in adata.obs."

    sc.tl.rank_genes_groups(adata, groupby, method=method)

    # Generate marker gene plot
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, show=False)
    plt.savefig(os.path.join(output_dir, "marker_genes.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Extract top markers per cluster
    result = adata.uns["rank_genes_groups"]
    clusters = result["names"].dtype.names

    report = f"Marker Gene Analysis (method={method}):\n"
    for cluster in clusters:
        top_genes = list(result["names"][cluster][:10])
        top_scores = list(result["scores"][cluster][:10])
        top_pvals = list(result["pvals_adj"][cluster][:10])
        report += f"\n  Cluster {cluster}:\n"
        for g, s, p in zip(top_genes[:5], top_scores[:5], top_pvals[:5]):
            report += f"    {g}: score={s:.3f}, padj={p:.2e}\n"

    report += f"\n  Marker gene plot saved to: {output_dir}/marker_genes.png\n"
    return report


def annotate_with_markers(adata, marker_dict: dict, groupby: str = "leiden",
                          output_dir: str = "./output") -> str:
    """Annotate clusters using known marker genes."""
    import scanpy as sc
    import numpy as np
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if groupby not in adata.obs:
        return f"Error: Group key '{groupby}' not found in adata.obs."

    # Filter markers to those present in the data
    valid_markers = {}
    for ct, genes in marker_dict.items():
        valid = [g for g in genes if g in adata.var_names]
        if valid:
            valid_markers[ct] = valid

    if not valid_markers:
        return "Error: No marker genes found in the dataset. Check gene names."

    # Generate dot plot
    sc.pl.dotplot(adata, valid_markers, groupby=groupby, show=False)
    plt.savefig(os.path.join(output_dir, "marker_dotplot.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Score each cluster for each cell type
    clusters = adata.obs[groupby].unique()
    cluster_annotations = {}

    for ct, genes in valid_markers.items():
        score_name = f"score_{ct.replace(' ', '_')}"
        sc.tl.score_genes(adata, genes, score_name=score_name)

    for cluster in sorted(clusters, key=lambda x: int(x) if x.isdigit() else x):
        mask = adata.obs[groupby] == cluster
        best_score = -np.inf
        best_ct = "Unknown"
        for ct, genes in valid_markers.items():
            score_name = f"score_{ct.replace(' ', '_')}"
            mean_score = adata.obs.loc[mask, score_name].mean()
            if mean_score > best_score:
                best_score = mean_score
                best_ct = ct
        cluster_annotations[cluster] = (best_ct, best_score)

    # Apply annotations
    adata.obs["cell_type_marker"] = adata.obs[groupby].map(
        {k: v[0] for k, v in cluster_annotations.items()}
    )

    report = "Marker-Based Annotation Results:\n"
    for cluster, (ct, score) in sorted(
        cluster_annotations.items(), key=lambda x: int(x[0]) if x[0].isdigit() else x[0]
    ):
        n_cells = (adata.obs[groupby] == cluster).sum()
        report += f"  Cluster {cluster}: {ct} (score={score:.3f}, n={n_cells})\n"

    report += f"\n  Dot plot saved to: {output_dir}/marker_dotplot.png\n"
    report += "  Annotations stored in adata.obs['cell_type_marker']\n"
    return report


def annotate_with_llm(adata, tissue_info: str, groupby: str = "leiden",
                      n_genes: int = 15) -> str:
    """Annotate cell types using LLM."""
    import scanpy as sc
    import json
    from cellagent.llm import llm_chat

    # Ensure DEGs are computed
    if groupby not in adata.obs:
        return f"Error: Group key '{groupby}' not found in adata.obs."

    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby, method="wilcoxon")

    result = adata.uns["rank_genes_groups"]
    clusters = result["names"].dtype.names

    # Build marker gene info for each cluster
    cluster_info = {}
    for cluster in clusters:
        top_genes = list(result["names"][cluster][:n_genes])
        top_scores = list(result["scores"][cluster][:n_genes])
        cluster_info[cluster] = {
            "top_genes": top_genes,
            "scores": [f"{s:.2f}" for s in top_scores],
            "n_cells": int((adata.obs[groupby] == cluster).sum()),
        }

    prompt = f"""You are an expert single-cell biologist. Based on the following marker genes 
for each cluster from a {tissue_info} sample, predict the cell type for each cluster.

Cluster information:
{json.dumps(cluster_info, indent=2)}

For each cluster, provide:
1. The predicted cell type name
2. Confidence level (high/medium/low)
3. Key marker genes supporting your prediction

Respond in JSON format:
{{
  "annotations": {{
    "cluster_id": {{
      "cell_type": "predicted type",
      "confidence": "high/medium/low",
      "key_markers": ["gene1", "gene2"]
    }}
  }}
}}
"""

    response = llm_chat(
        messages=[{"role": "user", "content": prompt}],
        temperature=0.1,
        max_tokens=4096,
    )

    # Parse LLM response
    try:
        # Try to extract JSON from the response
        json_match = response
        if "```json" in response:
            json_match = response.split("```json")[1].split("```")[0]
        elif "```" in response:
            json_match = response.split("```")[1].split("```")[0]

        annotations = json.loads(json_match)
        if "annotations" in annotations:
            annotations = annotations["annotations"]

        # Apply annotations
        cluster_to_type = {}
        report = f"LLM-Based Cell Type Annotation ({tissue_info}):\n"
        for cluster_id, info in annotations.items():
            ct = info.get("cell_type", "Unknown")
            conf = info.get("confidence", "unknown")
            markers = info.get("key_markers", [])
            cluster_to_type[str(cluster_id)] = ct
            n_cells = cluster_info.get(str(cluster_id), {}).get("n_cells", 0)
            report += (
                f"  Cluster {cluster_id}: {ct} "
                f"(confidence={conf}, n={n_cells})\n"
                f"    Key markers: {', '.join(markers[:5])}\n"
            )

        adata.obs["cell_type_llm"] = adata.obs[groupby].astype(str).map(cluster_to_type)
        report += "\n  Annotations stored in adata.obs['cell_type_llm']\n"

    except (json.JSONDecodeError, KeyError, IndexError) as e:
        report = f"LLM annotation completed but parsing failed: {e}\n"
        report += f"Raw LLM response:\n{response[:1000]}\n"

    return report


def score_cell_types(adata, gene_sets: dict, ctrl_size: int = 50) -> str:
    """Score cells for cell type signatures."""
    import scanpy as sc

    report = "Cell Type Scoring Results:\n"
    for name, genes in gene_sets.items():
        valid_genes = [g for g in genes if g in adata.var_names]
        if not valid_genes:
            report += f"  {name}: No valid genes found\n"
            continue

        score_name = f"score_{name}"
        sc.tl.score_genes(adata, valid_genes, score_name=score_name, ctrl_size=ctrl_size)

        mean_score = adata.obs[score_name].mean()
        std_score = adata.obs[score_name].std()
        report += (
            f"  {name}: mean={mean_score:.4f}, std={std_score:.4f}, "
            f"genes_used={len(valid_genes)}/{len(genes)}\n"
        )

    return report


def compare_annotations(adata, key1: str, key2: str) -> str:
    """Compare two annotation sets."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    import pandas as pd

    missing = [key for key in (key1, key2) if key not in adata.obs]
    if missing:
        return f"Error: Annotation key(s) not found in adata.obs: {missing}"

    labels1 = adata.obs[key1].astype(str)
    labels2 = adata.obs[key2].astype(str)

    ari = adjusted_rand_score(labels1, labels2)
    nmi = normalized_mutual_info_score(labels1, labels2)

    # Confusion matrix
    ct = pd.crosstab(labels1, labels2, margins=True)

    report = (
        f"Annotation Comparison: {key1} vs {key2}\n"
        f"  Adjusted Rand Index: {ari:.4f} (1.0 = perfect agreement)\n"
        f"  Normalized Mutual Information: {nmi:.4f} (1.0 = perfect agreement)\n"
        f"  Unique labels in {key1}: {labels1.nunique()}\n"
        f"  Unique labels in {key2}: {labels2.nunique()}\n"
        f"\n  Confusion Matrix (top 10x10):\n{ct.iloc[:10, :10].to_string()}\n"
    )
    return report


TOOL_FUNCTIONS = {
    "find_marker_genes": find_marker_genes,
    "annotate_with_markers": annotate_with_markers,
    "annotate_with_llm": annotate_with_llm,
    "score_cell_types": score_cell_types,
    "compare_annotations": compare_annotations,
}
