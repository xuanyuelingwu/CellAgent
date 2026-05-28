"""Trajectory and pseudotime analysis tools.

Implements trajectory inference, pseudotime ordering, and RNA velocity
analysis for understanding cellular dynamics and differentiation.
"""

TOOL_DESCRIPTIONS = [
    {
        "name": "run_diffusion_map",
        "description": (
            "Compute diffusion map embedding for trajectory analysis. "
            "Diffusion maps capture the global structure of the data and are "
            "well-suited for identifying continuous developmental trajectories."
        ),
        "module": "cellagent.tools.trajectory",
        "category": "trajectory",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with neighbors computed"},
        ],
        "optional_parameters": [
            {"name": "n_comps", "type": "int", "default": 15,
             "description": "Number of diffusion components"},
        ],
        "returns": "str",
        "example": 'result = run_diffusion_map(adata, n_comps=15)',
    },
    {
        "name": "compute_pseudotime_dpt",
        "description": (
            "Compute diffusion pseudotime (DPT) ordering of cells. "
            "Requires a root cell to be specified (typically a stem/progenitor cell). "
            "DPT measures the developmental distance from the root cell."
        ),
        "module": "cellagent.tools.trajectory",
        "category": "trajectory",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with diffusion map"},
            {"name": "root_cell_type", "type": "str",
             "description": "Cell type or cluster to use as root (e.g., 'HSC', 'stem cells')"},
        ],
        "optional_parameters": [
            {"name": "cell_type_key", "type": "str", "default": "cell_type_llm",
             "description": "Key in adata.obs for cell type labels"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save pseudotime plots"},
        ],
        "returns": "str",
        "example": 'result = compute_pseudotime_dpt(adata, root_cell_type="HSC")',
    },
    {
        "name": "run_paga",
        "description": (
            "Run PAGA (Partition-based Graph Abstraction) to infer trajectory topology. "
            "PAGA provides a coarse-grained map of the data manifold, showing "
            "connectivity between cell clusters. Useful for identifying lineage relationships."
        ),
        "module": "cellagent.tools.trajectory",
        "category": "trajectory",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with clustering"},
        ],
        "optional_parameters": [
            {"name": "groups", "type": "str", "default": "leiden",
             "description": "Clustering key for PAGA"},
            {"name": "threshold", "type": "float", "default": 0.05,
             "description": "Threshold for PAGA connectivity"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save PAGA plots"},
        ],
        "returns": "str",
        "example": 'result = run_paga(adata, groups="leiden")',
    },
    {
        "name": "identify_trajectory_genes",
        "description": (
            "Identify genes that vary along the pseudotime trajectory. "
            "Fits gene expression as a function of pseudotime and identifies "
            "genes with significant temporal patterns."
        ),
        "module": "cellagent.tools.trajectory",
        "category": "trajectory",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with pseudotime"},
        ],
        "optional_parameters": [
            {"name": "n_genes", "type": "int", "default": 50,
             "description": "Number of top trajectory genes to return"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save trajectory gene plots"},
        ],
        "returns": "str",
        "example": 'result = identify_trajectory_genes(adata, n_genes=50)',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def run_diffusion_map(adata, n_comps: int = 15) -> str:
    """Compute diffusion map."""
    import scanpy as sc

    if "neighbors" not in adata.uns:
        return "Error: Neighbor graph not found. Run compute_neighbors before run_diffusion_map."

    sc.tl.diffmap(adata, n_comps=n_comps)

    report = (
        f"Diffusion Map Computed:\n"
        f"  Components: {n_comps}\n"
        f"  Embedding shape: {adata.obsm['X_diffmap'].shape}\n"
    )
    return report


def compute_pseudotime_dpt(adata, root_cell_type: str,
                           cell_type_key: str = "cell_type_llm",
                           output_dir: str = "./output") -> str:
    """Compute diffusion pseudotime."""
    import scanpy as sc
    import numpy as np
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if "X_diffmap" not in adata.obsm:
        return "Error: Diffusion map not found. Run run_diffusion_map before compute_pseudotime_dpt."

    # Find root cell
    if cell_type_key in adata.obs.columns:
        root_mask = adata.obs[cell_type_key].str.contains(root_cell_type, case=False, na=False)
        if root_mask.sum() == 0:
            # Try matching with cluster labels
            for key in ["leiden", "louvain"]:
                if key in adata.obs.columns:
                    root_mask = adata.obs[key] == root_cell_type
                    if root_mask.sum() > 0:
                        break

        if root_mask.sum() > 0:
            root_idx = np.where(root_mask)[0]
            # Pick the cell closest to the centroid of root cells
            root_embedding = adata.obsm["X_diffmap"][root_idx]
            centroid = root_embedding.mean(axis=0)
            distances = np.linalg.norm(root_embedding - centroid, axis=1)
            adata.uns["iroot"] = root_idx[np.argmin(distances)]
        else:
            adata.uns["iroot"] = 0
            return f"Warning: Could not find root cell type '{root_cell_type}'. Using cell 0 as root.\n"
    else:
        adata.uns["iroot"] = 0

    sc.tl.dpt(adata)

    # Plot pseudotime on UMAP
    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color="dpt_pseudotime", show=False)
        plt.savefig(os.path.join(output_dir, "pseudotime_umap.png"), dpi=150, bbox_inches="tight")
        plt.close()

    report = (
        f"Diffusion Pseudotime Computed:\n"
        f"  Root cell type: {root_cell_type}\n"
        f"  Root cell index: {adata.uns['iroot']}\n"
        f"  Pseudotime range: [{adata.obs['dpt_pseudotime'].min():.4f}, "
        f"{adata.obs['dpt_pseudotime'].max():.4f}]\n"
        f"  Pseudotime stored in: adata.obs['dpt_pseudotime']\n"
    )
    if "X_umap" in adata.obsm:
        report += f"  Pseudotime UMAP saved to: {output_dir}/pseudotime_umap.png\n"
    return report


def run_paga(adata, groups: str = "leiden", threshold: float = 0.05,
             output_dir: str = "./output") -> str:
    """Run PAGA trajectory analysis."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if groups not in adata.obs:
        return f"Error: Group key '{groups}' not found in adata.obs."
    if "neighbors" not in adata.uns:
        return "Error: Neighbor graph not found. Run compute_neighbors before run_paga."

    sc.tl.paga(adata, groups=groups)

    # Plot PAGA graph
    sc.pl.paga(adata, threshold=threshold, show=False)
    plt.savefig(os.path.join(output_dir, "paga_graph.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # PAGA-initialized UMAP
    sc.tl.umap(adata, init_pos="paga")
    sc.pl.umap(adata, color=groups, show=False)
    plt.savefig(os.path.join(output_dir, "paga_umap.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Extract connectivity
    connectivity = adata.uns["paga"]["connectivities"].toarray()
    n_clusters = connectivity.shape[0]

    report = f"PAGA Trajectory Analysis:\n"
    report += f"  Groups: {groups}\n"
    report += f"  Number of groups: {n_clusters}\n"
    report += f"  Strong connections (>{threshold}):\n"

    for i in range(n_clusters):
        for j in range(i + 1, n_clusters):
            if connectivity[i, j] > threshold:
                report += f"    {i} <-> {j}: {connectivity[i, j]:.3f}\n"

    report += f"\n  PAGA graph saved to: {output_dir}/paga_graph.png\n"
    report += f"  PAGA-UMAP saved to: {output_dir}/paga_umap.png\n"
    return report


def identify_trajectory_genes(adata, n_genes: int = 50,
                               output_dir: str = "./output") -> str:
    """Identify genes varying along pseudotime."""
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr

    os.makedirs(output_dir, exist_ok=True)

    if "dpt_pseudotime" not in adata.obs.columns:
        return "Error: Pseudotime not computed. Run compute_pseudotime_dpt first."

    pseudotime = adata.obs["dpt_pseudotime"].values
    valid = ~np.isinf(pseudotime) & ~np.isnan(pseudotime)

    # Compute correlation of each gene with pseudotime
    if adata.raw is not None:
        gene_expr = adata.raw[valid].X
        gene_names = adata.raw.var_names
    else:
        gene_expr = adata[valid].X
        gene_names = adata.var_names

    if hasattr(gene_expr, "toarray"):
        gene_expr = gene_expr.toarray()

    pt_valid = pseudotime[valid]
    correlations = []
    for i in range(gene_expr.shape[1]):
        if gene_expr[:, i].std() > 0:
            corr, pval = spearmanr(pt_valid, gene_expr[:, i])
            correlations.append({"gene": gene_names[i], "correlation": corr, "pvalue": pval})

    df = pd.DataFrame(correlations)
    df = df.dropna()
    df["abs_corr"] = df["correlation"].abs()
    df = df.sort_values("abs_corr", ascending=False).head(n_genes)

    # Save results
    csv_path = os.path.join(output_dir, "trajectory_genes.csv")
    df.to_csv(csv_path, index=False)

    report = f"Trajectory Gene Analysis:\n"
    report += f"  Genes tested: {len(correlations)}\n"
    report += f"  Top {n_genes} trajectory genes:\n\n"

    for _, row in df.head(20).iterrows():
        direction = "+" if row["correlation"] > 0 else "-"
        report += (
            f"    {row['gene']}: corr={row['correlation']:.4f} ({direction}), "
            f"p={row['pvalue']:.2e}\n"
        )

    report += f"\n  Full results saved to: {csv_path}\n"
    return report


TOOL_FUNCTIONS = {
    "run_diffusion_map": run_diffusion_map,
    "compute_pseudotime_dpt": compute_pseudotime_dpt,
    "run_paga": run_paga,
    "identify_trajectory_genes": identify_trajectory_genes,
}
