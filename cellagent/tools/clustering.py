"""Clustering tools for single-cell RNA-seq data.

Implements neighborhood graph construction, clustering algorithms,
and dimensionality reduction for visualization.
"""

TOOL_DESCRIPTIONS = [
    {
        "name": "compute_neighbors",
        "description": (
            "Compute the k-nearest neighbor graph on PCA embeddings. "
            "This is a prerequisite for clustering and UMAP/tSNE visualization. "
            "Uses approximate nearest neighbors for scalability."
        ),
        "module": "cellagent.tools.clustering",
        "category": "clustering",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with PCA computed"},
        ],
        "optional_parameters": [
            {"name": "n_neighbors", "type": "int", "default": 15,
             "description": "Number of nearest neighbors"},
            {"name": "n_pcs", "type": "int", "default": 30,
             "description": "Number of PCs to use"},
        ],
        "returns": "str",
        "example": 'result = compute_neighbors(adata, n_neighbors=15, n_pcs=30)',
    },
    {
        "name": "run_umap",
        "description": (
            "Compute UMAP (Uniform Manifold Approximation and Projection) embedding "
            "for 2D visualization. Requires neighbor graph to be computed first."
        ),
        "module": "cellagent.tools.clustering",
        "category": "clustering",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with neighbors computed"},
        ],
        "optional_parameters": [
            {"name": "min_dist", "type": "float", "default": 0.5,
             "description": "Minimum distance between points in UMAP"},
            {"name": "spread", "type": "float", "default": 1.0,
             "description": "Spread of UMAP embedding"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save UMAP plot"},
        ],
        "returns": "str",
        "example": 'result = run_umap(adata, min_dist=0.3)',
    },
    {
        "name": "run_tsne",
        "description": (
            "Compute t-SNE (t-distributed Stochastic Neighbor Embedding) for "
            "2D visualization. Alternative to UMAP, better for preserving local structure."
        ),
        "module": "cellagent.tools.clustering",
        "category": "clustering",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with PCA computed"},
        ],
        "optional_parameters": [
            {"name": "n_pcs", "type": "int", "default": 30,
             "description": "Number of PCs to use"},
            {"name": "perplexity", "type": "float", "default": 30,
             "description": "Perplexity parameter for t-SNE"},
        ],
        "returns": "str",
        "example": 'result = run_tsne(adata, perplexity=30)',
    },
    {
        "name": "run_leiden_clustering",
        "description": (
            "Perform Leiden clustering on the neighbor graph. "
            "Leiden is the recommended clustering algorithm for scRNA-seq data, "
            "providing better-connected communities than Louvain. "
            "Resolution parameter controls cluster granularity."
        ),
        "module": "cellagent.tools.clustering",
        "category": "clustering",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with neighbors computed"},
        ],
        "optional_parameters": [
            {"name": "resolution", "type": "float", "default": 0.5,
             "description": "Resolution parameter (higher = more clusters)"},
            {"name": "key_added", "type": "str", "default": "leiden",
             "description": "Key to store clustering results in adata.obs"},
        ],
        "returns": "str",
        "example": 'result = run_leiden_clustering(adata, resolution=0.5)',
    },
    {
        "name": "run_louvain_clustering",
        "description": (
            "Perform Louvain clustering on the neighbor graph. "
            "Alternative to Leiden clustering. Generally Leiden is preferred."
        ),
        "module": "cellagent.tools.clustering",
        "category": "clustering",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with neighbors computed"},
        ],
        "optional_parameters": [
            {"name": "resolution", "type": "float", "default": 0.5,
             "description": "Resolution parameter"},
            {"name": "key_added", "type": "str", "default": "louvain",
             "description": "Key to store clustering results"},
        ],
        "returns": "str",
        "example": 'result = run_louvain_clustering(adata, resolution=1.0)',
    },
    {
        "name": "evaluate_clustering_quality",
        "description": (
            "Evaluate clustering quality using multiple metrics: "
            "silhouette score, Calinski-Harabasz index, and cluster size distribution. "
            "Helps determine if the clustering resolution is appropriate."
        ),
        "module": "cellagent.tools.clustering",
        "category": "clustering",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with clustering"},
        ],
        "optional_parameters": [
            {"name": "cluster_key", "type": "str", "default": "leiden",
             "description": "Key in adata.obs containing cluster labels"},
            {"name": "embedding_key", "type": "str", "default": "X_pca",
             "description": "Key in adata.obsm for embedding to evaluate"},
        ],
        "returns": "str",
        "example": 'result = evaluate_clustering_quality(adata, cluster_key="leiden")',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def compute_neighbors(adata, n_neighbors: int = 15, n_pcs: int = 30) -> str:
    """Compute nearest neighbor graph."""
    import scanpy as sc

    if "X_pca" not in adata.obsm:
        return "Error: PCA embedding not found. Run run_pca before compute_neighbors."

    max_pcs = adata.obsm["X_pca"].shape[1]
    n_pcs = min(n_pcs, max_pcs)
    n_neighbors = min(n_neighbors, max(1, adata.n_obs - 1))

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    report = (
        f"Neighbor Graph Computed:\n"
        f"  n_neighbors: {n_neighbors}\n"
        f"  n_pcs: {n_pcs}\n"
        f"  Cells: {adata.n_obs}\n"
    )
    return report


def run_umap(adata, min_dist: float = 0.5, spread: float = 1.0,
             output_dir: str = "./output") -> str:
    """Compute UMAP embedding."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    if "neighbors" not in adata.uns:
        return "Error: Neighbor graph not found. Run compute_neighbors before run_umap."

    sc.tl.umap(adata, min_dist=min_dist, spread=spread)

    # Save UMAP plot
    sc.pl.umap(adata, show=False)
    plt.savefig(os.path.join(output_dir, "umap.png"), dpi=150, bbox_inches="tight")
    plt.close()

    report = (
        f"UMAP Computed:\n"
        f"  min_dist: {min_dist}\n"
        f"  spread: {spread}\n"
        f"  Embedding shape: {adata.obsm['X_umap'].shape}\n"
        f"  Plot saved to: {output_dir}/umap.png\n"
    )
    return report


def run_tsne(adata, n_pcs: int = 30, perplexity: float = 30) -> str:
    """Compute t-SNE embedding."""
    import scanpy as sc

    if "X_pca" not in adata.obsm:
        return "Error: PCA embedding not found. Run run_pca before run_tsne."
    n_pcs = min(n_pcs, adata.obsm["X_pca"].shape[1])
    perplexity = min(perplexity, max(1, (adata.n_obs - 1) / 3))

    sc.tl.tsne(adata, n_pcs=n_pcs, perplexity=perplexity)
    report = (
        f"t-SNE Computed:\n"
        f"  n_pcs: {n_pcs}\n"
        f"  perplexity: {perplexity}\n"
        f"  Embedding shape: {adata.obsm['X_tsne'].shape}\n"
    )
    return report


def run_leiden_clustering(adata, resolution: float = 0.5,
                          key_added: str = "leiden") -> str:
    """Run Leiden clustering."""
    import scanpy as sc

    if "neighbors" not in adata.uns:
        return "Error: Neighbor graph not found. Run compute_neighbors before clustering."

    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    n_clusters = adata.obs[key_added].nunique()
    cluster_sizes = adata.obs[key_added].value_counts().sort_index()

    report = (
        f"Leiden Clustering Results:\n"
        f"  Resolution: {resolution}\n"
        f"  Number of clusters: {n_clusters}\n"
        f"  Cluster sizes:\n"
    )
    for cluster, size in cluster_sizes.items():
        report += f"    Cluster {cluster}: {size} cells ({size/adata.n_obs*100:.1f}%)\n"
    return report


def run_louvain_clustering(adata, resolution: float = 0.5,
                           key_added: str = "louvain") -> str:
    """Run Louvain clustering."""
    import scanpy as sc

    if "neighbors" not in adata.uns:
        return "Error: Neighbor graph not found. Run compute_neighbors before clustering."

    sc.tl.louvain(adata, resolution=resolution, key_added=key_added)
    n_clusters = adata.obs[key_added].nunique()
    report = (
        f"Louvain Clustering Results:\n"
        f"  Resolution: {resolution}\n"
        f"  Number of clusters: {n_clusters}\n"
    )
    return report


def evaluate_clustering_quality(adata, cluster_key: str = "leiden",
                                 embedding_key: str = "X_pca") -> str:
    """Evaluate clustering quality."""
    import numpy as np
    from sklearn.metrics import silhouette_score, calinski_harabasz_score

    if cluster_key not in adata.obs:
        return f"Error: Cluster key '{cluster_key}' not found in adata.obs."
    if embedding_key not in adata.obsm:
        return f"Error: Embedding '{embedding_key}' not found in adata.obsm."

    from sklearn.preprocessing import LabelEncoder

    labels = LabelEncoder().fit_transform(adata.obs[cluster_key].astype(str).values)
    embedding = adata.obsm[embedding_key][:, :30]  # Use first 30 PCs

    n_clusters = len(np.unique(labels))
    if n_clusters < 2:
        return (
            f"Clustering Quality Evaluation:\n"
            f"  Cluster key: {cluster_key}\n"
            f"  Number of clusters: {n_clusters}\n"
            f"  Metrics skipped: at least 2 clusters are required.\n"
        )

    # Subsample for large datasets
    if adata.n_obs > 10000:
        idx = np.random.choice(adata.n_obs, 10000, replace=False)
        embedding_sub = embedding[idx]
        labels_sub = labels[idx]
    else:
        embedding_sub = embedding
        labels_sub = labels

    sil_score = silhouette_score(embedding_sub, labels_sub)
    ch_score = calinski_harabasz_score(embedding_sub, labels_sub)

    cluster_sizes = np.bincount(labels)
    size_cv = cluster_sizes.std() / cluster_sizes.mean()

    report = (
        f"Clustering Quality Evaluation:\n"
        f"  Cluster key: {cluster_key}\n"
        f"  Number of clusters: {n_clusters}\n"
        f"  Silhouette Score: {sil_score:.4f} (range: -1 to 1, higher is better)\n"
        f"  Calinski-Harabasz Index: {ch_score:.2f} (higher is better)\n"
        f"  Cluster size CV: {size_cv:.3f} (lower = more balanced)\n"
        f"  Min cluster size: {cluster_sizes.min()}\n"
        f"  Max cluster size: {cluster_sizes.max()}\n"
    )

    if sil_score > 0.3:
        report += "  Assessment: Good clustering quality\n"
    elif sil_score > 0.1:
        report += "  Assessment: Moderate clustering quality, consider adjusting resolution\n"
    else:
        report += "  Assessment: Poor clustering quality, try different resolution or preprocessing\n"

    return report


TOOL_FUNCTIONS = {
    "compute_neighbors": compute_neighbors,
    "run_umap": run_umap,
    "run_tsne": run_tsne,
    "run_leiden_clustering": run_leiden_clustering,
    "run_louvain_clustering": run_louvain_clustering,
    "evaluate_clustering_quality": evaluate_clustering_quality,
}
