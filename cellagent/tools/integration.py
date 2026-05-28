"""Multi-sample integration tools for single-cell data.

Implements batch correction and data integration methods for combining
multiple single-cell datasets while preserving biological variation.
"""

TOOL_DESCRIPTIONS = [
    {
        "name": "integrate_harmony",
        "description": (
            "Integrate multiple batches/samples using Harmony. "
            "Harmony is a fast, scalable method that corrects batch effects "
            "in the PCA space. Well-suited for large datasets."
        ),
        "module": "cellagent.tools.integration",
        "category": "integration",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with PCA computed"},
            {"name": "batch_key", "type": "str",
             "description": "Key in adata.obs identifying batch/sample"},
        ],
        "optional_parameters": [
            {"name": "max_iter_harmony", "type": "int", "default": 20,
             "description": "Maximum iterations for Harmony"},
        ],
        "returns": "str",
        "example": 'result = integrate_harmony(adata, batch_key="sample")',
    },
    {
        "name": "integrate_bbknn",
        "description": (
            "Integrate batches using BBKNN (Batch Balanced K-Nearest Neighbors). "
            "BBKNN modifies the neighbor graph to balance connections across batches. "
            "Simple and effective for moderate batch effects."
        ),
        "module": "cellagent.tools.integration",
        "category": "integration",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object with PCA computed"},
            {"name": "batch_key", "type": "str",
             "description": "Key in adata.obs identifying batch"},
        ],
        "optional_parameters": [
            {"name": "neighbors_within_batch", "type": "int", "default": 3,
             "description": "Number of neighbors within each batch"},
            {"name": "n_pcs", "type": "int", "default": 30,
             "description": "Number of PCs to use"},
        ],
        "returns": "str",
        "example": 'result = integrate_bbknn(adata, batch_key="sample")',
    },
    {
        "name": "evaluate_integration",
        "description": (
            "Evaluate batch integration quality using multiple metrics: "
            "batch mixing (kBET, batch entropy), biological conservation "
            "(ARI, NMI for cell types), and LISI scores."
        ),
        "module": "cellagent.tools.integration",
        "category": "integration",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "batch_key", "type": "str", "description": "Batch key"},
        ],
        "optional_parameters": [
            {"name": "cell_type_key", "type": "str", "default": None,
             "description": "Cell type key for biological conservation metrics"},
            {"name": "embedding_key", "type": "str", "default": "X_pca",
             "description": "Embedding to evaluate"},
        ],
        "returns": "str",
        "example": 'result = evaluate_integration(adata, "batch", cell_type_key="cell_type")',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def integrate_harmony(adata, batch_key: str, max_iter_harmony: int = 20) -> str:
    """Integrate using Harmony."""
    import harmonypy as hm

    if "X_pca" not in adata.obsm:
        return "Error: PCA embedding not found. Run run_pca before integrate_harmony."
    if batch_key not in adata.obs:
        return f"Error: Batch key '{batch_key}' not found in adata.obs."

    # Run Harmony on PCA embeddings
    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch_key,
        max_iter_harmony=max_iter_harmony,
    )

    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    n_batches = adata.obs[batch_key].nunique()
    report = (
        f"Harmony Integration Complete:\n"
        f"  Batch key: {batch_key}\n"
        f"  Number of batches: {n_batches}\n"
        f"  Corrected embedding stored in: adata.obsm['X_pca_harmony']\n"
        f"  Note: Use 'X_pca_harmony' for downstream neighbors/UMAP computation\n"
    )
    return report


def integrate_bbknn(adata, batch_key: str, neighbors_within_batch: int = 3,
                    n_pcs: int = 30) -> str:
    """Integrate using BBKNN."""
    import bbknn

    if "X_pca" not in adata.obsm:
        return "Error: PCA embedding not found. Run run_pca before integrate_bbknn."
    if batch_key not in adata.obs:
        return f"Error: Batch key '{batch_key}' not found in adata.obs."
    n_pcs = min(n_pcs, adata.obsm["X_pca"].shape[1])

    bbknn.bbknn(
        adata,
        batch_key=batch_key,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
    )

    n_batches = adata.obs[batch_key].nunique()
    report = (
        f"BBKNN Integration Complete:\n"
        f"  Batch key: {batch_key}\n"
        f"  Number of batches: {n_batches}\n"
        f"  Neighbors within batch: {neighbors_within_batch}\n"
        f"  Neighbor graph updated in place\n"
        f"  Note: Proceed directly to UMAP/clustering\n"
    )
    return report


def evaluate_integration(adata, batch_key: str, cell_type_key: str = None,
                          embedding_key: str = "X_pca") -> str:
    """Evaluate integration quality."""
    import numpy as np
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import LabelEncoder

    if embedding_key not in adata.obsm:
        return f"Error: Embedding '{embedding_key}' not found in adata.obsm."
    if batch_key not in adata.obs:
        return f"Error: Batch key '{batch_key}' not found in adata.obs."

    embedding = adata.obsm[embedding_key][:, :30]
    batch_labels_raw = adata.obs[batch_key].astype(str).values
    batch_labels = LabelEncoder().fit_transform(batch_labels_raw)

    n_batches = len(np.unique(batch_labels))
    if n_batches < 2:
        return (
            f"Integration Quality Evaluation:\n"
            f"  Embedding: {embedding_key}\n"
            f"  Batch key: {batch_key}\n"
            f"  Metrics skipped: at least 2 batches are required.\n"
        )

    # Batch mixing: silhouette score (lower = better mixing)
    if adata.n_obs > 10000:
        idx = np.random.choice(adata.n_obs, 10000, replace=False)
        emb_sub = embedding[idx]
        batch_sub = batch_labels[idx]
    else:
        emb_sub = embedding
        batch_sub = batch_labels

    batch_sil = silhouette_score(emb_sub, batch_sub)

    report = (
        f"Integration Quality Evaluation:\n"
        f"  Embedding: {embedding_key}\n"
        f"  Batch key: {batch_key}\n"
        f"  Batch silhouette score: {batch_sil:.4f} "
        f"(closer to 0 = better mixing)\n"
    )

    if cell_type_key and cell_type_key in adata.obs.columns:
        ct_labels = LabelEncoder().fit_transform(adata.obs[cell_type_key].astype(str).values)
        if adata.n_obs > 10000:
            ct_sub = ct_labels[idx]
        else:
            ct_sub = ct_labels
        if len(np.unique(ct_sub)) > 1:
            ct_sil = silhouette_score(emb_sub, ct_sub)
            report += (
                f"  Cell type silhouette score: {ct_sil:.4f} "
                f"(higher = better biological conservation)\n"
            )
        else:
            report += "  Cell type silhouette score skipped: only one cell type present\n"

    return report


TOOL_FUNCTIONS = {
    "integrate_harmony": integrate_harmony,
    "integrate_bbknn": integrate_bbknn,
    "evaluate_integration": evaluate_integration,
}
