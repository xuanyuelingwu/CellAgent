"""Preprocessing tools for single-cell RNA-seq data.

Implements quality control, filtering, normalization, and feature selection
steps following best practices from the single-cell community.
"""

# Tool descriptions (declarative metadata, following BioMNI pattern)
TOOL_DESCRIPTIONS = [
    {
        "name": "load_h5ad",
        "description": (
            "Load a single-cell dataset from an h5ad file into an AnnData object. "
            "Supports .h5ad format (Scanpy/AnnData native format). "
            "Returns the loaded AnnData object."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "file_path", "type": "str", "description": "Path to the .h5ad file"},
        ],
        "optional_parameters": [],
        "returns": "AnnData",
        "example": 'adata = load_h5ad("/data/pbmc3k.h5ad")',
    },
    {
        "name": "load_10x_mtx",
        "description": (
            "Load single-cell data from 10X Genomics Cell Ranger output directory "
            "containing matrix.mtx, barcodes.tsv, and features/genes.tsv files. "
            "Returns an AnnData object with the expression matrix."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "data_dir", "type": "str", "description": "Path to the 10X output directory"},
        ],
        "optional_parameters": [
            {"name": "var_names", "type": "str", "default": "gene_symbols",
             "description": "Use 'gene_symbols' or 'gene_ids' for variable names"},
        ],
        "returns": "AnnData",
        "example": 'adata = load_10x_mtx("/data/filtered_feature_bc_matrix/")',
    },
    {
        "name": "calculate_qc_metrics",
        "description": (
            "Calculate quality control metrics for single-cell data including: "
            "number of genes per cell (n_genes_by_counts), total UMI counts per cell "
            "(total_counts), percentage of mitochondrial genes (pct_counts_mt), "
            "and percentage of ribosomal genes (pct_counts_ribo). "
            "Generates violin plots and scatter plots for QC visualization."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object to analyze"},
        ],
        "optional_parameters": [
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save QC plots"},
        ],
        "returns": "str",
        "example": 'qc_report = calculate_qc_metrics(adata, output_dir="./qc_plots")',
    },
    {
        "name": "filter_cells_and_genes",
        "description": (
            "Filter low-quality cells and genes based on QC thresholds. "
            "Removes cells with too few/many genes, extreme UMI counts, "
            "or high mitochondrial content. Also removes genes expressed in too few cells."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object to filter"},
        ],
        "optional_parameters": [
            {"name": "min_genes", "type": "int", "default": 200,
             "description": "Minimum number of genes per cell"},
            {"name": "max_genes", "type": "int", "default": 5000,
             "description": "Maximum number of genes per cell"},
            {"name": "min_cells", "type": "int", "default": 3,
             "description": "Minimum number of cells per gene"},
            {"name": "max_pct_mt", "type": "float", "default": 20.0,
             "description": "Maximum percentage of mitochondrial genes"},
            {"name": "max_counts", "type": "int", "default": None,
             "description": "Maximum total counts per cell"},
        ],
        "returns": "str",
        "example": 'result = filter_cells_and_genes(adata, min_genes=200, max_pct_mt=15)',
    },
    {
        "name": "detect_doublets",
        "description": (
            "Detect and optionally remove doublets (cell multiplets) using Scrublet. "
            "Doublets are artifacts where two cells are captured in a single droplet, "
            "leading to mixed transcriptional profiles."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
        ],
        "optional_parameters": [
            {"name": "expected_doublet_rate", "type": "float", "default": 0.06,
             "description": "Expected doublet rate (typically 0.8% per 1000 cells)"},
            {"name": "remove", "type": "bool", "default": False,
             "description": "Whether to remove detected doublets"},
        ],
        "returns": "str",
        "example": 'result = detect_doublets(adata, expected_doublet_rate=0.06, remove=True)',
    },
    {
        "name": "normalize_and_log_transform",
        "description": (
            "Normalize total counts per cell to a target sum (default 1e4) and "
            "apply log1p transformation. This is the standard normalization approach "
            "for scRNA-seq data. Optionally stores raw counts in adata.raw."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object to normalize"},
        ],
        "optional_parameters": [
            {"name": "target_sum", "type": "float", "default": 1e4,
             "description": "Target sum for normalization"},
            {"name": "save_raw", "type": "bool", "default": True,
             "description": "Whether to save raw counts in adata.raw"},
        ],
        "returns": "str",
        "example": 'result = normalize_and_log_transform(adata, target_sum=1e4)',
    },
    {
        "name": "select_highly_variable_genes",
        "description": (
            "Identify highly variable genes (HVGs) for downstream analysis. "
            "Uses the Seurat v3 method by default. HVGs capture the most "
            "biologically meaningful variation in the data."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
        ],
        "optional_parameters": [
            {"name": "n_top_genes", "type": "int", "default": 2000,
             "description": "Number of highly variable genes to select"},
            {"name": "flavor", "type": "str", "default": "seurat_v3",
             "description": "Method for HVG selection: 'seurat', 'seurat_v3', or 'cell_ranger'"},
            {"name": "batch_key", "type": "str", "default": None,
             "description": "Batch key for batch-aware HVG selection"},
        ],
        "returns": "str",
        "example": 'result = select_highly_variable_genes(adata, n_top_genes=2000)',
    },
    {
        "name": "run_pca",
        "description": (
            "Perform Principal Component Analysis (PCA) on the data. "
            "Scales the data, computes PCA, and generates an elbow plot "
            "to help determine the optimal number of PCs."
        ),
        "module": "cellagent.tools.preprocessing",
        "category": "preprocessing",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
        ],
        "optional_parameters": [
            {"name": "n_comps", "type": "int", "default": 50,
             "description": "Number of principal components to compute"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save PCA plots"},
        ],
        "returns": "str",
        "example": 'result = run_pca(adata, n_comps=50)',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def _ensure_qc_metrics(adata) -> None:
    """Compute QC metrics when downstream filters need them."""
    import scanpy as sc

    if "n_genes_by_counts" in adata.obs and "total_counts" in adata.obs:
        return

    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def load_h5ad(file_path: str):
    """Load an h5ad file and return an AnnData object."""
    import scanpy as sc

    adata = sc.read_h5ad(file_path)
    adata.uns["cellagent_source"] = file_path
    return adata


def load_10x_mtx(data_dir: str, var_names: str = "gene_symbols"):
    """Load 10X Genomics data and return an AnnData object."""
    import scanpy as sc

    adata = sc.read_10x_mtx(data_dir, var_names=var_names, cache=True)
    adata.var_names_make_unique()
    adata.uns["cellagent_source"] = data_dir
    return adata


def calculate_qc_metrics(adata, output_dir: str = "./output") -> str:
    """Calculate QC metrics and generate plots."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")

    os.makedirs(output_dir, exist_ok=True)

    # Annotate mitochondrial and ribosomal genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )

    # Generate QC violin plots
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    import matplotlib.pyplot as plt
    plt.savefig(os.path.join(output_dir, "qc_violin.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Generate scatter plots
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", show=False)
    plt.savefig(os.path.join(output_dir, "qc_scatter.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Compute summary statistics
    report = (
        f"QC Metrics Summary:\n"
        f"  Total cells: {adata.n_obs}\n"
        f"  Total genes: {adata.n_vars}\n"
        f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}\n"
        f"  Median UMI/cell: {adata.obs['total_counts'].median():.0f}\n"
        f"  Median %MT: {adata.obs['pct_counts_mt'].median():.2f}%\n"
        f"  Cells with >20% MT: {(adata.obs['pct_counts_mt'] > 20).sum()}\n"
        f"  QC plots saved to: {output_dir}/\n"
    )
    return report


def filter_cells_and_genes(
    adata,
    min_genes: int = 200,
    max_genes: int = 5000,
    min_cells: int = 3,
    max_pct_mt: float = 20.0,
    max_counts: int | None = None,
) -> str:
    """Filter cells and genes based on QC thresholds."""
    import scanpy as sc

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    _ensure_qc_metrics(adata)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if max_genes is not None:
        adata._inplace_subset_obs(adata.obs["n_genes_by_counts"] < max_genes)

    if "pct_counts_mt" in adata.obs.columns:
        adata._inplace_subset_obs(adata.obs["pct_counts_mt"] < max_pct_mt)

    if max_counts is not None:
        adata._inplace_subset_obs(adata.obs["total_counts"] < max_counts)

    report = (
        f"Filtering Results:\n"
        f"  Cells: {n_cells_before} -> {adata.n_obs} "
        f"(removed {n_cells_before - adata.n_obs})\n"
        f"  Genes: {n_genes_before} -> {adata.n_vars} "
        f"(removed {n_genes_before - adata.n_vars})\n"
        f"  Thresholds: min_genes={min_genes}, max_genes={max_genes}, "
        f"min_cells={min_cells}, max_pct_mt={max_pct_mt}\n"
    )
    return report


def detect_doublets(
    adata, expected_doublet_rate: float = 0.06, remove: bool = False
) -> str:
    """Detect doublets using Scrublet."""
    import scrublet as scr
    import numpy as np

    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    n_doublets = predicted_doublets.sum()
    n_cells_before = adata.n_obs

    if remove:
        adata._inplace_subset_obs(~adata.obs["predicted_doublet"])

    report = (
        f"Doublet Detection Results:\n"
        f"  Detected doublets: {n_doublets} ({n_doublets/n_cells_before*100:.1f}%)\n"
        f"  Expected rate: {expected_doublet_rate*100:.1f}%\n"
        f"  Threshold: {scrub.threshold_:.3f}\n"
    )
    if remove:
        report += f"  Removed doublets: cells {n_cells_before} -> {adata.n_obs}\n"
    else:
        report += "  Doublets marked in adata.obs['predicted_doublet'] (not removed)\n"
    return report


def normalize_and_log_transform(
    adata, target_sum: float = 1e4, save_raw: bool = True
) -> str:
    """Normalize and log-transform the data."""
    import scanpy as sc

    if save_raw:
        adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    report = (
        f"Normalization Complete:\n"
        f"  Method: Total-count normalization + log1p\n"
        f"  Target sum: {target_sum:.0f}\n"
        f"  Raw data saved: {save_raw}\n"
        f"  Shape: {adata.shape}\n"
    )
    return report


def select_highly_variable_genes(
    adata,
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    batch_key: str | None = None,
) -> str:
    """Select highly variable genes."""
    import scanpy as sc

    requested_flavor = flavor

    try:
        # seurat_v3 expects counts. If raw counts were saved before
        # normalization, run HVG selection on that count matrix and map the
        # selected genes back onto the working AnnData object.
        if flavor == "seurat_v3" and adata.raw is not None:
            raw_adata = adata.raw.to_adata()
            sc.pp.highly_variable_genes(
                raw_adata, n_top_genes=n_top_genes, flavor=flavor,
                batch_key=batch_key
            )
            selected = set(raw_adata.var_names[raw_adata.var["highly_variable"]])
            adata.var["highly_variable"] = adata.var_names.isin(selected)
        else:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top_genes, flavor=flavor,
                batch_key=batch_key
            )
    except ImportError:
        flavor = "seurat"
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top_genes, flavor="seurat",
            batch_key=batch_key
        )

    n_hvg = adata.var["highly_variable"].sum()
    report = (
        f"Highly Variable Gene Selection:\n"
        f"  Method: {flavor}"
        f"{' (fallback from seurat_v3)' if requested_flavor != flavor else ''}\n"
        f"  HVGs selected: {n_hvg}\n"
        f"  Total genes: {adata.n_vars}\n"
    )
    if batch_key:
        report += f"  Batch-aware: Yes (key={batch_key})\n"
    return report


def run_pca(adata, n_comps: int = 50, output_dir: str = "./output") -> str:
    """Run PCA and generate elbow plot."""
    import scanpy as sc
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    max_comps = min(adata.n_obs, adata.n_vars) - 1
    if max_comps < 1:
        return "Error: PCA requires at least 2 cells and 2 genes after filtering."
    if n_comps > max_comps:
        n_comps = max_comps

    # Scale data
    sc.pp.scale(adata, max_value=10)

    # Run PCA
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")

    # Elbow plot
    sc.pl.pca_variance_ratio(adata, n_pcs=n_comps, show=False)
    plt.savefig(os.path.join(output_dir, "pca_elbow.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Determine suggested n_pcs (elbow heuristic)
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    cumvar = variance_ratio.cumsum()
    suggested_pcs = int((cumvar < 0.9).sum()) + 1
    suggested_pcs = min(suggested_pcs, n_comps)
    top_n = min(30, len(cumvar))

    report = (
        f"PCA Results:\n"
        f"  Components computed: {n_comps}\n"
        f"  Variance explained (top 10): "
        f"{', '.join(f'{v:.3f}' for v in variance_ratio[:10])}\n"
        f"  Cumulative variance (top {top_n}): {cumvar[top_n - 1]:.3f}\n"
        f"  Suggested n_pcs (90% variance): {suggested_pcs}\n"
        f"  Elbow plot saved to: {output_dir}/pca_elbow.png\n"
    )
    return report


# Map of function names to implementations
TOOL_FUNCTIONS = {
    "load_h5ad": load_h5ad,
    "load_10x_mtx": load_10x_mtx,
    "calculate_qc_metrics": calculate_qc_metrics,
    "filter_cells_and_genes": filter_cells_and_genes,
    "detect_doublets": detect_doublets,
    "normalize_and_log_transform": normalize_and_log_transform,
    "select_highly_variable_genes": select_highly_variable_genes,
    "run_pca": run_pca,
}
