# CellAgent

CellAgent is an alpha-stage Python package for single-cell RNA-seq analysis. It combines a ReAct-style LLM agent, a persistent Python execution environment, a small curated knowledge base, and 39 Scanpy-centered analysis tools.

The project is useful as a prototype or teaching scaffold for LLM-assisted single-cell workflows. It is not yet a fully validated production analysis platform, and its multi-omics support is currently limited to guidance documents plus a few batch-integration tools.

## What It Does

- Runs an LLM-guided analysis loop that can plan, write Python code, execute it, observe results, and continue.
- Provides direct Python tools for common scRNA-seq workflows: loading, QC, filtering, normalization, HVG selection, PCA, neighbors, clustering, UMAP, marker discovery, annotation, differential expression, trajectory analysis, integration, and visualization.
- Keeps analysis state in a persistent execution namespace, so variables such as `adata` survive across ReAct steps.
- Loads five Markdown knowledge documents for QC, marker genes, workflow guidance, troubleshooting, and multi-omics concepts.
- Supports OpenAI-compatible chat completion APIs.
- Applies basic executor hardening: obvious shell/subprocess calls are blocked and Python-level execution timeout is enforced.

## Current Scope

Implemented well enough for experimentation:

- Standard scRNA-seq workflows using AnnData and Scanpy.
- PBMC-style marker-based annotation.
- LLM-assisted annotation from cluster marker genes.
- Harmony and BBKNN batch integration wrappers.
- Exploratory differential expression with `scanpy.tl.rank_genes_groups`.
- Markdown summary reports and common plots.

Not fully implemented yet:

- Full multi-omics pipelines for CITE-seq, scATAC-seq, spatial transcriptomics, or 10X Multiome.
- Reference mapping tools such as CellTypist, Azimuth, scArches, or scANVI.
- RNA velocity, ligand-receptor communication, MOFA+, ArchR, Signac, or CellChat/LIANA backends.
- Production-grade sandboxing for arbitrary generated code.
- Statistically rigorous multi-sample condition testing such as pseudobulk DE with donor-level replication.

## Installation

```bash
cd BioInformatics/CellAgent
pip install -e .
```

Install optional tools used by the full tool catalog:

```bash
pip install -e ".[full]"
```

For development:

```bash
pip install -e ".[dev,full]"
```

Core dependencies are declared in `pyproject.toml`: `openai`, `scanpy`, `anndata`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `scipy`, and `scikit-learn`.

Optional dependencies include `harmonypy`, `bbknn`, `scrublet`, `gseapy`, and `leidenalg`.

## Configuration

Environment variables:

| Variable | Purpose |
| --- | --- |
| `OPENAI_API_KEY` | API key for OpenAI-compatible LLM calls |
| `CELLAGENT_BASE_URL` | Optional custom OpenAI-compatible API base URL |
| `CELLAGENT_LLM_MODEL` | Default model override |
| `CELLAGENT_LLM_SOURCE` | Provider label, currently informational |
| `CELLAGENT_OUTPUT_PATH` | Default output path in config |
| `CELLAGENT_DATA_PATH` | Default data path in config |
| `CELLAGENT_TIMEOUT` | Python execution timeout in seconds |

Python configuration:

```python
from cellagent import CellAgent, CellAgentConfig

config = CellAgentConfig(
    llm_model="gpt-4.1-mini",
    temperature=0.1,
    max_iterations=20,
    timeout_seconds=600,
    use_resource_retriever=True,
    verbose=True,
)

agent = CellAgent(config=config, output_dir="./output")
```

## Quick Start: Agent Mode

```python
from cellagent import CellAgent, CellAgentConfig

config = CellAgentConfig(
    llm_model="gpt-4.1-mini",
    temperature=0.1,
    max_iterations=20,
    verbose=True,
)

agent = CellAgent(config=config, output_dir="./output/pbmc_agent")

result = agent.run(
    "Load the PBMC3k dataset from scanpy, run QC, filter low-quality cells, "
    "normalize, select highly variable genes, run PCA, neighbors, Leiden "
    "clustering, UMAP, marker discovery, and marker-based PBMC annotation."
)

print(result)
```

With a custom `.h5ad` file:

```python
result = agent.run(
    query="Load this human liver biopsy dataset, run QC, clustering, UMAP, and annotation.",
    data_path="./data/liver_sample.h5ad",
)
```

When `data_path` is provided, the path is also injected into the execution namespace as `DATA_PATH`. The agent is instructed to load it into an AnnData variable named `adata`.

## Quick Start: Direct Tool Mode

Direct tool mode is more deterministic than agent mode and is recommended for reproducible analysis scripts.

```python
import scanpy as sc

from cellagent.tools.preprocessing import (
    calculate_qc_metrics,
    filter_cells_and_genes,
    normalize_and_log_transform,
    select_highly_variable_genes,
    run_pca,
)
from cellagent.tools.clustering import (
    compute_neighbors,
    run_leiden_clustering,
    run_umap,
)
from cellagent.tools.annotation import find_marker_genes, annotate_with_markers
from cellagent.tools.visualization import plot_umap_colored, generate_analysis_summary

output_dir = "./output/direct_tools"

adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()

print(calculate_qc_metrics(adata, output_dir=output_dir))
print(filter_cells_and_genes(adata, min_genes=200, max_genes=2500, max_pct_mt=5))
print(normalize_and_log_transform(adata))
print(select_highly_variable_genes(adata, n_top_genes=2000, flavor="seurat"))

adata = adata[:, adata.var["highly_variable"]].copy()

print(run_pca(adata, n_comps=50, output_dir=output_dir))
print(compute_neighbors(adata, n_neighbors=10, n_pcs=40))
print(run_leiden_clustering(adata, resolution=0.5))
print(run_umap(adata, min_dist=0.3, output_dir=output_dir))
print(find_marker_genes(adata, groupby="leiden", output_dir=output_dir))

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

print(annotate_with_markers(adata, pbmc_markers, output_dir=output_dir))
print(plot_umap_colored(adata, ["leiden", "cell_type_marker"], output_dir=output_dir))
print(generate_analysis_summary(adata, output_dir=output_dir, title="PBMC3k Analysis"))
```

## Command Line

```bash
cellagent --query "Analyze PBMC3k dataset" --output ./results
cellagent --query "Run QC and clustering" --data ./data/sample.h5ad --output ./results
cellagent --interactive
cellagent --query "Cluster the cells" --model gpt-4.1-mini --max-iterations 20
```

Useful flags:

| Flag | Description |
| --- | --- |
| `--query`, `-q` | Single natural-language analysis request |
| `--data`, `-d` | Path to input `.h5ad` file |
| `--output`, `-o` | Output directory |
| `--model`, `-m` | LLM model name |
| `--interactive`, `-i` | Start chat-style interactive mode |
| `--no-retriever` | Disable LLM resource retrieval and expose all tools |
| `--max-iterations` | Maximum ReAct loop iterations |

## Architecture

```text
User query
   |
   v
ResourceRetriever
   | selects tools, knowledge, and libraries when available
   v
CellAgent ReAct loop
   | THINK -> CODE -> EXECUTE -> OBSERVE
   v
CodeExecutor
   | persistent Python namespace with tool functions injected
   v
AnnData objects, plots, CSV files, and Markdown reports
```

Important implementation details:

- `CellAgent` lives in `cellagent/agent/cell_agent.py`.
- `CodeExecutor` lives in `cellagent/agent/executor.py`.
- Tool schemas and registration live in `cellagent/tools/tool_registry.py` and `cellagent/tools/__init__.py`.
- Knowledge loading lives in `cellagent/knowledge/loader.py`.
- Resource selection lives in `cellagent/retriever/resource_retriever.py`.
- LLM client creation lives in `cellagent/llm.py`.

## Tool Catalog

CellAgent currently registers 39 tools across 8 categories.

| Category | Tools |
| --- | --- |
| Preprocessing | `load_h5ad`, `load_10x_mtx`, `calculate_qc_metrics`, `filter_cells_and_genes`, `detect_doublets`, `normalize_and_log_transform`, `select_highly_variable_genes`, `run_pca` |
| Clustering | `compute_neighbors`, `run_umap`, `run_tsne`, `run_leiden_clustering`, `run_louvain_clustering`, `evaluate_clustering_quality` |
| Annotation | `find_marker_genes`, `annotate_with_markers`, `annotate_with_llm`, `score_cell_types`, `compare_annotations` |
| Differential expression | `differential_expression`, `volcano_plot`, `gene_set_enrichment`, `gsea_prerank` |
| Trajectory | `run_diffusion_map`, `compute_pseudotime_dpt`, `run_paga`, `identify_trajectory_genes` |
| Visualization | `plot_umap_colored`, `plot_dotplot`, `plot_heatmap`, `plot_stacked_violin`, `plot_cell_composition`, `generate_analysis_summary` |
| Integration | `integrate_harmony`, `integrate_bbknn`, `evaluate_integration` |
| Gene analysis | `gene_correlation_network`, `cell_cycle_scoring`, `query_gene_info` |

## Knowledge Base

The bundled knowledge documents are Markdown files under `cellagent/knowledge/`.

| File | Purpose |
| --- | --- |
| `qc_best_practices.md` | QC thresholds and filtering guidance |
| `cell_markers.md` | Canonical marker genes for common tissues and species |
| `workflow_guide.md` | Step-by-step scRNA-seq workflow recommendations |
| `troubleshooting.md` | Common problems and suggested fixes |
| `multiomics_integration.md` | Conceptual guidance for multi-sample and multi-omics integration |

Some knowledge documents mention methods that are not implemented as tools yet, such as scVI, CellTypist, RNA velocity, MOFA+, ArchR, Signac, CellChat, and LIANA. Treat those sections as expert guidance for future extension, not as a guarantee of current built-in functionality.

## Project Layout

```text
CellAgent/
├── .gitignore
├── LICENSE
├── README.md
├── pyproject.toml
├── cellagent/
│   ├── __init__.py
│   ├── cli.py
│   ├── config.py
│   ├── llm.py
│   ├── agent/
│   │   ├── __init__.py
│   │   ├── cell_agent.py
│   │   └── executor.py
│   ├── knowledge/
│   │   ├── __init__.py
│   │   ├── loader.py
│   │   ├── cell_markers.md
│   │   ├── multiomics_integration.md
│   │   ├── qc_best_practices.md
│   │   ├── troubleshooting.md
│   │   └── workflow_guide.md
│   ├── retriever/
│   │   ├── __init__.py
│   │   └── resource_retriever.py
│   └── tools/
│       ├── __init__.py
│       ├── annotation.py
│       ├── clustering.py
│       ├── differential.py
│       ├── gene_analysis.py
│       ├── integration.py
│       ├── preprocessing.py
│       ├── tool_registry.py
│       ├── trajectory.py
│       └── visualization.py
├── examples/
│   ├── quickstart.py
│   ├── direct_tools.py
│   ├── test_report.md
│   └── demo_output/
│       ├── cell_composition.png
│       ├── deg_results.csv
│       ├── dotplot.png
│       ├── marker_dotplot.png
│       ├── marker_genes.png
│       ├── pca_elbow.png
│       ├── qc_scatter.png
│       ├── qc_violin.png
│       ├── umap.png
│       ├── umap_annotated.png
│       └── volcano_0.png
└── tests/
    └── test_components.py
```

Generated directories such as `__pycache__/` and `.pyc` files are intentionally omitted from the layout and should not be committed. They are already covered by `.gitignore`.

## Examples And Demo Artifacts

- `examples/quickstart.py` demonstrates agent mode.
- `examples/direct_tools.py` demonstrates deterministic direct tool usage.
- `examples/test_report.md` is a historical PBMC3k end-to-end report snapshot.
- `examples/demo_output/` contains static images and a CSV used by that report.

The report and demo outputs are useful for documentation, but they are not automatically regenerated by the test suite. If tool behavior changes, rerun the example and refresh those artifacts before treating them as current benchmark results.

## Testing

Run the component tests:

```bash
python tests/test_components.py
```

Or, if using pytest:

```bash
pytest tests
```

The tests cover:

- Tool registry construction.
- Knowledge document loading.
- Persistent code execution and unsafe-call blocking.
- Resource retriever prompt parsing and invalid-index filtering.
- A synthetic Scanpy pipeline when explicitly enabled.

By default, the heavier Scanpy/AnnData execution checks are skipped so that core component tests stay fast and stable across CI environments. Enable them with:

```bash
CELLAGENT_RUN_SCANPY_TESTS=1 python tests/test_components.py
```

## Known Limitations

- Agent mode depends on LLM output quality and can still produce invalid Python.
- The execution environment is hardened but not a secure sandbox for untrusted code.
- Many tools mutate AnnData objects in place.
- Output filenames are mostly fixed, so repeated runs can overwrite files in the same output directory.
- Differential expression is exploratory and cluster-centric; it should not replace careful multi-sample statistical analysis.
- The CLI currently defaults to verbose output.
- The knowledge base is small and static.

## Development Notes

- Keep tool metadata and implementations synchronized. Every tool should appear in `TOOL_DESCRIPTIONS`, `TOOL_FUNCTIONS`, and the default registry.
- Prefer adding precondition checks that return clear error messages over letting Scanpy exceptions surface without context.
- When adding optional dependencies, update both `pyproject.toml` and this README.
- If a knowledge document mentions a method that is not implemented, label it as guidance rather than a built-in feature.
- Refresh `examples/test_report.md` and `examples/demo_output/` only after rerunning the corresponding example.

## License

MIT License. See `LICENSE`.

## Acknowledgments

CellAgent is inspired by the BioMNI agent architecture and builds on the Python single-cell ecosystem, especially Scanpy and AnnData.
