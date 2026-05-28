"""Component tests for CellAgent.

Tests individual components: tool registry, knowledge loader,
code executor, and the full agent pipeline.
"""

import sys
import os
import importlib.util
import json
import tempfile

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _run_scanpy_tests() -> bool:
    """Gate heavier Scanpy/AnnData tests behind an explicit opt-in."""
    return os.environ.get("CELLAGENT_RUN_SCANPY_TESTS", "").lower() in {
        "1", "true", "yes"
    }


def test_tool_registry():
    """Test tool registry creation and lookup."""
    from cellagent.tools import build_default_registry

    registry = build_default_registry()

    # Check total tools registered
    n_tools = len(registry)
    print(f"[PASS] Tool registry: {n_tools} tools registered")
    assert n_tools >= 30, f"Expected >= 30 tools, got {n_tools}"

    # Check categories
    categories = registry.list_categories()
    expected_cats = ["preprocessing", "clustering", "annotation",
                     "differential_expression", "trajectory", "visualization",
                     "integration", "gene_analysis"]
    for cat in expected_cats:
        assert cat in categories, f"Missing category: {cat}"
    print(f"[PASS] All {len(expected_cats)} categories present: {categories}")

    # Check specific tools
    tool = registry.get_tool("load_h5ad")
    assert tool is not None
    assert tool.fn is not None
    print(f"[PASS] Tool lookup works: {tool.name} ({tool.category})")

    # Check prompt generation
    prompt = registry.get_prompt_description()
    assert len(prompt) > 500
    print(f"[PASS] Prompt description generated: {len(prompt)} chars")

    return True


def test_knowledge_loader():
    """Test knowledge base loading."""
    from cellagent.knowledge.loader import KnowledgeLoader

    loader = KnowledgeLoader()
    docs = loader.get_all_documents()

    print(f"[PASS] Knowledge loader: {len(docs)} documents loaded")
    assert len(docs) >= 4, f"Expected >= 4 docs, got {len(docs)}"

    # Check document structure
    for doc in docs:
        assert "id" in doc
        assert "name" in doc
        assert "description" in doc
        assert "content" in doc
        print(f"  - {doc['id']}: {doc['name'][:50]}")

    # Check specific document
    qc_doc = loader.get_document_by_id("qc_best_practices")
    assert qc_doc is not None
    assert "mitochondrial" in qc_doc["content"].lower()
    print(f"[PASS] QC best practices document loaded and contains expected content")

    guardrail_doc = loader.get_document_by_id("study_design_guardrails")
    assert guardrail_doc is not None
    assert "sample-level" in guardrail_doc["content"].lower()
    print(f"[PASS] Study design guardrails document loaded")

    return True


def test_code_executor():
    """Test persistent code execution."""
    from cellagent.agent.executor import CodeExecutor

    executor = CodeExecutor(output_dir="/tmp/cellagent_test")

    # Test basic execution
    result = executor.execute("x = 42\nprint(f'x = {x}')")
    assert result["success"], f"Execution failed: {result['error']}"
    assert "x = 42" in result["output"]
    print(f"[PASS] Basic code execution works")

    # Test persistence
    result = executor.execute("y = x * 2\nprint(f'y = {y}')")
    assert result["success"], f"Persistence failed: {result['error']}"
    assert "y = 84" in result["output"]
    print(f"[PASS] Variable persistence works")

    # Test error handling
    result = executor.execute("1/0")
    assert not result["success"]
    assert "ZeroDivisionError" in result["error"]
    print(f"[PASS] Error handling works")

    # Test variable listing
    variables = executor.list_variables()
    assert "x" in variables
    assert "y" in variables
    print(f"[PASS] Variable listing works: {variables}")

    # Test unsafe operation blocking
    result = executor.execute("import os\nos.system('echo unsafe')")
    assert not result["success"]
    assert "Blocked unsafe call" in result["error"]
    print(f"[PASS] Unsafe shell calls are blocked")

    # Test scanpy in executor
    if not _run_scanpy_tests():
        print("[SKIP] Scanpy/AnnData execution test skipped; set CELLAGENT_RUN_SCANPY_TESTS=1 to enable")
        return True
    if importlib.util.find_spec("scanpy") is None or importlib.util.find_spec("anndata") is None:
        print("[SKIP] Scanpy/AnnData execution test skipped; dependencies not installed")
        return True

    result = executor.execute("""
import scanpy as sc
import anndata as ad
import numpy as np

# Create a small test dataset
np.random.seed(42)
adata = ad.AnnData(
    X=np.random.poisson(1, size=(100, 200)).astype(float),
    obs={'batch': ['A']*50 + ['B']*50},
)
adata.var_names = [f'Gene_{i}' for i in range(200)]
adata.obs_names = [f'Cell_{i}' for i in range(100)]
print(f'Created test AnnData: {adata.shape}')
""")
    assert result["success"], f"Scanpy test failed: {result['error']}"
    print(f"[PASS] Scanpy execution in executor works")

    # Test AnnData summary
    summary = executor.get_adata_summary()
    assert summary is not None
    assert "100 cells" in summary
    print(f"[PASS] AnnData summary works")

    return True


def test_resource_retriever():
    """Test resource retrieval (mock without LLM)."""
    from cellagent.retriever.resource_retriever import ResourceRetriever

    retriever = ResourceRetriever()

    # Test prompt building
    resources = {
        "tools": [
            {"name": "load_h5ad", "description": "Load h5ad file"},
            {"name": "run_pca", "description": "Run PCA"},
        ],
        "knowledge": [
            {"id": "qc", "name": "QC Guide", "description": "QC best practices"},
        ],
        "libraries": [
            {"name": "scanpy", "description": "Single-cell analysis"},
        ],
    }

    prompt = retriever._build_prompt("Analyze PBMC data", resources)
    assert "PBMC" in prompt
    assert "load_h5ad" in prompt
    print(f"[PASS] Resource retriever prompt building works ({len(prompt)} chars)")

    # Test response parsing
    mock_response = "TOOLS: [0, 1]\nKNOWLEDGE: [0]\nLIBRARIES: [0]"
    indices = retriever._parse_response(mock_response)
    assert indices["tools"] == [0, 1]
    assert indices["knowledge"] == [0]
    print(f"[PASS] Response parsing works: {indices}")

    # Test invalid index filtering
    selected = retriever._select_resources(
        resources,
        {"tools": [-1, 0, 99], "knowledge": [0], "libraries": [0]},
    )
    assert [tool["name"] for tool in selected["tools"]] == ["load_h5ad"]
    print(f"[PASS] Invalid retrieval indices are ignored")

    return True


def test_provenance_recorder():
    """Test run manifest creation and redaction."""
    from cellagent.agent.provenance import ProvenanceRecorder

    with tempfile.TemporaryDirectory() as tmpdir:
        recorder = ProvenanceRecorder(tmpdir, run_id="test-run")
        selected_resources = {
            "tools": [{"name": "load_h5ad", "category": "preprocessing"}],
            "knowledge": [{"id": "qc_best_practices", "name": "QC Guide"}],
            "libraries": [{"name": "scanpy"}],
        }

        recorder.start_run(
            query="Analyze PBMC data",
            data_path="pbmc.h5ad",
            config={"llm_model": "test-model", "api_key": "secret"},
            selected_resources=selected_resources,
        )
        recorder.record_iteration(
            iteration=1,
            response="THOUGHT: test\nCODE:\n```python\nx = 1\n```",
            code="x = 1",
            execution_result={
                "success": True,
                "duration_seconds": 0.01,
                "output": "ok",
                "error": "",
                "variables": ["x"],
            },
        )

        output_file = os.path.join(recorder.run_dir, "plot.png")
        with open(output_file, "wb") as f:
            f.write(b"fake")

        manifest_path = recorder.finish_run(
            status="completed",
            final_answer="Done",
            final_adata_summary="AnnData: adata",
        )

        with open(manifest_path, encoding="utf-8") as f:
            manifest = json.load(f)

        assert manifest["run_id"] == "test-run"
        assert manifest["config"]["api_key"] == "***"
        assert manifest["selected_resources"]["tools"][0]["name"] == "load_h5ad"
        assert manifest["iterations"][0]["code"] == "x = 1"
        assert manifest["iterations"][0]["execution"]["duration_seconds"] == 0.01
        assert manifest["output_files"][0]["path"] == "plot.png"
        print(f"[PASS] Provenance manifest created: {manifest_path}")

    return True


def test_agent_run_directory():
    """Test CellAgent uses isolated run output directories."""
    from cellagent.config import CellAgentConfig
    from cellagent.agent.cell_agent import CellAgent

    with tempfile.TemporaryDirectory() as tmpdir:
        config = CellAgentConfig(verbose=False, use_resource_retriever=False)
        agent = CellAgent(config=config, output_dir=tmpdir, run_id="demo-run")

        assert agent.provenance.run_id == "demo-run"
        assert agent.output_dir == os.path.join(os.path.abspath(tmpdir), "demo-run")
        assert agent.executor.get_variable("OUTPUT_DIR") == agent.output_dir
        print(f"[PASS] Agent output directory isolated: {agent.output_dir}")

    return True


def test_synthetic_pipeline():
    """Test a synthetic analysis pipeline using tool functions directly."""
    if not _run_scanpy_tests():
        print("[SKIP] Synthetic pipeline skipped; set CELLAGENT_RUN_SCANPY_TESTS=1 to enable")
        return True
    if importlib.util.find_spec("scanpy") is None or importlib.util.find_spec("anndata") is None:
        print("[SKIP] Synthetic pipeline skipped; Scanpy/AnnData not installed")
        return True

    import numpy as np
    import anndata as ad
    import scanpy as sc

    from cellagent.tools.preprocessing import (
        calculate_qc_metrics, filter_cells_and_genes,
        normalize_and_log_transform, select_highly_variable_genes, run_pca
    )
    from cellagent.tools.clustering import (
        compute_neighbors, run_umap, run_leiden_clustering
    )
    from cellagent.tools.annotation import find_marker_genes

    output_dir = "/tmp/cellagent_test_pipeline"
    os.makedirs(output_dir, exist_ok=True)

    # Create synthetic dataset
    np.random.seed(42)
    n_cells, n_genes = 500, 1000

    # Create 3 distinct cell populations
    counts = np.zeros((n_cells, n_genes))
    # Population 1: high expression of genes 0-100
    counts[:200, :100] = np.random.poisson(5, size=(200, 100))
    # Population 2: high expression of genes 100-200
    counts[200:350, 100:200] = np.random.poisson(5, size=(150, 100))
    # Population 3: high expression of genes 200-300
    counts[350:, 200:300] = np.random.poisson(5, size=(150, 100))
    # Background expression
    counts += np.random.poisson(0.5, size=(n_cells, n_genes))

    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    # Add some "MT-" genes
    gene_names[990] = "MT-CO1"
    gene_names[991] = "MT-CO2"
    gene_names[992] = "MT-ND1"

    adata = ad.AnnData(
        X=counts.astype(float),
        obs={"batch": (["A"] * 250 + ["B"] * 250)},
    )
    adata.var_names = gene_names
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    print(f"Created synthetic dataset: {adata.shape}")

    # Step 1: QC
    result = calculate_qc_metrics(adata, output_dir=output_dir)
    assert "QC Metrics Summary" in result
    print(f"[PASS] QC metrics calculated")

    # Step 2: Filter
    result = filter_cells_and_genes(adata, min_genes=10, min_cells=1, max_pct_mt=50)
    assert "Filtering Results" in result
    print(f"[PASS] Filtering completed: {adata.shape}")

    # Step 3: Normalize
    result = normalize_and_log_transform(adata)
    assert "Normalization Complete" in result
    print(f"[PASS] Normalization completed")

    # Step 4: HVG
    result = select_highly_variable_genes(adata, n_top_genes=500, flavor="seurat")
    assert "Highly Variable Gene" in result
    print(f"[PASS] HVG selection completed")

    # Step 5: PCA
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    result = run_pca(adata_hvg, n_comps=30, output_dir=output_dir)
    assert "PCA Results" in result
    print(f"[PASS] PCA completed")

    # Copy PCA results back
    adata.obsm["X_pca"] = np.zeros((adata.n_obs, 30))
    adata.obsm["X_pca"][:adata_hvg.n_obs, :] = adata_hvg.obsm["X_pca"]

    # Step 6: Neighbors
    result = compute_neighbors(adata_hvg, n_neighbors=15, n_pcs=20)
    assert "Neighbor Graph" in result
    print(f"[PASS] Neighbor graph computed")

    # Step 7: Clustering
    result = run_leiden_clustering(adata_hvg, resolution=0.5)
    assert "Leiden Clustering" in result
    n_clusters = adata_hvg.obs["leiden"].nunique()
    print(f"[PASS] Leiden clustering: {n_clusters} clusters")

    # Step 8: UMAP
    result = run_umap(adata_hvg, output_dir=output_dir)
    assert "UMAP Computed" in result
    print(f"[PASS] UMAP computed")

    # Step 9: Marker genes
    result = find_marker_genes(adata_hvg, groupby="leiden", output_dir=output_dir)
    assert "Marker Gene" in result
    print(f"[PASS] Marker genes found")

    # Check output files
    output_files = os.listdir(output_dir)
    print(f"\nGenerated files: {output_files}")
    assert any("qc" in f for f in output_files), "Missing QC plots"
    assert any("umap" in f for f in output_files), "Missing UMAP plot"
    print(f"[PASS] All expected output files generated")

    return True


def run_all_tests():
    """Run all component tests."""
    print("=" * 60)
    print("CellAgent Component Tests")
    print("=" * 60)

    tests = [
        ("Tool Registry", test_tool_registry),
        ("Knowledge Loader", test_knowledge_loader),
        ("Code Executor", test_code_executor),
        ("Resource Retriever", test_resource_retriever),
        ("Provenance Recorder", test_provenance_recorder),
        ("Agent Run Directory", test_agent_run_directory),
        ("Synthetic Pipeline", test_synthetic_pipeline),
    ]

    results = []
    for name, test_fn in tests:
        print(f"\n--- Testing: {name} ---")
        try:
            success = test_fn()
            results.append((name, True, None))
        except Exception as e:
            import traceback
            results.append((name, False, str(e)))
            print(f"[FAIL] {name}: {e}")
            traceback.print_exc()

    print("\n" + "=" * 60)
    print("Test Summary:")
    print("=" * 60)
    for name, success, error in results:
        status = "PASS" if success else "FAIL"
        print(f"  [{status}] {name}" + (f" - {error}" if error else ""))

    n_passed = sum(1 for _, s, _ in results if s)
    n_total = len(results)
    print(f"\n{n_passed}/{n_total} tests passed")

    return n_passed == n_total


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
