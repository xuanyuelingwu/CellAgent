"""Gene-level analysis tools for single-cell data.

Implements gene correlation analysis, cell-cycle scoring, and LLM-assisted
gene information lookup.
"""

TOOL_DESCRIPTIONS = [
    {
        "name": "gene_correlation_network",
        "description": (
            "Compute gene-gene correlation network for a set of genes. "
            "Identifies co-expressed gene modules and hub genes. "
            "Uses Spearman correlation on the expression matrix."
        ),
        "module": "cellagent.tools.gene_analysis",
        "category": "gene_analysis",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
            {"name": "genes", "type": "list", "description": "List of genes to analyze"},
        ],
        "optional_parameters": [
            {"name": "threshold", "type": "float", "default": 0.3,
             "description": "Correlation threshold for network edges"},
            {"name": "output_dir", "type": "str", "default": "./output",
             "description": "Directory to save network plots"},
        ],
        "returns": "str",
        "example": 'result = gene_correlation_network(adata, ["CD3D", "CD3E", "CD4", "CD8A"])',
    },
    {
        "name": "cell_cycle_scoring",
        "description": (
            "Score cells for cell cycle phase (G1, S, G2M) using known "
            "cell cycle gene signatures. Assigns each cell a predicted phase."
        ),
        "module": "cellagent.tools.gene_analysis",
        "category": "gene_analysis",
        "required_parameters": [
            {"name": "adata", "type": "AnnData", "description": "AnnData object"},
        ],
        "optional_parameters": [
            {"name": "organism", "type": "str", "default": "human",
             "description": "Organism: 'human' or 'mouse'"},
        ],
        "returns": "str",
        "example": 'result = cell_cycle_scoring(adata, organism="human")',
    },
    {
        "name": "query_gene_info",
        "description": (
            "Query information about a gene using the LLM's knowledge base. "
            "Returns gene function, associated pathways, disease associations, "
            "and known cell type markers."
        ),
        "module": "cellagent.tools.gene_analysis",
        "category": "gene_analysis",
        "required_parameters": [
            {"name": "gene_name", "type": "str", "description": "Gene symbol to query"},
        ],
        "optional_parameters": [
            {"name": "context", "type": "str", "default": None,
             "description": "Additional context (e.g., tissue type, disease)"},
        ],
        "returns": "str",
        "example": 'result = query_gene_info("TP53", context="breast cancer")',
    },
]


# ============================================================
# Tool implementations
# ============================================================

def gene_correlation_network(adata, genes: list, threshold: float = 0.3,
                              output_dir: str = "./output") -> str:
    """Compute gene correlation network."""
    import numpy as np
    import pandas as pd
    from scipy.stats import spearmanr
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    valid_genes = [g for g in genes if g in adata.var_names]
    if len(valid_genes) < 2:
        return "Error: Need at least 2 valid genes for correlation network."

    # Extract expression
    if adata.raw is not None:
        expr = adata.raw[:, valid_genes].X
    else:
        expr = adata[:, valid_genes].X

    if hasattr(expr, "toarray"):
        expr = expr.toarray()

    # Compute correlation matrix
    corr_matrix, p_matrix = spearmanr(expr)
    if len(valid_genes) == 2:
        corr_matrix = np.array([[1, corr_matrix], [corr_matrix, 1]])

    # Create correlation heatmap
    fig, ax = plt.subplots(figsize=(max(8, len(valid_genes)), max(6, len(valid_genes) * 0.8)))
    im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(valid_genes)))
    ax.set_yticks(range(len(valid_genes)))
    ax.set_xticklabels(valid_genes, rotation=45, ha="right")
    ax.set_yticklabels(valid_genes)
    plt.colorbar(im, label="Spearman Correlation")
    ax.set_title("Gene Correlation Network")
    plt.tight_layout()

    path = os.path.join(output_dir, "gene_correlation.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    # Find strong correlations
    edges = []
    for i in range(len(valid_genes)):
        for j in range(i + 1, len(valid_genes)):
            if abs(corr_matrix[i, j]) > threshold:
                edges.append((valid_genes[i], valid_genes[j], corr_matrix[i, j]))

    report = f"Gene Correlation Network:\n"
    report += f"  Genes analyzed: {len(valid_genes)}\n"
    report += f"  Strong correlations (|r|>{threshold}): {len(edges)}\n"
    for g1, g2, r in sorted(edges, key=lambda x: abs(x[2]), reverse=True)[:20]:
        report += f"    {g1} <-> {g2}: r={r:.4f}\n"
    report += f"  Correlation plot saved to: {path}\n"
    return report


def cell_cycle_scoring(adata, organism: str = "human") -> str:
    """Score cells for cell cycle phase."""
    import scanpy as sc

    # Cell cycle genes (Tirosh et al. 2016)
    s_genes_human = [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
        "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
        "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
        "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
        "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
        "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
    ]
    g2m_genes_human = [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
        "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
        "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
        "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
        "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
        "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
        "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
        "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
    ]

    if organism.lower() == "mouse":
        s_genes = [g[0] + g[1:].lower() for g in s_genes_human]
        g2m_genes = [g[0] + g[1:].lower() for g in g2m_genes_human]
    else:
        s_genes = s_genes_human
        g2m_genes = g2m_genes_human

    # Filter to genes present in data
    s_valid = [g for g in s_genes if g in adata.var_names]
    g2m_valid = [g for g in g2m_genes if g in adata.var_names]

    if not s_valid or not g2m_valid:
        return (
            "Error: Not enough cell cycle genes found in the dataset. "
            f"S genes: {len(s_valid)}/{len(s_genes)}, "
            f"G2M genes: {len(g2m_valid)}/{len(g2m_genes)}. "
            "Check organism and gene symbol capitalization."
        )

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_valid, g2m_genes=g2m_valid)

    phase_counts = adata.obs["phase"].value_counts()
    report = (
        f"Cell Cycle Scoring:\n"
        f"  S-phase genes used: {len(s_valid)}/{len(s_genes)}\n"
        f"  G2M-phase genes used: {len(g2m_valid)}/{len(g2m_genes)}\n"
        f"  Phase distribution:\n"
    )
    for phase, count in phase_counts.items():
        report += f"    {phase}: {count} cells ({count/adata.n_obs*100:.1f}%)\n"
    report += "  Results stored in: adata.obs['phase'], adata.obs['S_score'], adata.obs['G2M_score']\n"
    return report


def query_gene_info(gene_name: str, context: str = None) -> str:
    """Query gene information using LLM."""
    from cellagent.llm import llm_chat

    prompt = f"""You are an expert molecular biologist. Provide a comprehensive summary of the gene {gene_name}.

Include:
1. Full gene name and aliases
2. Gene function and biological role
3. Associated pathways (KEGG, Reactome)
4. Known disease associations
5. Cell type specificity (which cell types express this gene)
6. Relevance as a marker gene in single-cell analysis
"""
    if context:
        prompt += f"\nAdditional context: {context}\n"

    prompt += "\nProvide a concise but informative response."

    response = llm_chat(
        messages=[{"role": "user", "content": prompt}],
        temperature=0.1,
        max_tokens=2048,
    )
    return f"Gene Information for {gene_name}:\n{response}"


TOOL_FUNCTIONS = {
    "gene_correlation_network": gene_correlation_network,
    "cell_cycle_scoring": cell_cycle_scoring,
    "query_gene_info": query_gene_info,
}
