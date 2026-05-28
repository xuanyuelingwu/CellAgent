"""CellAgent: AI Agent for Single-Cell Multi-Omics Analysis.

An autonomous AI agent specialized in single-cell RNA-seq data analysis,
inspired by Stanford's BioMNI framework. Features include:

- ReAct reasoning loop for autonomous analysis
- 39 specialized bioinformatics tools
- Domain-specific knowledge base
- LLM-driven resource retrieval
- Persistent code execution environment
- LLM-assisted cell type annotation
"""

__version__ = "0.1.0"
__author__ = "CellAgent Team"

from cellagent.config import CellAgentConfig

__all__ = ["CellAgent", "CellAgentConfig"]


def __getattr__(name):
    """Load the LLM-backed agent only when it is explicitly requested."""
    if name == "CellAgent":
        from cellagent.agent.cell_agent import CellAgent

        return CellAgent
    raise AttributeError(f"module 'cellagent' has no attribute {name!r}")
