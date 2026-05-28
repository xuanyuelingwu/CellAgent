"""CellAgent core agent module."""

from cellagent.agent.executor import CodeExecutor

__all__ = ["CellAgent", "CodeExecutor"]


def __getattr__(name):
    """Load the LLM-backed agent only when requested."""
    if name == "CellAgent":
        from cellagent.agent.cell_agent import CellAgent

        return CellAgent
    raise AttributeError(f"module 'cellagent.agent' has no attribute {name!r}")
