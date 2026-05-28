"""CellAgent tool system.

Provides a unified registry of all available single-cell analysis tools
with declarative metadata and implementations.
"""

from cellagent.tools.tool_registry import ToolRegistry, ToolSchema

from cellagent.tools import preprocessing
from cellagent.tools import clustering
from cellagent.tools import annotation
from cellagent.tools import differential
from cellagent.tools import trajectory
from cellagent.tools import visualization
from cellagent.tools import integration
from cellagent.tools import gene_analysis


def build_default_registry() -> ToolRegistry:
    """Build and return the default tool registry with all tools registered."""
    registry = ToolRegistry()

    modules = [
        (preprocessing.TOOL_DESCRIPTIONS, preprocessing.TOOL_FUNCTIONS),
        (clustering.TOOL_DESCRIPTIONS, clustering.TOOL_FUNCTIONS),
        (annotation.TOOL_DESCRIPTIONS, annotation.TOOL_FUNCTIONS),
        (differential.TOOL_DESCRIPTIONS, differential.TOOL_FUNCTIONS),
        (trajectory.TOOL_DESCRIPTIONS, trajectory.TOOL_FUNCTIONS),
        (visualization.TOOL_DESCRIPTIONS, visualization.TOOL_FUNCTIONS),
        (integration.TOOL_DESCRIPTIONS, integration.TOOL_FUNCTIONS),
        (gene_analysis.TOOL_DESCRIPTIONS, gene_analysis.TOOL_FUNCTIONS),
    ]

    for descriptions, functions in modules:
        for desc in descriptions:
            fn = functions.get(desc["name"])
            registry.register_from_dict(desc, fn=fn)

    return registry


ALL_TOOL_FUNCTIONS = {}
for _mod in [preprocessing, clustering, annotation, differential,
             trajectory, visualization, integration, gene_analysis]:
    ALL_TOOL_FUNCTIONS.update(_mod.TOOL_FUNCTIONS)
