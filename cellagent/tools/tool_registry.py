"""Tool Registry for CellAgent.

Manages registration, lookup, and metadata of all available analysis tools.
Following BioMNI's pattern of declarative tool descriptions separated from
implementation, enabling dynamic tool discovery and LLM-driven selection.
"""

from __future__ import annotations

import json
from typing import Any, Callable


class ToolSchema:
    """Schema describing a single tool's interface."""

    def __init__(
        self,
        name: str,
        description: str,
        module: str,
        category: str,
        required_parameters: list[dict],
        optional_parameters: list[dict] | None = None,
        returns: str = "str",
        example: str | None = None,
        fn: Callable | None = None,
    ):
        self.name = name
        self.description = description
        self.module = module
        self.category = category
        self.required_parameters = required_parameters
        self.optional_parameters = optional_parameters or []
        self.returns = returns
        self.example = example
        self.fn = fn

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "description": self.description,
            "module": self.module,
            "category": self.category,
            "required_parameters": self.required_parameters,
            "optional_parameters": self.optional_parameters,
            "returns": self.returns,
            "example": self.example,
        }

    def to_prompt_string(self) -> str:
        """Format tool info for inclusion in LLM prompts."""
        lines = [
            f"Function: {self.name}",
            f"  Module: {self.module}",
            f"  Description: {self.description}",
        ]
        if self.required_parameters:
            lines.append("  Required Parameters:")
            for p in self.required_parameters:
                lines.append(
                    f"    - {p['name']} ({p.get('type', 'Any')}): {p.get('description', '')}"
                )
        if self.optional_parameters:
            lines.append("  Optional Parameters:")
            for p in self.optional_parameters:
                default = p.get("default", "None")
                lines.append(
                    f"    - {p['name']} ({p.get('type', 'Any')}, default={default}): "
                    f"{p.get('description', '')}"
                )
        if self.example:
            lines.append(f"  Example: {self.example}")
        return "\n".join(lines)


class ToolRegistry:
    """Central registry for all CellAgent tools."""

    def __init__(self):
        self._tools: dict[str, ToolSchema] = {}
        self._categories: dict[str, list[str]] = {}

    def register(self, schema: ToolSchema) -> None:
        """Register a tool schema."""
        self._tools[schema.name] = schema
        cat = schema.category
        if cat not in self._categories:
            self._categories[cat] = []
        if schema.name not in self._categories[cat]:
            self._categories[cat].append(schema.name)

    def register_from_dict(self, d: dict, fn: Callable | None = None) -> None:
        """Register a tool from a dictionary description."""
        schema = ToolSchema(
            name=d["name"],
            description=d["description"],
            module=d.get("module", ""),
            category=d.get("category", "general"),
            required_parameters=d.get("required_parameters", []),
            optional_parameters=d.get("optional_parameters", []),
            returns=d.get("returns", "str"),
            example=d.get("example"),
            fn=fn,
        )
        self.register(schema)

    def get_tool(self, name: str) -> ToolSchema | None:
        """Get a tool by name."""
        return self._tools.get(name)

    def get_tools_by_category(self, category: str) -> list[ToolSchema]:
        """Get all tools in a category."""
        names = self._categories.get(category, [])
        return [self._tools[n] for n in names if n in self._tools]

    def list_tools(self) -> list[str]:
        """List all registered tool names."""
        return list(self._tools.keys())

    def list_categories(self) -> list[str]:
        """List all categories."""
        return list(self._categories.keys())

    def get_all_tools(self) -> list[ToolSchema]:
        """Get all registered tools."""
        return list(self._tools.values())

    def get_prompt_description(self, tools: list[str] | None = None) -> str:
        """Generate a formatted description of tools for LLM prompts.

        Args:
            tools: Optional list of tool names to include. If None, includes all.

        Returns:
            Formatted string describing the tools.
        """
        if tools is None:
            schemas = self.get_all_tools()
        else:
            schemas = [self._tools[n] for n in tools if n in self._tools]

        # Group by category
        by_category: dict[str, list[ToolSchema]] = {}
        for s in schemas:
            if s.category not in by_category:
                by_category[s.category] = []
            by_category[s.category].append(s)

        lines = []
        for cat in sorted(by_category.keys()):
            lines.append(f"\n=== {cat.upper()} ===")
            for s in by_category[cat]:
                lines.append(s.to_prompt_string())
                lines.append("")

        return "\n".join(lines)

    def get_summary_for_retrieval(self) -> list[dict]:
        """Get a list of tool summaries for the resource retriever."""
        return [
            {
                "name": s.name,
                "description": s.description,
                "category": s.category,
                "module": s.module,
            }
            for s in self._tools.values()
        ]

    def remove_tool(self, name: str) -> bool:
        """Remove a tool by name."""
        if name in self._tools:
            cat = self._tools[name].category
            del self._tools[name]
            if cat in self._categories and name in self._categories[cat]:
                self._categories[cat].remove(name)
            return True
        return False

    def __len__(self) -> int:
        return len(self._tools)

    def __contains__(self, name: str) -> bool:
        return name in self._tools
