"""Resource retriever for CellAgent.

Uses LLM to dynamically select the most relevant tools, knowledge documents,
and libraries for a given user query. Follows the BioMNI retriever pattern.
"""

import contextlib
import re

from cellagent.config import CellAgentConfig, default_config
from cellagent.llm import llm_chat


class ResourceRetriever:
    """Retrieve relevant resources for a user query using LLM-based selection."""

    def __init__(self, config: CellAgentConfig | None = None):
        self.config = config or default_config

    def retrieve(self, query: str, resources: dict) -> dict:
        """Select relevant resources for the query.

        Args:
            query: User's analysis query.
            resources: Dictionary with keys 'tools', 'knowledge', 'libraries',
                      each containing a list of available resources.

        Returns:
            Dictionary with same keys but filtered to relevant resources.
        """
        prompt = self._build_prompt(query, resources)

        try:
            response = llm_chat(
                messages=[{"role": "user", "content": prompt}],
                model=self.config.llm_model,
                temperature=0.0,
                max_tokens=1024,
                base_url=self.config.base_url,
                api_key=self.config.api_key,
            )
        except Exception:
            # Resource retrieval is an optimization. If the LLM is unavailable
            # or returns an unexpected response, keep the agent usable by
            # falling back to the complete resource set.
            return resources

        selected_indices = self._parse_response(response)
        selected = self._select_resources(resources, selected_indices)

        if not any(selected.get(key) for key in ("tools", "knowledge", "libraries")):
            return resources

        return selected

    def _build_prompt(self, query: str, resources: dict) -> str:
        """Build the retrieval prompt."""
        sections = [
            f"""You are an expert single-cell bioinformatics assistant. Select the most relevant resources to help answer a user's query.

USER QUERY: {query}

AVAILABLE TOOLS:
{self._format_items(resources.get("tools", []))}

AVAILABLE KNOWLEDGE DOCUMENTS:
{self._format_items(resources.get("knowledge", []))}

AVAILABLE LIBRARIES:
{self._format_items(resources.get("libraries", []))}

For each category, respond with ONLY the indices of relevant items:
TOOLS: [list of indices]
KNOWLEDGE: [list of indices]
LIBRARIES: [list of indices]

Example:
TOOLS: [0, 2, 5]
KNOWLEDGE: [0, 1]
LIBRARIES: [1, 3]

GUIDELINES:
1. Be generous - include resources that might be useful
2. Always include preprocessing tools for raw data queries
3. Include relevant knowledge documents for best practices
4. Include visualization tools when analysis results need to be shown
5. For annotation queries, include both marker-based and LLM-based tools
6. When in doubt, include the resource
"""
        ]
        return "\n".join(sections)

    def _format_items(self, items: list) -> str:
        """Format items for the prompt."""
        if not items:
            return "None available"
        lines = []
        for i, item in enumerate(items):
            if isinstance(item, dict):
                name = item.get("name", f"Item {i}")
                desc = item.get("description", "")
                lines.append(f"{i}. {name}: {desc}")
            else:
                lines.append(f"{i}. {item}")
        return "\n".join(lines)

    def _parse_response(self, response: str) -> dict:
        """Parse LLM response to extract selected indices."""
        indices = {"tools": [], "knowledge": [], "libraries": []}

        patterns = {
            "tools": r"TOOLS:\s*\[(.*?)\]",
            "knowledge": r"KNOWLEDGE:\s*\[(.*?)\]",
            "libraries": r"LIBRARIES:\s*\[(.*?)\]",
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, response, re.IGNORECASE)
            if match and match.group(1).strip():
                with contextlib.suppress(ValueError):
                    indices[key] = [
                        int(idx.strip())
                        for idx in match.group(1).split(",")
                        if idx.strip()
                    ]

        return indices

    def _select_resources(self, resources: dict, indices: dict) -> dict:
        """Select resources based on parsed indices."""
        selected = {}
        for key in ["tools", "knowledge", "libraries"]:
            items = resources.get(key, [])
            selected_idx = indices.get(key, [])
            selected[key] = [
                items[i] for i in selected_idx
                if 0 <= i < len(items)
            ]
        return selected
