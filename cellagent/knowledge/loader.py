"""Knowledge base loader for CellAgent.

Loads and manages domain-specific know-how documents following the
BioMNI pattern of markdown-based knowledge management.
"""

import glob
import os
from pathlib import Path


class KnowledgeLoader:
    """Load and manage know-how documents for single-cell analysis."""

    def __init__(self, knowledge_dir: str | None = None):
        """Initialize the knowledge loader.

        Args:
            knowledge_dir: Directory containing knowledge documents.
                          If None, uses the default knowledge directory.
        """
        if knowledge_dir is None:
            current_dir = Path(__file__).parent
            knowledge_dir = str(current_dir)

        self.knowledge_dir = knowledge_dir
        self.documents: dict[str, dict] = {}
        self._load_documents()

    def _load_documents(self):
        """Load all markdown documents from the knowledge directory."""
        pattern = os.path.join(self.knowledge_dir, "*.md")
        md_files = glob.glob(pattern)

        for filepath in md_files:
            filename = os.path.basename(filepath)
            filename_no_ext = os.path.splitext(filename)[0]

            # Skip meta files
            if filename.upper() in ["README.MD"] or filename_no_ext.isupper():
                continue

            with open(filepath) as f:
                content = f.read()

            title, description, metadata = self._extract_metadata(content, filename)

            if "short_description" in metadata and metadata["short_description"]:
                description = metadata["short_description"]

            doc_id = filename_no_ext
            self.documents[doc_id] = {
                "id": doc_id,
                "name": title,
                "description": description,
                "content": content,
                "content_without_metadata": self._strip_metadata(content),
                "filepath": filepath,
                "metadata": metadata,
            }

    def _extract_metadata(self, content: str, filename: str) -> tuple[str, str, dict]:
        """Extract title, description, and metadata from markdown content."""
        lines = content.split("\n")

        # Extract title
        title = None
        for line in lines:
            if line.startswith("# "):
                title = line[2:].strip()
                break
        if title is None:
            title = filename.replace("_", " ").replace(".md", "").title()

        # Extract metadata section
        metadata = {}
        in_metadata = False
        current_field = None

        for line in lines:
            if line.startswith("## Metadata"):
                in_metadata = True
                continue
            elif in_metadata:
                if line.startswith("##") and "Metadata" not in line:
                    break
                elif line.startswith("**") and "**:" in line:
                    field_match = line.split("**")[1]
                    current_field = field_match.lower().replace(" ", "_")
                    colon_idx = line.find("**:")
                    if colon_idx != -1:
                        value_part = line[colon_idx + 3:].strip()
                        metadata[current_field] = value_part
                elif current_field and line.strip() and not line.startswith("---"):
                    if current_field not in metadata:
                        metadata[current_field] = ""
                    if line.startswith("- "):
                        if metadata[current_field]:
                            metadata[current_field] += ", " + line[2:].strip()
                        else:
                            metadata[current_field] = line[2:].strip()

        # Extract description from Overview section
        description = ""
        in_overview = False
        overview_lines = []

        for line in lines:
            if line.startswith("## Overview"):
                in_overview = True
                continue
            elif in_overview:
                if line.startswith("##"):
                    break
                elif line.strip():
                    overview_lines.append(line.strip())

        if overview_lines:
            description = " ".join(overview_lines)
        else:
            found_title = False
            for line in lines:
                if line.startswith("# "):
                    found_title = True
                    continue
                if found_title and line.strip() and not line.startswith("#"):
                    description = line.strip()
                    break

        if len(description) > 200:
            description = description[:197] + "..."

        return title, description, metadata

    def _strip_metadata(self, content: str) -> str:
        """Strip the metadata section from document content."""
        lines = content.split("\n")
        result_lines = []
        in_metadata = False

        for line in lines:
            if line.startswith("## Metadata"):
                in_metadata = True
                continue
            if in_metadata:
                if line.strip() == "---":
                    in_metadata = False
                    continue
                if line.startswith("##") and "Metadata" not in line:
                    in_metadata = False
                    result_lines.append(line)
                continue
            result_lines.append(line)

        result = "\n".join(result_lines)
        while "\n\n\n\n" in result:
            result = result.replace("\n\n\n\n", "\n\n\n")
        return result.strip()

    def get_all_documents(self) -> list[dict]:
        """Get all knowledge documents as a list."""
        return list(self.documents.values())

    def get_document_by_id(self, doc_id: str) -> dict | None:
        """Get a specific document by ID."""
        return self.documents.get(doc_id)

    def get_document_summaries(self) -> list[dict]:
        """Get summaries of all documents (without full content)."""
        return [
            {"id": doc["id"], "name": doc["name"], "description": doc["description"]}
            for doc in self.documents.values()
        ]

    def add_custom_document(self, doc_id: str, name: str, description: str,
                            content: str, metadata: dict | None = None):
        """Add a custom knowledge document programmatically."""
        if metadata is None:
            metadata = {}
        self.documents[doc_id] = {
            "id": doc_id,
            "name": name,
            "description": description,
            "content": content,
            "content_without_metadata": content,
            "filepath": None,
            "metadata": metadata,
        }

    def reload(self):
        """Reload all documents from disk."""
        self.documents = {}
        self._load_documents()
