"""Run provenance recording for CellAgent analyses.

The manifest is intentionally plain JSON so analysis runs can be audited,
shared, diffed, and regenerated without depending on CellAgent internals.
"""

from __future__ import annotations

import json
import os
import platform
import uuid
from datetime import datetime, timezone
from importlib import metadata
from typing import Any


class ProvenanceRecorder:
    """Collect and persist metadata for one CellAgent analysis run."""

    schema_version = "1.0"

    def __init__(self, base_output_dir: str, run_id: str | None = None):
        self.base_output_dir = os.path.abspath(base_output_dir)
        self.run_id = run_id or self._make_run_id()
        self.run_dir = os.path.join(self.base_output_dir, self.run_id)
        os.makedirs(self.run_dir, exist_ok=True)

        self.manifest: dict[str, Any] = {
            "schema_version": self.schema_version,
            "run_id": self.run_id,
            "created_at": self._now(),
            "completed_at": None,
            "status": "running",
            "base_output_dir": self.base_output_dir,
            "run_dir": self.run_dir,
            "query": None,
            "data_path": None,
            "config": {},
            "environment": self._environment_snapshot(),
            "selected_resources": {},
            "iterations": [],
            "final_answer": None,
            "final_adata_summary": None,
            "output_files": [],
        }

    def start_run(
        self,
        query: str,
        data_path: str | None,
        config: dict,
        selected_resources: dict,
    ) -> None:
        """Record immutable run inputs and selected domain resources."""
        self.manifest["query"] = query
        self.manifest["data_path"] = os.path.abspath(data_path) if data_path else None
        self.manifest["config"] = self._redact_config(config)
        self.manifest["selected_resources"] = self._summarize_resources(selected_resources)
        self.save()

    def record_iteration(
        self,
        iteration: int,
        response: str,
        code: str | None = None,
        execution_result: dict | None = None,
    ) -> None:
        """Append one ReAct iteration to the manifest."""
        entry: dict[str, Any] = {
            "iteration": iteration,
            "timestamp": self._now(),
            "assistant_response": response,
            "code": code,
        }
        if execution_result is not None:
            entry["execution"] = {
                "success": execution_result.get("success"),
                "duration_seconds": execution_result.get("duration_seconds"),
                "output": execution_result.get("output", ""),
                "error": execution_result.get("error", ""),
                "variables": execution_result.get("variables", []),
            }
        self.manifest["iterations"].append(entry)
        self.save()

    def finish_run(
        self,
        status: str,
        final_answer: str,
        final_adata_summary: str | None = None,
    ) -> str:
        """Finalize and write the manifest.

        Returns:
            Absolute path to the JSON manifest.
        """
        self.manifest["status"] = status
        self.manifest["completed_at"] = self._now()
        self.manifest["final_answer"] = final_answer
        self.manifest["final_adata_summary"] = final_adata_summary
        self.manifest["output_files"] = self._list_output_files()
        return self.save()

    def save(self) -> str:
        """Persist the manifest and return its path."""
        path = self.manifest_path
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.manifest, f, indent=2, ensure_ascii=False)
        return path

    @property
    def manifest_path(self) -> str:
        return os.path.join(self.run_dir, "run_manifest.json")

    def _list_output_files(self) -> list[dict[str, Any]]:
        files = []
        for root, _, names in os.walk(self.run_dir):
            for name in names:
                path = os.path.join(root, name)
                rel_path = os.path.relpath(path, self.run_dir)
                if rel_path == "run_manifest.json":
                    continue
                files.append({
                    "path": rel_path.replace(os.sep, "/"),
                    "size_bytes": os.path.getsize(path),
                })
        return sorted(files, key=lambda item: item["path"])

    @staticmethod
    def _make_run_id() -> str:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
        return f"run-{timestamp}-{uuid.uuid4().hex[:8]}"

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat()

    @staticmethod
    def _redact_config(config: dict) -> dict:
        redacted = dict(config)
        if redacted.get("api_key"):
            redacted["api_key"] = "***"
        return redacted

    @staticmethod
    def _summarize_resources(resources: dict) -> dict:
        return {
            "tools": [
                {"name": item.get("name"), "category": item.get("category")}
                for item in resources.get("tools", [])
            ],
            "knowledge": [
                {"id": item.get("id"), "name": item.get("name")}
                for item in resources.get("knowledge", [])
            ],
            "libraries": [
                {"name": item.get("name")}
                for item in resources.get("libraries", [])
            ],
        }

    @staticmethod
    def _environment_snapshot() -> dict:
        packages = {}
        for name in [
            "cellagent",
            "scanpy",
            "anndata",
            "numpy",
            "pandas",
            "scipy",
            "scikit-learn",
            "matplotlib",
            "seaborn",
            "openai",
        ]:
            try:
                packages[name] = metadata.version(name)
            except metadata.PackageNotFoundError:
                packages[name] = None

        return {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "packages": packages,
        }
