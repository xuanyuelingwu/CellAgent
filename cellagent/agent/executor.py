"""Code execution engine for CellAgent.

Provides a persistent Python REPL environment for executing analysis code,
following the BioMNI pattern of maintaining state across execution steps.
"""

import io
import sys
import time
import traceback
from typing import Any


class CodeExecutor:
    """Persistent Python code execution environment.

    Maintains a shared namespace across executions, allowing variables
    (like AnnData objects) to persist between steps.
    """

    def __init__(self, output_dir: str = "./output", timeout_seconds: int = 600):
        """Initialize the executor with a fresh namespace.

        Args:
            output_dir: Default output directory for generated files.
        """
        self.namespace: dict[str, Any] = {}
        self.output_dir = output_dir
        self.timeout_seconds = timeout_seconds
        self.execution_history: list[dict] = []

        # Pre-import common libraries in the namespace
        self._setup_namespace()

    def _setup_namespace(self):
        """Set up the execution namespace with common imports."""
        output_dir_literal = repr(self.output_dir)
        base_setup_code = f"""
import os
import sys

# Set default output directory
OUTPUT_DIR = {output_dir_literal}
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Suppress warnings for cleaner output
import warnings
warnings.filterwarnings('ignore')
"""
        try:
            exec(base_setup_code, self.namespace)
        except Exception:
            pass  # Non-critical; execution will surface path issues later.

        optional_imports = [
            "import numpy as np",
            "import pandas as pd",
            "import matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt",
        ]
        for import_code in optional_imports:
            try:
                exec(import_code, self.namespace)
            except Exception:
                continue

    def execute(self, code: str) -> dict:
        """Execute Python code in the persistent namespace.

        Args:
            code: Python code string to execute.

        Returns:
            Dictionary with keys:
                - 'success': bool indicating if execution succeeded
                - 'output': captured stdout output
                - 'error': error message if failed
                - 'variables': list of new/modified variable names
        """
        # Capture stdout
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        captured_out = io.StringIO()
        captured_err = io.StringIO()

        result = {
            "success": False,
            "output": "",
            "error": "",
            "variables": [],
            "duration_seconds": None,
        }

        # Track existing variables
        existing_vars = set(self.namespace.keys())
        started = time.monotonic()

        try:
            self._validate_code(code)
            sys.stdout = captured_out
            sys.stderr = captured_err

            # Execute the code
            sys.settrace(self._build_timeout_tracer())
            exec(code, self.namespace)

            result["success"] = True
            result["output"] = captured_out.getvalue()

            # Track new/modified variables
            new_vars = set(self.namespace.keys()) - existing_vars
            result["variables"] = list(new_vars)

        except Exception as e:
            result["error"] = f"{type(e).__name__}: {str(e)}\n{traceback.format_exc()}"
            result["output"] = captured_out.getvalue()

        finally:
            sys.settrace(None)
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            result["duration_seconds"] = round(time.monotonic() - started, 6)

        # Record execution
        self.execution_history.append({
            "code": code,
            "success": result["success"],
            "duration_seconds": result["duration_seconds"],
            "output_preview": result["output"][:500] if result["output"] else "",
            "error_preview": result["error"][:500] if result["error"] else "",
        })

        return result

    def _build_timeout_tracer(self):
        """Create a lightweight line tracer that interrupts long Python code."""
        start = time.monotonic()

        def trace_calls(frame, event, arg):
            if event == "line" and time.monotonic() - start > self.timeout_seconds:
                raise TimeoutError(
                    f"Code execution exceeded {self.timeout_seconds} seconds"
                )
            return trace_calls

        return trace_calls

    def _validate_code(self, code: str) -> None:
        """Reject high-risk operations before executing LLM-generated code."""
        import ast

        blocked_calls = {
            ("os", "system"),
            ("os", "popen"),
            ("subprocess", "run"),
            ("subprocess", "Popen"),
            ("subprocess", "call"),
            ("subprocess", "check_call"),
            ("subprocess", "check_output"),
            ("shutil", "rmtree"),
        }
        blocked_names = {"eval", "exec", "compile", "__import__"}

        tree = ast.parse(code)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name) and node.func.id in blocked_names:
                    raise PermissionError(f"Blocked unsafe call: {node.func.id}")
                if isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name):
                    call = (node.func.value.id, node.func.attr)
                    if call in blocked_calls:
                        raise PermissionError(f"Blocked unsafe call: {'.'.join(call)}")

    def get_variable(self, name: str) -> Any:
        """Get a variable from the namespace."""
        return self.namespace.get(name)

    def set_variable(self, name: str, value: Any):
        """Set a variable in the namespace."""
        self.namespace[name] = value

    def list_variables(self) -> dict[str, str]:
        """List all user-defined variables with their types."""
        skip = {"__builtins__", "os", "sys", "np", "numpy", "pd", "pandas",
                "plt", "matplotlib", "warnings", "OUTPUT_DIR"}
        return {
            name: type(val).__name__
            for name, val in self.namespace.items()
            if name not in skip and not name.startswith("_")
        }

    def get_adata_summary(self, var_name: str = "adata") -> str | None:
        """Get a summary of an AnnData object in the namespace."""
        adata = self.namespace.get(var_name)
        if adata is None:
            return None

        try:
            summary = f"AnnData: {var_name}\n"
            summary += f"  Shape: {adata.n_obs} cells x {adata.n_vars} genes\n"
            summary += f"  Obs columns: {list(adata.obs.columns)}\n"
            summary += f"  Var columns: {list(adata.var.columns)}\n"
            summary += f"  Obsm keys: {list(adata.obsm.keys())}\n"
            summary += f"  Uns keys: {list(adata.uns.keys())[:20]}\n"
            if adata.raw is not None:
                summary += f"  Raw: {adata.raw.shape}\n"
            return summary
        except Exception:
            return f"Variable '{var_name}' exists but is not a valid AnnData object."

    def reset(self):
        """Reset the execution environment."""
        self.namespace.clear()
        self.execution_history.clear()
        self._setup_namespace()
