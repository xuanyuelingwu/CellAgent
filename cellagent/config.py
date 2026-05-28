"""CellAgent configuration management.

Centralized configuration for the single-cell multi-omics AI agent,
following the BioMNI pattern of dataclass-based config with environment
variable overrides.
"""

import os
from dataclasses import dataclass, field


@dataclass
class CellAgentConfig:
    """Global configuration for CellAgent."""

    # Data paths
    data_path: str = "./data"
    output_path: str = "./output"

    # LLM settings
    llm_model: str = "gpt-4.1-mini"
    llm_source: str = "OpenAI"
    temperature: float = 0.1
    base_url: str | None = None
    api_key: str | None = None

    # Agent settings
    timeout_seconds: int = 600
    max_iterations: int = 50
    use_resource_retriever: bool = True

    # Execution settings
    verbose: bool = True

    # Single-cell specific defaults
    default_min_genes: int = 200
    default_min_cells: int = 3
    default_n_top_genes: int = 2000
    default_resolution: float = 0.5
    default_n_neighbors: int = 15
    default_n_pcs: int = 50

    def __post_init__(self):
        """Override defaults with environment variables if set."""
        env_mappings = {
            "CELLAGENT_DATA_PATH": "data_path",
            "CELLAGENT_OUTPUT_PATH": "output_path",
            "CELLAGENT_LLM_MODEL": "llm_model",
            "CELLAGENT_LLM_SOURCE": "llm_source",
            "CELLAGENT_BASE_URL": "base_url",
            "CELLAGENT_TIMEOUT": "timeout_seconds",
            "OPENAI_API_KEY": "api_key",
        }

        for env_var, attr in env_mappings.items():
            val = os.environ.get(env_var)
            if val is not None:
                current = getattr(self, attr)
                if isinstance(current, int):
                    setattr(self, attr, int(val))
                elif isinstance(current, float):
                    setattr(self, attr, float(val))
                elif isinstance(current, bool):
                    setattr(self, attr, val.lower() in ("true", "1", "yes"))
                else:
                    setattr(self, attr, val)

    def to_dict(self) -> dict:
        """Serialize key runtime settings."""
        return {
            "data_path": self.data_path,
            "output_path": self.output_path,
            "llm_model": self.llm_model,
            "llm_source": self.llm_source,
            "temperature": self.temperature,
            "base_url": self.base_url,
            "timeout_seconds": self.timeout_seconds,
            "max_iterations": self.max_iterations,
            "use_resource_retriever": self.use_resource_retriever,
        }


# Global default config instance
default_config = CellAgentConfig()
