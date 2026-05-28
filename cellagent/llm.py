"""LLM factory for CellAgent.

Provides a unified interface to create LLM instances from various providers,
following the BioMNI pattern but simplified for the OpenAI-compatible API
available in this environment.
"""

from __future__ import annotations

import os
from typing import Literal


SourceType = Literal["OpenAI", "Custom"]


def get_llm_client(
    model: str | None = None,
    source: SourceType | None = None,
    base_url: str | None = None,
    api_key: str | None = None,
) -> tuple[OpenAI, str]:
    """Create an OpenAI client and return (client, model_name).

    This uses the pre-configured OpenAI-compatible API available in the
    sandbox environment. The OPENAI_API_KEY and base_url are pre-configured.

    Args:
        model: Model name to use. Defaults to gpt-4.1-mini.
        source: Provider source type.
        base_url: Custom base URL for the API.
        api_key: API key override.

    Returns:
        Tuple of (OpenAI client, model name string).
    """
    try:
        from openai import OpenAI
    except ImportError as exc:
        raise ImportError(
            "The 'openai' package is required for LLM-backed CellAgent features. "
            "Install the project with `pip install -e .` or install `openai`."
        ) from exc

    from cellagent.config import default_config

    if model is None:
        model = default_config.llm_model
    if api_key is None:
        api_key = default_config.api_key or os.environ.get("OPENAI_API_KEY")

    kwargs = {}
    if api_key:
        kwargs["api_key"] = api_key
    if base_url:
        kwargs["base_url"] = base_url

    client = OpenAI(**kwargs)
    return client, model


def llm_chat(
    messages: list[dict],
    model: str | None = None,
    temperature: float | None = None,
    max_tokens: int = 4096,
    client=None,
    stop: list[str] | None = None,
    base_url: str | None = None,
    api_key: str | None = None,
) -> str:
    """Send a chat completion request and return the response text.

    Args:
        messages: List of message dicts with 'role' and 'content'.
        model: Model name.
        temperature: Sampling temperature.
        max_tokens: Maximum tokens in response.
        client: Pre-created OpenAI client.
        stop: Stop sequences.

    Returns:
        The assistant's response text.
    """
    from cellagent.config import default_config

    if model is None:
        model = default_config.llm_model
    if temperature is None:
        temperature = default_config.temperature
    if client is None:
        client, model = get_llm_client(model=model, base_url=base_url, api_key=api_key)

    kwargs = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
        "max_tokens": max_tokens,
    }
    if stop:
        kwargs["stop"] = stop

    response = client.chat.completions.create(**kwargs)
    return response.choices[0].message.content
