"""CellAgent Quick Start Example.

Demonstrates the basic usage of CellAgent for analyzing a PBMC dataset.
This example uses the Scanpy built-in PBMC3k dataset.
"""

import os
import sys

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cellagent import CellAgent, CellAgentConfig


def main():
    """Run a quick start analysis on PBMC3k data."""

    # Configure the agent
    config = CellAgentConfig(
        llm_model="gpt-4.1-mini",
        temperature=0.1,
        max_iterations=20,
        verbose=True,
    )

    # Initialize the agent
    agent = CellAgent(config=config, output_dir="./output/quickstart")

    # Run analysis with a natural language query
    result = agent.run(
        query="""
        Please perform a complete analysis of the PBMC3k dataset:
        1. Load the PBMC3k dataset from scanpy.datasets
        2. Run quality control and filter low-quality cells
        3. Normalize and select highly variable genes
        4. Run PCA, compute neighbors, and perform Leiden clustering
        5. Compute UMAP and visualize clusters
        6. Find marker genes for each cluster
        7. Annotate cell types using known PBMC markers
        8. Generate a summary report
        """,
    )

    print(result)


def example_with_custom_data():
    """Example: Analyze a custom h5ad file."""

    config = CellAgentConfig(verbose=True, max_iterations=15)
    agent = CellAgent(config=config, output_dir="./output/custom")

    # Point to your data file
    result = agent.run(
        query="Load the data, run QC, normalize, cluster, and annotate cell types. "
              "This is a human liver biopsy sample.",
        data_path="/path/to/your/data.h5ad",
    )
    print(result)


def example_follow_up():
    """Example: Interactive follow-up queries."""

    config = CellAgentConfig(verbose=True)
    agent = CellAgent(config=config, output_dir="./output/interactive")

    # Initial analysis
    agent.run(
        query="Load the PBMC3k dataset and perform basic preprocessing "
              "(QC, normalize, HVG, PCA, neighbors, clustering, UMAP)."
    )

    # Follow-up questions
    response = agent.chat("What are the top marker genes for cluster 0?")
    print(response)

    response = agent.chat("Can you try a higher clustering resolution of 1.0?")
    print(response)

    response = agent.chat("Generate a dot plot of T cell and B cell markers.")
    print(response)


if __name__ == "__main__":
    main()
