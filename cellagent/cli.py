"""Command-line interface for CellAgent.

Provides both interactive chat mode and single-query execution mode.
"""

import argparse
import os
import sys


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="CellAgent: AI Agent for Single-Cell Multi-Omics Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  cellagent --interactive

  # Single query with data file
  cellagent --query "Analyze this PBMC dataset" --data pbmc3k.h5ad

  # Specify output directory
  cellagent --query "Run QC on the data" --data data.h5ad --output ./results

  # Use specific model
  cellagent --query "Cluster the cells" --model gpt-4.1-mini
""",
    )

    parser.add_argument(
        "--query", "-q", type=str, help="Analysis query to execute"
    )
    parser.add_argument(
        "--data", "-d", type=str, help="Path to input data file (.h5ad)"
    )
    parser.add_argument(
        "--output", "-o", type=str, default="./output",
        help="Output directory (default: ./output)"
    )
    parser.add_argument(
        "--model", "-m", type=str, default=None,
        help="LLM model to use (default: gpt-4.1-mini)"
    )
    parser.add_argument(
        "--interactive", "-i", action="store_true",
        help="Run in interactive chat mode"
    )
    parser.add_argument(
        "--no-retriever", action="store_true",
        help="Disable resource retriever (use all tools)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", default=True,
        help="Enable verbose output"
    )
    parser.add_argument(
        "--max-iterations", type=int, default=50,
        help="Maximum ReAct iterations (default: 50)"
    )

    args = parser.parse_args()

    # Import here to avoid slow startup for --help
    from cellagent.config import CellAgentConfig
    from cellagent.agent.cell_agent import CellAgent

    # Build config
    config = CellAgentConfig(
        output_path=args.output,
        verbose=args.verbose,
        max_iterations=args.max_iterations,
        use_resource_retriever=not args.no_retriever,
    )
    if args.model:
        config.llm_model = args.model

    # Initialize agent
    agent = CellAgent(config=config, output_dir=args.output)

    if args.interactive:
        _run_interactive(agent)
    elif args.query:
        result = agent.run(args.query, data_path=args.data)
        print(result)
    else:
        parser.print_help()
        sys.exit(1)


def _run_interactive(agent):
    """Run interactive chat mode."""
    print("\n" + "=" * 60)
    print("  CellAgent Interactive Mode")
    print("  Type 'quit' or 'exit' to end the session")
    print("  Type 'reset' to clear the conversation")
    print("  Type 'status' to see current state")
    print("=" * 60 + "\n")

    first_query = True
    while True:
        try:
            user_input = input("You: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye!")
            break

        if not user_input:
            continue

        if user_input.lower() in ("quit", "exit"):
            print("Goodbye!")
            break

        if user_input.lower() == "reset":
            agent.reset()
            print("Agent state reset.\n")
            first_query = True
            continue

        if user_input.lower() == "status":
            variables = agent.executor.list_variables()
            print(f"Variables: {variables}")
            adata_summary = agent.executor.get_adata_summary()
            if adata_summary:
                print(adata_summary)
            print(f"Iterations: {agent.iteration_count}")
            print()
            continue

        if first_query:
            result = agent.run(user_input)
            first_query = False
        else:
            result = agent.chat(user_input)

        print(f"\nCellAgent: {result}\n")


if __name__ == "__main__":
    main()
