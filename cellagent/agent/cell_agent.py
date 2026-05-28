"""CellAgent: AI Agent for Single-Cell Multi-Omics Analysis.

The core agent class implementing a ReAct (Reasoning + Acting) loop
for autonomous single-cell data analysis. Inspired by Stanford BioMNI's
architecture with domain-specific adaptations for single-cell genomics.
"""

import json
import os
import re
import time
from typing import Any

from cellagent.config import CellAgentConfig, default_config
from cellagent.llm import llm_chat
from cellagent.tools import build_default_registry, ALL_TOOL_FUNCTIONS
from cellagent.tools.tool_registry import ToolRegistry
from cellagent.knowledge.loader import KnowledgeLoader
from cellagent.retriever.resource_retriever import ResourceRetriever
from cellagent.agent.executor import CodeExecutor
from cellagent.agent.provenance import ProvenanceRecorder


# Available Python libraries for single-cell analysis
AVAILABLE_LIBRARIES = [
    {"name": "scanpy", "description": "Core single-cell analysis library (preprocessing, clustering, DE, visualization)"},
    {"name": "anndata", "description": "Data structure for annotated single-cell data (AnnData objects)"},
    {"name": "numpy", "description": "Numerical computing (arrays, linear algebra, statistics)"},
    {"name": "pandas", "description": "Data manipulation and analysis (DataFrames, CSV I/O)"},
    {"name": "matplotlib", "description": "Plotting and visualization library"},
    {"name": "seaborn", "description": "Statistical data visualization built on matplotlib"},
    {"name": "scipy", "description": "Scientific computing (statistics, sparse matrices, optimization)"},
    {"name": "sklearn", "description": "Machine learning (clustering metrics, PCA, classification)"},
    {"name": "harmonypy", "description": "Harmony batch correction for single-cell data"},
    {"name": "gseapy", "description": "Gene set enrichment analysis (Enrichr, GSEA)"},
    {"name": "scrublet", "description": "Doublet detection for single-cell RNA-seq"},
]


class CellAgent:
    """AI Agent for single-cell multi-omics analysis.

    Implements a ReAct loop that:
    1. Understands the user's analysis request
    2. Retrieves relevant tools and domain knowledge
    3. Generates and executes Python code iteratively
    4. Returns results with explanations

    Attributes:
        config: Agent configuration.
        registry: Tool registry with all available analysis tools.
        knowledge: Knowledge base loader.
        retriever: Resource retriever for dynamic tool selection.
        executor: Persistent code execution environment.
    """

    SYSTEM_PROMPT = """You are CellAgent, an expert AI assistant specialized in single-cell multi-omics data analysis.

You have deep expertise in:
- Single-cell RNA-seq (scRNA-seq) data analysis using Scanpy
- Quality control, normalization, and preprocessing
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Cell clustering (Leiden, Louvain)
- Cell type annotation using marker genes and LLM-assisted methods
- Differential expression analysis
- Trajectory and pseudotime analysis
- Multi-sample integration and batch correction
- Gene set enrichment analysis
- Publication-quality visualization

IMPORTANT RULES:
1. You operate in a ReAct (Reasoning + Acting) loop
2. For each step, first THINK about what to do, then write executable Python CODE
3. Your code runs in a persistent Python environment - variables persist between steps
4. Always use the provided tools and follow best practices from the knowledge base
5. Generate clear, well-commented code
6. After each code execution, analyze the output and decide the next step
7. Save all plots to the output directory
8. Provide clear explanations of your analysis choices

SINGLE-CELL BIOMEDICINE GUARDRAILS:
1. Distinguish exploratory cell-level marker discovery from sample/donor-level
   condition differential expression. Do not present cluster markers as
   condition-level biological conclusions.
2. Before condition comparisons, check whether sample, donor, batch, tissue,
   and condition metadata are available and whether replication is adequate.
3. Prefer explicit parameters, random seeds, and saved intermediate artifacts
   so analyses can be inspected and rerun.
4. When cell type annotation is uncertain, report marker evidence and caveats
   instead of overstating a label.
5. Treat LLM-generated analysis as a draft that needs domain review for
   experimental design, QC thresholds, and biological interpretation.

RESPONSE FORMAT:
For each step, respond with:

THOUGHT: [Your reasoning about what to do next]

CODE:
```python
[Your Python code here]
```

When the analysis is complete, respond with:

THOUGHT: [Summary of what was accomplished]

ANSWER: [Final summary of results and findings]

AVAILABLE TOOLS AND FUNCTIONS:
{tool_descriptions}

DOMAIN KNOWLEDGE:
{knowledge_context}
"""

    def __init__(
        self,
        config: CellAgentConfig | None = None,
        knowledge_dir: str | None = None,
        output_dir: str = "./output",
        run_id: str | None = None,
        create_run_dir: bool = True,
    ):
        """Initialize CellAgent.

        Args:
            config: Agent configuration. Uses default if None.
            knowledge_dir: Custom knowledge directory path.
            output_dir: Base directory for output files and plots.
            run_id: Optional stable run identifier for provenance.
            create_run_dir: If True, create an isolated run subdirectory.
        """
        self.config = config or default_config
        self.base_output_dir = output_dir
        self.provenance = ProvenanceRecorder(output_dir, run_id=run_id)
        self.output_dir = self.provenance.run_dir if create_run_dir else os.path.abspath(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

        # Initialize components
        self.registry = build_default_registry()
        self.knowledge = KnowledgeLoader(knowledge_dir)
        self.retriever = ResourceRetriever(config=self.config)
        self.executor = CodeExecutor(
            output_dir=self.output_dir,
            timeout_seconds=self.config.timeout_seconds,
        )

        # Inject tool functions into executor namespace
        for name, fn in ALL_TOOL_FUNCTIONS.items():
            self.executor.set_variable(name, fn)

        # Conversation history
        self.messages: list[dict] = []
        self.iteration_count = 0

        if self.config.verbose:
            print(f"CellAgent initialized with {len(self.registry)} tools "
                  f"and {len(self.knowledge.documents)} knowledge documents")
            print(f"Run ID: {self.provenance.run_id}")

    def configure(self, query: str) -> dict:
        """Configure the agent for a specific query by retrieving relevant resources.

        Args:
            query: User's analysis query.

        Returns:
            Dictionary of selected resources.
        """
        # Build resource lists for retrieval
        resources = {
            "tools": self.registry.get_summary_for_retrieval(),
            "knowledge": self.knowledge.get_document_summaries(),
            "libraries": AVAILABLE_LIBRARIES,
        }

        if self.config.use_resource_retriever:
            selected = self.retriever.retrieve(query, resources)
        else:
            # Use all resources
            selected = resources

        return selected

    def _build_system_prompt(self, selected_resources: dict) -> str:
        """Build the system prompt with selected tools and knowledge."""
        # Format tool descriptions
        tool_names = [t["name"] for t in selected_resources.get("tools", [])]
        tool_desc = self.registry.get_prompt_description(tool_names if tool_names else None)

        # Format knowledge context
        knowledge_parts = []
        for doc_summary in selected_resources.get("knowledge", []):
            doc = self.knowledge.get_document_by_id(doc_summary["id"])
            if doc:
                knowledge_parts.append(
                    f"### {doc['name']}\n{doc['content_without_metadata'][:2000]}"
                )

        knowledge_context = "\n\n".join(knowledge_parts) if knowledge_parts else "No specific knowledge documents selected."

        return self.SYSTEM_PROMPT.format(
            tool_descriptions=tool_desc,
            knowledge_context=knowledge_context,
        )

    def _extract_code(self, response: str) -> str | None:
        """Extract Python code from the LLM response."""
        # Try to find code block
        patterns = [
            r"```python\n(.*?)```",
            r"```\n(.*?)```",
            r"CODE:\n```python\n(.*?)```",
            r"CODE:\n```\n(.*?)```",
        ]

        for pattern in patterns:
            match = re.search(pattern, response, re.DOTALL)
            if match:
                return match.group(1).strip()

        return None

    def _extract_answer(self, response: str) -> str | None:
        """Extract final answer from the LLM response."""
        match = re.search(r"ANSWER:\s*(.*?)$", response, re.DOTALL)
        if match:
            return match.group(1).strip()
        return None

    def _format_execution_result(self, result: dict) -> str:
        """Format code execution result for the conversation."""
        if result["success"]:
            output = result["output"]
            if len(output) > 3000:
                output = output[:3000] + "\n... [output truncated]"
            msg = f"EXECUTION RESULT (success):\n{output}"
            if result["variables"]:
                msg += f"\nNew variables: {result['variables']}"
        else:
            msg = f"EXECUTION ERROR:\n{result['error']}"
            if result["output"]:
                msg += f"\nPartial output:\n{result['output'][:1000]}"

        # Add current state info
        adata_summary = self.executor.get_adata_summary()
        if adata_summary:
            msg += f"\n\nCurrent AnnData state:\n{adata_summary}"

        return msg

    def run(self, query: str, data_path: str | None = None) -> str:
        """Run the agent on a user query.

        This is the main entry point. It sets up the agent, then runs
        the ReAct loop until completion or max iterations.

        Args:
            query: User's analysis request.
            data_path: Optional path to the data file.

        Returns:
            Final analysis result as a string.
        """
        print(f"\n{'='*60}")
        print(f"CellAgent: Processing query...")
        print(f"{'='*60}\n")

        # Configure resources for this query
        selected_resources = self.configure(query)

        if self.config.verbose:
            n_tools = len(selected_resources.get("tools", []))
            n_knowledge = len(selected_resources.get("knowledge", []))
            print(f"Selected {n_tools} tools and {n_knowledge} knowledge documents")
            print(f"Output directory: {self.output_dir}")

        self.provenance.start_run(
            query=query,
            data_path=data_path,
            config=self.config.to_dict(),
            selected_resources=selected_resources,
        )

        # Build system prompt
        system_prompt = self._build_system_prompt(selected_resources)

        # Initialize conversation
        self.messages = [{"role": "system", "content": system_prompt}]

        # Build user message
        user_msg = query
        if data_path:
            user_msg += (
                f"\n\nData file path: {data_path}"
                "\nThe same path is available in the Python environment as DATA_PATH. "
                "Load it into an AnnData variable named `adata` before analysis."
            )
            self.executor.set_variable("DATA_PATH", data_path)

        # Add current environment state
        env_state = f"\nOutput directory: {self.output_dir}"
        variables = self.executor.list_variables()
        if variables:
            env_state += f"\nExisting variables: {variables}"
        user_msg += env_state

        self.messages.append({"role": "user", "content": user_msg})

        # ReAct loop
        final_answer = None
        self.iteration_count = 0

        while self.iteration_count < self.config.max_iterations:
            self.iteration_count += 1

            if self.config.verbose:
                print(f"\n--- Iteration {self.iteration_count} ---")

            # Get LLM response
            try:
                response = llm_chat(
                    messages=self.messages,
                    model=self.config.llm_model,
                    temperature=self.config.temperature,
                    max_tokens=4096,
                    base_url=self.config.base_url,
                    api_key=self.config.api_key,
                )
            except Exception as e:
                print(f"LLM error: {e}")
                break

            if self.config.verbose:
                # Print thought
                thought_match = re.search(r"THOUGHT:\s*(.*?)(?=CODE:|ANSWER:|$)", response, re.DOTALL)
                if thought_match:
                    print(f"THOUGHT: {thought_match.group(1).strip()[:200]}")

            # Check for final answer
            answer = self._extract_answer(response)
            if answer:
                final_answer = answer
                self.provenance.record_iteration(
                    iteration=self.iteration_count,
                    response=response,
                )
                if self.config.verbose:
                    print(f"\nFINAL ANSWER: {answer[:300]}...")
                break

            # Extract and execute code
            code = self._extract_code(response)
            if code:
                if self.config.verbose:
                    print(f"Executing code ({len(code)} chars)...")

                exec_result = self.executor.execute(code)
                result_msg = self._format_execution_result(exec_result)
                self.provenance.record_iteration(
                    iteration=self.iteration_count,
                    response=response,
                    code=code,
                    execution_result=exec_result,
                )

                if self.config.verbose:
                    status = "SUCCESS" if exec_result["success"] else "ERROR"
                    print(f"Execution: {status}")
                    if exec_result["output"]:
                        print(f"Output preview: {exec_result['output'][:200]}")

                # Add to conversation
                self.messages.append({"role": "assistant", "content": response})
                self.messages.append({"role": "user", "content": result_msg})
            else:
                # No code and no answer - might be just thinking
                self.provenance.record_iteration(
                    iteration=self.iteration_count,
                    response=response,
                )
                self.messages.append({"role": "assistant", "content": response})
                self.messages.append({
                    "role": "user",
                    "content": "Please continue with the analysis. Write Python code for the next step."
                })

        if final_answer is None:
            final_answer = "Analysis completed (max iterations reached). Check the output directory for results."

        manifest_path = self.provenance.finish_run(
            status="completed",
            final_answer=final_answer,
            final_adata_summary=self.executor.get_adata_summary(),
        )

        # Generate summary
        summary = self._generate_summary(query, final_answer, manifest_path)
        return summary

    def _generate_summary(
        self,
        query: str,
        answer: str,
        manifest_path: str | None = None,
    ) -> str:
        """Generate a formatted summary of the analysis."""
        lines = [
            "=" * 60,
            "CELLAGENT ANALYSIS REPORT",
            "=" * 60,
            f"\nQuery: {query}",
            f"Run ID: {self.provenance.run_id}",
            f"Iterations: {self.iteration_count}",
            f"Output directory: {self.output_dir}",
        ]
        if manifest_path:
            lines.append(f"Run manifest: {manifest_path}")
        lines.append("")

        # List generated files
        if os.path.exists(self.output_dir):
            files = os.listdir(self.output_dir)
            if files:
                lines.append("Generated files:")
                for f in sorted(files):
                    fpath = os.path.join(self.output_dir, f)
                    size = os.path.getsize(fpath)
                    lines.append(f"  - {f} ({size:,} bytes)")
                lines.append("")

        # Add AnnData summary
        adata_summary = self.executor.get_adata_summary()
        if adata_summary:
            lines.append("Final dataset state:")
            lines.append(adata_summary)
            lines.append("")

        lines.append("Analysis Result:")
        lines.append(answer)

        return "\n".join(lines)

    def chat(self, message: str) -> str:
        """Send a follow-up message in the current conversation.

        Args:
            message: Follow-up question or instruction.

        Returns:
            Agent's response.
        """
        if not self.messages:
            return self.run(message)

        self.messages.append({"role": "user", "content": message})

        response = llm_chat(
            messages=self.messages,
            model=self.config.llm_model,
            temperature=self.config.temperature,
            max_tokens=4096,
            base_url=self.config.base_url,
            api_key=self.config.api_key,
        )

        # Check for code to execute
        code = self._extract_code(response)
        if code:
            exec_result = self.executor.execute(code)
            result_msg = self._format_execution_result(exec_result)
            self.provenance.record_iteration(
                iteration=self.iteration_count + 1,
                response=response,
                code=code,
                execution_result=exec_result,
            )
            self.messages.append({"role": "assistant", "content": response})
            self.messages.append({"role": "user", "content": result_msg})

            # Get follow-up response
            follow_up = llm_chat(
                messages=self.messages,
                model=self.config.llm_model,
                temperature=self.config.temperature,
                max_tokens=4096,
                base_url=self.config.base_url,
                api_key=self.config.api_key,
            )
            self.messages.append({"role": "assistant", "content": follow_up})
            self.provenance.record_iteration(
                iteration=self.iteration_count + 2,
                response=follow_up,
            )
            return follow_up
        else:
            self.messages.append({"role": "assistant", "content": response})
            self.provenance.record_iteration(
                iteration=self.iteration_count + 1,
                response=response,
            )
            return response

    def reset(self):
        """Reset the agent state."""
        self.messages.clear()
        self.executor.reset()
        self.iteration_count = 0
        # Re-inject tool functions
        for name, fn in ALL_TOOL_FUNCTIONS.items():
            self.executor.set_variable(name, fn)
