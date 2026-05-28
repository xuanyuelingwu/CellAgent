# Single-Cell Study Design Guardrails

## Metadata

**short_description**: Practical guardrails for reproducible and statistically cautious single-cell biomedical analysis.

---

## Overview

Single-cell workflows should separate exploratory cell-level analysis from
sample-level biological claims. Clusters, marker genes, and UMAP structure are
useful for hypothesis generation, but condition-level conclusions require
appropriate sample or donor replication, metadata checks, and statistical models
that respect the experimental design.

## Key Checks Before Analysis

- Confirm whether the dataset has sample, donor, batch, tissue, condition, and
  technology metadata.
- Check whether condition and batch are confounded before interpreting disease,
  treatment, or tissue differences.
- Record QC thresholds, normalization choices, highly variable gene settings,
  dimensionality reduction parameters, clustering resolution, random seeds, and
  package versions.
- Save both intermediate AnnData objects and final plots when a workflow will be
  reviewed or rerun.

## Differential Expression Guidance

Cluster marker discovery with cell-level methods such as rank_genes_groups is
best treated as exploratory evidence for annotation. It should not be used as
the sole basis for condition-level claims when cells come from replicated
samples or donors.

For condition differential expression, prefer a sample-aware workflow:

- Aggregate counts by sample or donor within each cell type or cluster.
- Model condition effects at the sample level, while considering batch and other
  covariates.
- Report the number of biological replicates per group.
- Clearly label results as exploratory when replication is missing or metadata
  are incomplete.

## Annotation Guidance

Cell type labels should be supported by multiple marker genes, reference
mapping, or external biological evidence. When marker evidence is mixed, report
candidate labels and uncertainty instead of forcing a single confident label.

## Agent Behavior

When metadata are missing, the agent should ask for them or state the limitation.
When a requested analysis would overstate the evidence, the agent should offer a
safer exploratory alternative and explain what extra information would be needed
for a stronger conclusion.
