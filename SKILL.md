---
name: analyze-bulk-transcriptomics
description: Analyze bulk RNA-seq gene-count matrices with a reproducible PyDESeq2 workflow, gene-symbol-annotated outputs, QC, PCA, differential expression, volcano plots, and separate up/down pathway enrichment. Use for featureCounts or tabular gene counts; do not use for single-cell RNA-seq or FASTQ alignment/quantification.
---

# Analyze bulk transcriptomics

Use the bundled code instead of reimplementing statistical methods. Preserve defaults unless the user asks to change them.

## Workflow

1. Inspect count-matrix orientation, metadata columns, sample identifiers, and available gene annotations.
2. Read [references/input-and-options.md](references/input-and-options.md) for schemas and switches. Read [references/methods.md](references/methods.md) when reporting or changing methods.
3. Reuse `scripts/transcriptomics.py` and `scripts/run_transcriptomics.py`. The runner uses the bundled Ensembl gene-name and biotype mappings by default. If copying the scripts into a project, copy both mapping pickle files beside the `scripts/` directory or configure explicit mapping paths.
4. Create a JSON configuration from `references/example-config.json`. Prefer explicit paths and contrast levels.
5. Run `python run_transcriptomics.py --config analysis.json` from the project.
6. Inspect `run_summary.json`, sample alignment, retained-gene count, PCA, library sizes, detected genes, correlation heatmap, DE diagnostics, and the always-exported top-ten up/down enrichment tables and plots before interpreting biology.
7. Record every non-default choice in the response and generated `methods.txt`.

## Statistical guardrails

- Supply raw, non-negative integer counts to PyDESeq2 differential expression. Never run DESeq2 on VST, TPM, CPM, or log-transformed values.
- Treat TPM as visualization/expression normalization only. Require gene lengths for TPM; continue to use raw counts for differential expression.
- For paired analysis, include subject before condition, for example `~subject + condition`. Verify each retained subject has the requested levels; report incomplete pairs.
- Define contrast direction explicitly as `[factor, test, reference]`; positive log2 fold change means higher in `test`.
- Do not silently remove samples. Extra count columns may be excluded when absent from metadata, but report them; error when a metadata sample has no counts.
- Do not silently change notebook-derived defaults: count >=10 in at least 3 samples, protein-coding genes, VST with `use_design=False`, 500 most variable genes for PCA, per-gene standardization before PCA, Cook's refitting enabled, Cook's filtering disabled for the contrast, and independent filtering disabled.
- If the model matrix is not full rank, stop and explain which design variables are confounded. Do not drop a covariate without user approval.

## Default visual language

Use editable vector PDF plus PNG. Keep PDF fonts editable (`matplotlib.rc("pdf", fonttype=42)`). Use seaborn white-grid styling, 50-point PCA markers with alpha 0.8, variance percentages on PCA axes, condition-consistent colors across figures, and significant/non-significant volcano colors based on adjusted p-values even though the volcano y-axis shows nominal p-values. Accept named palettes or category-to-color mappings.

## Modular execution

Run the complete workflow by default. When the user asks for only part, import functions from `scripts/transcriptomics.py` rather than running unrelated stages. Common entry points are `read_counts`, `align_samples`, `filter_low_counts`, `filter_biotype`, `normalize_vst`, `normalize_tpm`, `compute_pca`, `plot_qc`, `fit_deseq2`, and `plot_volcano`.

Refresh the bundled human Ensembl mappings only when requested by running `scripts/update_ensembl_mappings.py`; read [references/input-and-options.md](references/input-and-options.md) first. Treat pickle files as trusted executable data and never substitute an untrusted pickle silently.

## Deliverables

Return the configuration, scripts or notebook that called them, tables, figures, `methods.txt`, and `run_summary.json`. State whether the run used notebook defaults or list overrides. Flag design limitations, outlying samples, incomplete pairs, low library sizes, and mismatches between metadata and counts without overstating exclusion decisions.
