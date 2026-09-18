# Analyze Bulk Transcriptomics

A Codex skill for reproducible bulk RNA-seq analysis from raw gene-count matrices. It runs a PyDESeq2 workflow covering metadata sample annotation, expression filtering, normalization, quality control, PCA, differential expression, and separate pathway enrichment for upregulated and downregulated genes.

## What it does

* Reads featureCounts output or tabular gene-count matrices.
* Validates raw counts and aligns samples explicitly to metadata.
* Reports count columns without metadata instead of silently including them.
* Filters genes at count ≥10 in at least three samples by default.
* Retains protein-coding genes using bundled human Ensembl annotations.
* Uses design-blind DESeq2 VST for QC and PCA, or TPM when requested.
* Supports unpaired, covariate-adjusted, and paired designs.
* Runs explicit PyDESeq2 contrasts on filtered raw counts.
* Exports gene-level tables with Ensembl IDs as the index and gene symbols in a separate column.
* Produces library-depth, detected-gene, correlation, PCA, volcano, and pathway figures as PNG and editable-text PDF.
* Runs MSigDB Hallmark over-representation analysis separately for significant positive- and negative-fold-change genes.
* Exports complete pathway results and the top ten pathways per direction, and reports FDR.
* Records configuration choices, package versions, metadata sample annotation, and run summaries.

## Inputs

The workflow requires:

1. A raw, non-negative integer gene-count matrix.
2. A CSV or TSV metadata table with one row per analyzed sample.
3. An explicit contrast such as `condition: treated vs control`.

Supported count layouts:

* `featurecounts`: standard featureCounts output with genes as rows.
* `genes\_rows`: delimited table with genes as rows and samples as columns.
* `samples\_rows`: delimited table with samples as rows and genes as columns.

This skill does not perform FASTQ alignment or quantification and is not intended for single-cell RNA-seq.

## Installation

Clone the repository into your agent’s skills directory:

```bash
git clone https://github.com/<OWNER>/<REPOSITORY>.git \
  <SKILLS_DIRECTORY>/analyze-bulk-transcriptomics
```

Restart the agent or begin a new session to load the skill.

## Default analysis choices

|Step|Default|
|-|-|
|Low-count filter|Count ≥10 in at least 3 samples|
|Biotype|Protein-coding|
|QC normalization|DESeq2 VST with `use\_design=False`|
|PCA|500 most variable genes, standardized per gene|
|DE model|PyDESeq2 0.5.4, parametric fit|
|Cook's refitting|Enabled|
|Contrast Cook's filtering|Disabled|
|Independent filtering|Disabled|
|Significance|Adjusted p < 0.05|
|Pathways|MSigDB Hallmark `h.all`, release `2024.1.Hs`|
|Pathway display|Top 10 per direction, with FDR threshold marked|

## Outputs

Each run creates:

```text
transcriptomics\_results/
├── figures/
│   ├── library\_sizes.{pdf,png}
│   ├── detected\_genes.{pdf,png}
│   ├── sample\_correlation.{pdf,png}
│   ├── pca.{pdf,png}
│   ├── volcano\_<contrast>.{pdf,png}
│   └── hallmark\_<direction>\_<contrast>\*.{pdf,png}
├── tables/
│   ├── filtered\_raw\_counts.csv
│   ├── normalized\_vst.csv
│   ├── pca\_scores.csv
│   ├── pca\_genes.csv
│   ├── de\_<contrast>\_annotated.csv
│   ├── hallmark\_<direction>\_<contrast>\_all.csv
│   └── hallmark\_<direction>\_<contrast>\_top10.csv
├── methods.txt
├── run\_summary.json
└── sample\_alignment.json
```

The volcano y-axis shows `−log10(p-value)`, while point coloring, labels, and significance calls use adjusted p-values. Pathway plots show `−log10(adjusted p-value)`, color bars by FDR status, and mark adjusted p = 0.05 with a dashed line.

## Statistical guardrails

* Differential expression always uses raw integer counts—not VST, TPM, CPM, or log-transformed values.
* TPM is used only for visualization and requires valid gene lengths.
* Paired analysis includes subject identity before condition in the design.
* Every contrast direction is explicit.
* A rank-deficient design stops with an error rather than silently dropping a covariate.
* QC flags possible outliers but does not remove samples automatically.

## Reproducibility

The workflow records package versions and all non-default configuration choices. Numerical results are deterministic under the same inputs, code, configuration, and package versions. PDF files can differ at the byte level because they contain creation timestamps, even when they render identically. GSEApy may also return pathway-member gene strings in a different order without changing pathway membership or statistics.

## Repository structure

* [`SKILL.md`](SKILL.md): Codex workflow and guardrails.
* [`scripts/run\_transcriptomics.py`](scripts/run_transcriptomics.py): configuration-driven runner.
* [`scripts/transcriptomics.py`](scripts/transcriptomics.py): reusable analysis functions.
* [`references/`](references/): configuration and method documentation.
* `ensemble\_mappings\_\*.pkl`: bundled human Ensembl annotations.

