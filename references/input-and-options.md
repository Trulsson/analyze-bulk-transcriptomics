# Inputs and options

## Counts

Use raw gene-level integer counts.

- `featurecounts`: tab-separated featureCounts output. Ignore `#` lines, use `Geneid`, retain `Length`, drop the first six annotation columns, and simplify BAM-derived sample names.
- `genes_rows`: CSV/TSV with genes as rows and samples as columns; use the first column as gene identifier unless configured otherwise.
- `samples_rows`: CSV/TSV with samples as rows and genes as columns; use the first column as sample identifier.

Counts must be finite, non-negative, integer-valued, and have unique gene and sample identifiers. Sum duplicated genes upstream only when their biological meaning is clear; otherwise stop.

## Metadata

Supply one row per sample. `sample_id_column` becomes the analysis index. Every metadata sample must exist in counts. If count-matrix columns use sequencing IDs while analysis labels use biopsy or patient IDs, set `count_sample_id_column` to the metadata column holding the sequencing IDs; the code performs the notebook's explicit mapping. Exclude count columns absent from metadata and list them in `sample_alignment.json`.

The default condition column is `condition`. Add categorical covariates in `covariates`. For pairing, set `paired: true` and provide `subject_column`; the design becomes `~subject + covariate1 + condition`.

## Gene annotation

The notebook default keeps only `protein_coding` genes. The skill bundles Ensembl gene-ID-to-biotype and gene-ID-to-symbol dictionaries; `biotype_path: "bundled"` and `gene_name_path: "bundled"` use them automatically. Supply a trusted `.pkl` dictionary or CSV/TSV with `gene_id` plus the configured value column to override either mapping. Never load an untrusted pickle. Set `biotype: "all"` explicitly to skip biotype filtering. FeatureCounts `Length` supplies gene lengths for TPM; otherwise use `gene_length_path` with `gene_id` and `length_bp`.

Refresh both bundled mappings with `python scripts/update_ensembl_mappings.py`. This queries the human Ensembl BioMart dataset and replaces each mapping only after a non-empty response has been parsed. Install the updater-only dependency with `python -m pip install biomart`.

Strip Ensembl version suffixes only when annotation identifiers are unversioned. Record this transformation because versioned IDs are distinct strings.

## Notebook-derived defaults

| Option | Default |
| --- | --- |
| Minimum count | 10 |
| Minimum samples meeting count | 3 |
| Biotype | `protein_coding` |
| Visualization normalization | `vst` |
| VST design-aware dispersion fit | `false` |
| PCA gene set | 500 most variable genes |
| PCA scaling | Standardize each gene across samples |
| DE fit | PyDESeq2 parametric |
| Cook's refitting | `true` |
| Contrast Cook's filtering | `false` |
| Independent filtering | `false` |
| FDR threshold | 0.05 |
| Volcano label LFC threshold | 1.0 |
| Hallmark enrichment | Separate up/down sets, MSigDB `h.all`, release `2024.1.Hs` |
| Enrichment significance | Adjusted p-value < 0.05 |
| Enrichment plot and focused CSV | Top 10 terms per direction, regardless of significance |
| Figure formats | PDF and PNG |

## Useful overrides

- `normalization: "tpm"`: use log2(TPM + 1) for QC/PCA and raw counts for DE.
- `paired: true`, `subject_column: "patient"`: account for within-subject pairing.
- `covariates: ["sex", "batch"]`: adjust the DE design. Include only justified, non-confounded variables.
- `contrasts`: run several explicit comparisons without changing preprocessing.
- `palette`: a seaborn palette name such as `colorblind`, or a mapping such as `{"treated":"#D55E00","control":"#0072B2"}`.
- `pca_color_column`, `pca_shape_column`, `sample_label_column`: alter PCA encoding without changing the fit.
- `top_variable_genes: null`: run PCA using every retained gene, matching the notebook's alternate all-gene PCA.
- `run_de: false`: generate filtering, normalization, and QC only.
- `run_enrichment: false`: skip Hallmark over-representation analysis. Enrichment otherwise runs separately for significant positive- and negative-LFC genes, using all retained mapped protein-coding genes as background.
- `enrichment_category`, `enrichment_dbver`, `enrichment_top_terms`: select the MSigDB collection/release and the number of ranked terms in the focused CSV and plot without truncating the complete `_all.csv` table. Both tables include the FDR threshold and a significance flag.

## Installation

Use Python 3.11-3.13 and install:

```bash
python -m pip install "pydeseq2==0.5.4" pandas numpy scipy scikit-learn matplotlib seaborn gseapy
```

PyDESeq2 is a Python reimplementation and can differ slightly from Bioconductor DESeq2. Pin the environment and retain package versions with results.
