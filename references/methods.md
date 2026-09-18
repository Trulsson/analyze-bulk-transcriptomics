# Methods and interpretation

## Default preprocessing

Align samples by identifier using metadata as the analysis set. Retain genes with at least 10 raw reads in at least 3 samples. Map gene biotypes and retain protein-coding genes. These are the exact filtering choices used in the source notebook.

## Normalization and PCA

Estimate DESeq2 size factors and apply the variance-stabilizing transformation with a design-blind dispersion fit (`use_design=False`). For PCA, select the 500 genes with the largest sample standard deviation in VST space, standardize each selected gene to mean 0 and variance 1 across samples, and fit two principal components. Report variance explained on each axis.

When TPM is requested, divide each gene's count by length in kilobases, scale each sample's rates to one million, and use log2(TPM + 1) for PCA and correlation figures. Do not use TPM for differential expression.

## Differential expression

Fit a negative-binomial model to filtered raw counts using PyDESeq2 0.5.4. Use a design formula ending in the condition term. For paired data include subject identity before condition; for unpaired data include only requested covariates before condition. Enable Cook's outlier refitting in the dataset fit. Adjust p-values by the method implemented by PyDESeq2 and call genes significant at adjusted p < 0.05 unless overridden.

Positive log2 fold changes indicate greater expression in the named test level relative to the reference level.

## QC figures

- Library-size plot: total raw assigned reads per sample, ordered by depth and colored by requested metadata group.
- Detected-gene plot: genes with raw count >=10 in each sample.
- Correlation heatmap: sample-to-sample Pearson correlation in normalized expression space.
- PCA: condition-colored scatter plot with optional shapes and sample labels.
- Volcano plot: log2 fold change against -log10 nominal p-value; significance, coloring, and labeling remain based on adjusted p < the configured FDR threshold. Cap zero nominal p-values only for display and retain original values in the table.

## Output tables and pathway enrichment

Write every gene-level table with Ensembl gene IDs as the row index and `gene_symbol` as a separate column. Never use symbols as the index because several Ensembl IDs can map to one symbol. Transpose count and normalized-expression matrices for export so genes are rows and samples are columns. Write only the symbol-annotated DE table.

Run over-representation analysis separately for significant positive- and negative-log2-fold-change gene sets. Use the mapped symbols of all retained protein-coding genes as the background. By default use the MSigDB Hallmark collection (`h.all`, release `2024.1.Hs`) and call pathways significant at adjusted p < 0.05. For each direction, export every tested pathway to `_all.csv` and the ten best-ranked pathways to `_top10.csv`, including `fdr_threshold` and `significant` columns. Always plot the top ten pathways, color bars by FDR status, and mark adjusted p = 0.05 with a dashed line so non-significant trends remain visible for hypothesis generation. Plot the significant directional genes contributing to those top pathways. Add significant Hallmark membership back to the annotated DE table without changing its Ensembl index.

QC highlights potential problems but does not define automatic sample exclusions. Investigate possible outliers against sequencing and experimental metadata before removal.

## Source provenance

Defaults and visual conventions were extracted from `pyDEseq_transcriptome_baseline.ipynb` in `Trulsson/Transcriptomic_profiles_in_CHC_predict_HCC` at Git blob `2b9e81a872fb1d4e4c2bf36b72de501c4df113d4` (accessed 2026-08-24). The modular implementation targets the current PyDESeq2 formula API while retaining the notebook's statistical choices.
