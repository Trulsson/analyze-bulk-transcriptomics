"""Deterministic modules for bulk RNA-seq count analysis.

Matrix convention inside this module is samples x genes. Raw counts are retained
for differential expression; normalized values are used only for QC and PCA.
"""

from __future__ import annotations

import json
import pickle
import re
from pathlib import Path
from typing import Any, Iterable

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

matplotlib.rc("pdf", fonttype=42)
sns.set_theme(style="whitegrid")


def _separator(path: str | Path) -> str:
    return "\t" if str(path).lower().endswith((".tsv", ".txt", ".out")) else ","


def _validate_unique(df: pd.DataFrame, axis: int, label: str) -> None:
    idx = df.index if axis == 0 else df.columns
    duplicated = idx[idx.duplicated()].unique().tolist()
    if duplicated:
        raise ValueError(f"Duplicate {label}: {duplicated[:10]}")


def validate_counts(counts: pd.DataFrame) -> pd.DataFrame:
    """Validate and return integer samples x genes counts."""
    counts = counts.apply(pd.to_numeric, errors="raise")
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Counts contain NaN or infinite values")
    if (values < 0).any():
        raise ValueError("Counts contain negative values")
    if not np.allclose(values, np.rint(values)):
        raise ValueError("Counts must be integer-valued raw counts")
    _validate_unique(counts, 0, "sample identifiers")
    _validate_unique(counts, 1, "gene identifiers")
    return counts.astype(np.int64)


def read_counts(
    path: str | Path,
    counts_format: str = "featurecounts",
    id_column: str | int | None = None,
    strip_ensembl_version: bool = False,
) -> tuple[pd.DataFrame, pd.Series | None]:
    """Read counts and optional gene lengths; return samples x genes."""
    path = Path(path)
    sep = _separator(path)
    lengths = None
    if counts_format == "featurecounts":
        table = pd.read_csv(path, sep="\t", comment="#")
        if "Geneid" not in table or "Length" not in table:
            raise ValueError("featureCounts input must contain Geneid and Length")
        genes = table["Geneid"].astype(str)
        lengths = pd.Series(table["Length"].to_numpy(float), index=genes, name="length_bp")
        annotation = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
        sample_cols = [c for c in table.columns if c not in annotation]
        counts = table.loc[:, sample_cols].copy()
        counts.columns = [
            re.sub(r"Aligned\.sortedByCoord\.out\.bam$", "", Path(str(c)).name)
            for c in counts.columns
        ]
        counts.index = genes
        counts = counts.T
    elif counts_format in {"genes_rows", "samples_rows"}:
        table = pd.read_csv(path, sep=sep)
        selected = table.columns[0] if id_column is None else id_column
        if isinstance(selected, int):
            selected = table.columns[selected]
        table[selected] = table[selected].astype(str)
        table = table.set_index(selected)
        counts = table.T if counts_format == "genes_rows" else table
    else:
        raise ValueError("counts_format must be featurecounts, genes_rows, or samples_rows")

    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    if strip_ensembl_version:
        counts.columns = counts.columns.str.replace(r"\.\d+$", "", regex=True)
        if lengths is not None:
            lengths.index = lengths.index.str.replace(r"\.\d+$", "", regex=True)
    return validate_counts(counts), lengths


def read_metadata(path: str | Path, sample_id_column: str) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep=_separator(path), dtype={sample_id_column: str})
    if sample_id_column not in metadata:
        raise ValueError(f"Metadata lacks sample ID column {sample_id_column!r}")
    metadata = metadata.set_index(sample_id_column, drop=True)
    metadata.index = metadata.index.astype(str)
    _validate_unique(metadata, 0, "metadata sample identifiers")
    return metadata


def align_samples(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    count_sample_id_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    if count_sample_id_column:
        if count_sample_id_column not in metadata:
            raise ValueError(f"Metadata lacks count sample ID column {count_sample_id_column!r}")
        count_ids = metadata[count_sample_id_column].astype(str)
        if count_ids.duplicated().any():
            duplicated = count_ids[count_ids.duplicated()].unique().tolist()
            raise ValueError(f"Duplicate count sample identifiers in metadata: {duplicated}")
        missing_counts = count_ids[~count_ids.isin(counts.index)].tolist()
        if missing_counts:
            raise ValueError(f"Metadata count IDs missing from counts: {missing_counts}")
        excluded_counts = counts.index.difference(count_ids).tolist()
        rename = dict(zip(count_ids, metadata.index))
        counts = counts.loc[count_ids].rename(index=rename)
    else:
        missing_counts = metadata.index.difference(counts.index).tolist()
        if missing_counts:
            raise ValueError(f"Metadata samples missing from counts: {missing_counts}")
        excluded_counts = counts.index.difference(metadata.index).tolist()
        counts = counts.loc[metadata.index]
    order = metadata.index
    report = {
        "included_samples": order.tolist(),
        "count_samples_excluded_absent_from_metadata": excluded_counts,
        "metadata_samples_missing_counts": [],
    }
    return counts.loc[order].copy(), metadata.loc[order].copy(), report


def filter_low_counts(
    counts: pd.DataFrame, count_threshold: int = 10, min_samples: int = 3
) -> tuple[pd.DataFrame, pd.Series]:
    if count_threshold < 0 or min_samples < 1:
        raise ValueError("count_threshold must be >=0 and min_samples >=1")
    keep = counts.ge(count_threshold).sum(axis=0).ge(min_samples)
    return counts.loc[:, keep].copy(), keep


def read_gene_annotation(
    path: str | Path, gene_id_column: str = "gene_id", value_column: str = "biotype"
) -> pd.Series:
    path = Path(path)
    if path.suffix.lower() in {".pkl", ".pickle"}:
        with path.open("rb") as handle:
            mapping = pickle.load(handle)
        if not isinstance(mapping, dict):
            raise ValueError(f"Pickle mapping must contain a dictionary: {path}")
        if not all(isinstance(key, str) for key in mapping):
            raise ValueError(f"Pickle mapping keys must be gene-ID strings: {path}")
        return pd.Series(mapping, name=value_column, dtype=object)
    table = pd.read_csv(path, sep=_separator(path), dtype={gene_id_column: str})
    missing = {gene_id_column, value_column}.difference(table.columns)
    if missing:
        raise ValueError(f"Annotation lacks columns: {sorted(missing)}")
    if table[gene_id_column].duplicated().any():
        raise ValueError("Gene annotation contains duplicate gene identifiers")
    return table.set_index(gene_id_column)[value_column]


def filter_biotype(
    counts: pd.DataFrame, biotypes: pd.Series | None, biotype: str = "protein_coding"
) -> tuple[pd.DataFrame, dict[str, int]]:
    if biotype == "all":
        return counts.copy(), {"unmapped_genes": 0, "retained_genes": counts.shape[1]}
    if biotypes is None:
        raise ValueError("Notebook-default biotype filtering requires biotype_path; use biotype='all' to opt out")
    mapped = biotypes.reindex(counts.columns)
    keep = mapped.eq(biotype)
    report = {
        "unmapped_genes": int(mapped.isna().sum()),
        "retained_genes": int(keep.sum()),
    }
    return counts.loc[:, keep].copy(), report


def build_design(
    condition_column: str,
    covariates: Iterable[str] = (),
    paired: bool = False,
    subject_column: str | None = None,
) -> str:
    terms: list[str] = []
    if paired:
        if not subject_column:
            raise ValueError("paired=true requires subject_column")
        terms.append(subject_column)
    for item in covariates:
        if item not in terms and item != condition_column:
            terms.append(item)
    terms.append(condition_column)
    return "~" + " + ".join(terms)


def validate_design(metadata: pd.DataFrame, design: str, contrasts: list[list[str]]) -> None:
    terms = [x.strip() for x in design.lstrip("~").split("+")]
    invalid = [x for x in terms if not re.fullmatch(r"[A-Za-z_]\w*", x)]
    if invalid:
        raise ValueError(f"Design columns must be simple Python-style names: {invalid}")
    missing = [x for x in terms if x not in metadata]
    if missing:
        raise ValueError(f"Metadata lacks design columns: {missing}")
    if metadata[terms].isna().any().any():
        bad = metadata.index[metadata[terms].isna().any(axis=1)].tolist()
        raise ValueError(f"Missing design metadata for samples: {bad}")
    for factor, test, reference in contrasts:
        if factor not in terms:
            raise ValueError(f"Contrast factor {factor!r} is not present in design {design!r}")
        levels = set(metadata[factor].astype(str))
        absent = {str(test), str(reference)}.difference(levels)
        if absent:
            raise ValueError(f"Contrast {factor} lacks levels: {sorted(absent)}")


def incomplete_pairs(
    metadata: pd.DataFrame, subject_column: str, factor: str, test: str, reference: str
) -> list[str]:
    required = {str(test), str(reference)}
    observed = metadata.groupby(subject_column)[factor].agg(lambda x: set(x.astype(str)))
    return observed.index[~observed.map(required.issubset)].astype(str).tolist()


def normalize_vst(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    n_cpus: int = 4,
    use_design: bool = False,
) -> tuple[pd.DataFrame, Any]:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
    except ImportError as exc:
        raise ImportError("VST requires pydeseq2==0.5.4") from exc
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        inference=DefaultInference(n_cpus=n_cpus),
        quiet=True,
    )
    dds.vst(use_design=use_design)
    normalized = pd.DataFrame(dds.layers["vst_counts"], index=counts.index, columns=counts.columns)
    return normalized, dds


def normalize_tpm(counts: pd.DataFrame, lengths_bp: pd.Series) -> pd.DataFrame:
    lengths = pd.to_numeric(lengths_bp.reindex(counts.columns), errors="coerce")
    if lengths.isna().any() or (lengths <= 0).any():
        bad = lengths.index[lengths.isna() | (lengths <= 0)].tolist()
        raise ValueError(f"Missing or invalid gene lengths for TPM: {bad[:10]}")
    rates = counts.div(lengths / 1000.0, axis=1)
    scale = rates.sum(axis=1) / 1_000_000.0
    if (scale == 0).any():
        raise ValueError("Cannot compute TPM for a zero-depth sample")
    return rates.div(scale, axis=0)


def compute_pca(
    normalized: pd.DataFrame, top_variable_genes: int | None = 500
) -> tuple[pd.DataFrame, PCA, list[str]]:
    variances = normalized.std(axis=0, ddof=1).sort_values(ascending=False)
    if top_variable_genes is None:
        genes = variances.index.tolist()
    else:
        if top_variable_genes < 2:
            raise ValueError("top_variable_genes must be >=2 or null")
        genes = variances.index[: min(top_variable_genes, len(variances))].tolist()
    variable = variances.loc[genes].gt(0)
    genes = variable.index[variable].tolist()
    if len(genes) < 2:
        raise ValueError("PCA requires at least two variable genes")
    scaled = StandardScaler().fit_transform(normalized.loc[:, genes])
    pca = PCA(n_components=2)
    scores = pca.fit_transform(scaled)
    return pd.DataFrame(scores, index=normalized.index, columns=["PC1", "PC2"]), pca, genes


def _palette_for(metadata: pd.DataFrame, color_column: str, palette: str | dict[str, str]):
    levels = metadata[color_column].astype(str).drop_duplicates().tolist()
    if isinstance(palette, dict):
        missing = set(levels).difference(palette)
        if missing:
            raise ValueError(f"Palette lacks colors for: {sorted(missing)}")
        return palette
    return dict(zip(levels, sns.color_palette(palette, n_colors=len(levels))))


def save_figure(fig: plt.Figure, stem: str | Path, formats: Iterable[str] = ("pdf", "png")) -> None:
    for ext in formats:
        fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def annotate_gene_table(table: pd.DataFrame, gene_names: pd.Series | None) -> pd.DataFrame:
    """Return a gene-indexed table with a non-unique-safe gene_symbol column."""
    annotated = table.copy()
    symbols = gene_names.reindex(annotated.index) if gene_names is not None else pd.Series(index=annotated.index, dtype=object)
    annotated.insert(0, "gene_symbol", symbols.replace("", np.nan))
    annotated.index.name = "gene_id"
    return annotated


def plot_qc(
    counts: pd.DataFrame,
    normalized: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: str | Path,
    color_column: str,
    palette: str | dict[str, str] = "colorblind",
    detected_threshold: int = 10,
    formats: Iterable[str] = ("pdf", "png"),
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    colors = _palette_for(metadata, color_column, palette)
    qc = metadata[[color_column]].copy()
    qc["library_size"] = counts.sum(axis=1)
    qc["detected_genes"] = counts.ge(detected_threshold).sum(axis=1)
    qc = qc.sort_values("library_size", ascending=False)

    for metric, ylabel, name in [
        ("library_size", "Total assigned reads", "library_sizes"),
        ("detected_genes", f"Genes with count >= {detected_threshold}", "detected_genes"),
    ]:
        fig, ax = plt.subplots(figsize=(max(7, 0.25 * len(qc)), 5))
        ax.bar(qc.index, qc[metric], color=qc[color_column].astype(str).map(colors))
        ax.set(xlabel="Sample", ylabel=ylabel)
        ax.tick_params(axis="x", rotation=90)
        save_figure(fig, output_dir / name, formats)

    corr = normalized.T.corr(method="pearson")
    size = max(6, 0.22 * len(corr))
    fig, ax = plt.subplots(figsize=(size, size))
    sns.heatmap(corr, cmap="vlag", center=0, vmin=-1, vmax=1, square=True, ax=ax)
    ax.set_title("Sample correlation")
    save_figure(fig, output_dir / "sample_correlation", formats)


def plot_pca(
    scores: pd.DataFrame,
    pca: PCA,
    metadata: pd.DataFrame,
    output_stem: str | Path,
    color_column: str,
    shape_column: str | None = None,
    label_column: str | None = None,
    palette: str | dict[str, str] = "colorblind",
    formats: Iterable[str] = ("pdf", "png"),
) -> None:
    data = scores.join(metadata)
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(
        data=data,
        x="PC1",
        y="PC2",
        hue=color_column,
        style=shape_column,
        palette=_palette_for(metadata, color_column, palette),
        s=50,
        alpha=0.8,
        ax=ax,
    )
    if label_column:
        for sample, row in data.iterrows():
            label = sample if label_column == "__index__" else row[label_column]
            ax.annotate(str(label), (row.PC1, row.PC2), xytext=(3, 3), textcoords="offset points", fontsize=8)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
    save_figure(fig, output_stem, formats)


def fit_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrasts: list[list[str]],
    n_cpus: int = 4,
    alpha: float = 0.05,
) -> tuple[Any, dict[str, pd.DataFrame]]:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise ImportError("Differential expression requires pydeseq2==0.5.4") from exc
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    dds.deseq2(fit_type="parametric")
    results: dict[str, pd.DataFrame] = {}
    for factor, test, reference in contrasts:
        stats = DeseqStats(
            dds,
            contrast=[factor, test, reference],
            alpha=alpha,
            cooks_filter=False,
            independent_filter=False,
            inference=inference,
            quiet=True,
        )
        stats.summary()
        key = f"{factor}__{test}_vs_{reference}"
        results[key] = stats.results_df.copy().sort_values("padj", na_position="last")
    return dds, results


def plot_volcano(
    results: pd.DataFrame,
    output_stem: str | Path,
    title: str,
    alpha: float = 0.05,
    label_lfc: float = 1.0,
    gene_names: pd.Series | None = None,
    formats: Iterable[str] = ("pdf", "png"),
) -> None:
    data = results.copy()
    if gene_names is not None:
        data["label"] = gene_names.reindex(data.index).fillna(pd.Series(data.index, index=data.index))
    else:
        data["label"] = data.index.astype(str)
    positive = data["pvalue"].dropna()
    floor = max(np.nextafter(0, 1), positive[positive.gt(0)].min() / 10 if positive.gt(0).any() else 1e-300)
    data["minus_log10_pvalue"] = -np.log10(data["pvalue"].clip(lower=floor))
    data["significant"] = np.where(data["padj"].lt(alpha), "yes", "no")
    cmap = {"yes": sns.color_palette("tab10")[-1], "no": sns.color_palette("pastel")[7]}
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(
        data=data,
        x="log2FoldChange",
        y="minus_log10_pvalue",
        hue="significant",
        hue_order=["yes", "no"],
        palette=cmap,
        alpha=0.6,
        linewidth=0,
        ax=ax,
    )
    labels = data[data["padj"].lt(alpha) & data["log2FoldChange"].abs().ge(label_lfc)]
    for _, row in labels.iterrows():
        ax.annotate(str(row["label"]), (row["log2FoldChange"], row["minus_log10_pvalue"]), fontsize=8)
    ax.set(title=title, ylabel="-log10(p-value)")
    save_figure(fig, output_stem, formats)


def run_pathway_enrichment(
    results: pd.DataFrame,
    gene_names: pd.Series,
    gene_sets: dict[str, list[str]],
    background_ids: Iterable[str],
    output_stem: str | Path,
    all_table_path: str | Path,
    top_table_path: str | Path,
    direction: str,
    alpha: float = 0.05,
    top_terms: int = 10,
    formats: Iterable[str] = ("pdf", "png"),
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Run symbol-based over-representation analysis for one DE direction."""
    try:
        import gseapy as gp
    except ImportError as exc:
        raise ImportError("Pathway enrichment requires gseapy") from exc

    named = annotate_gene_table(results, gene_names)
    significant = named["padj"].lt(alpha)
    directional = named["log2FoldChange"].gt(0) if direction == "up" else named["log2FoldChange"].lt(0)
    gene_list = named.loc[significant & directional, "gene_symbol"].dropna().drop_duplicates().tolist()
    background = gene_names.reindex(list(background_ids)).replace("", np.nan).dropna().drop_duplicates().tolist()
    if not gene_list:
        empty = pd.DataFrame(columns=["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Genes", "fdr_threshold", "significant"])
        empty.to_csv(all_table_path, index=False)
        empty.to_csv(top_table_path, index=False)
        return empty, {}

    enrichment = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism="human",
        background=background,
        outdir=None,
        no_plot=True,
    ).results.sort_values("Adjusted P-value").reset_index(drop=True)
    enrichment["Gene_ratio"] = enrichment["Overlap"].str.split("/").str[0].astype(int) / enrichment["Overlap"].str.split("/").str[1].astype(int)
    enrichment["minus_log10_padj"] = -np.log10(enrichment["Adjusted P-value"].clip(lower=np.nextafter(0, 1)))
    enrichment["fdr_threshold"] = alpha
    enrichment["significant"] = enrichment["Adjusted P-value"].lt(alpha)
    enrichment.to_csv(all_table_path, index=False)

    significant_pathways = enrichment[enrichment["Adjusted P-value"].lt(alpha)]
    top = enrichment.head(top_terms)
    top.to_csv(top_table_path, index=False)
    plotted = top.sort_values("minus_log10_padj")
    if not plotted.empty:
        fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(plotted))))
        sns.barplot(
            data=plotted,
            y="Term",
            x="minus_log10_padj",
            hue="significant",
            palette={True: "#D55E00", False: "#9E9E9E"},
            dodge=False,
            ax=ax,
        )
        ax.axvline(-np.log10(alpha), color="black", linestyle="--", linewidth=1.2, label=f"FDR = {alpha:g}")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, title="FDR significance")
        ax.set(title=f"Top {len(top)} Hallmark pathways: {direction}regulated", xlabel="-log10(adjusted p-value)", ylabel="")
        save_figure(fig, output_stem, formats)

        pathway_genes = {
            row.Term: set(str(row.Genes).split(";"))
            for row in plotted.itertuples()
        }
        long = []
        for term, symbols in pathway_genes.items():
            values = named.loc[significant & directional & named["gene_symbol"].isin(symbols), ["gene_symbol", "log2FoldChange"]]
            long.extend((term, symbol, lfc) for symbol, lfc in values.itertuples(index=False))
        if long:
            pathway_lfc = pd.DataFrame(long, columns=["Term", "gene_symbol", "log2FoldChange"])
            fig, ax = plt.subplots(figsize=(9, max(4, 0.4 * len(pathway_genes))))
            sns.swarmplot(data=pathway_lfc, y="Term", x="log2FoldChange", orient="h", size=4, ax=ax)
            ax.set(title=f"Genes in top Hallmark pathways: {direction}regulated", ylabel="")
            save_figure(fig, f"{output_stem}_gene_fold_changes", formats)

    pathway_map: dict[str, list[str]] = {}
    for row in significant_pathways.itertuples():
        for symbol in str(row.Genes).split(";"):
            pathway_map.setdefault(symbol, []).append(row.Term)
    return enrichment, pathway_map


def write_json(data: Any, path: str | Path) -> None:
    Path(path).write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")
