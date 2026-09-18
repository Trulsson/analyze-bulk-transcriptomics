#!/usr/bin/env python3
"""Run the notebook-derived bulk transcriptomics workflow from JSON config."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd

from transcriptomics import (
    align_samples,
    annotate_gene_table,
    build_design,
    compute_pca,
    filter_biotype,
    filter_low_counts,
    fit_deseq2,
    incomplete_pairs,
    normalize_tpm,
    normalize_vst,
    plot_pca,
    plot_qc,
    plot_volcano,
    read_counts,
    read_gene_annotation,
    read_metadata,
    run_pathway_enrichment,
    validate_design,
    write_json,
)

SKILL_ROOT = Path(__file__).resolve().parent.parent
BUNDLED_BIOTYPES = SKILL_ROOT / "ensemble_mappings_biotype.pkl"
BUNDLED_GENE_NAMES = SKILL_ROOT / "ensemble_mappings_gene.pkl"

DEFAULTS = {
    "counts_format": "featurecounts",
    "sample_id_column": "sample_id",
    "count_sample_id_column": None,
    "condition_column": "condition",
    "biotype": "protein_coding",
    "biotype_path": "bundled",
    "gene_name_path": "bundled",
    "count_threshold": 10,
    "min_samples": 3,
    "normalization": "vst",
    "vst_use_design": False,
    "top_variable_genes": 500,
    "paired": False,
    "subject_column": None,
    "covariates": [],
    "contrasts": [],
    "palette": "colorblind",
    "pca_color_column": None,
    "pca_shape_column": None,
    "sample_label_column": None,
    "run_de": True,
    "run_enrichment": True,
    "enrichment_category": "h.all",
    "enrichment_dbver": "2024.1.Hs",
    "enrichment_top_terms": 10,
    "n_cpus": 4,
    "alpha": 0.05,
    "figure_formats": ["pdf", "png"],
    "strip_ensembl_version": False,
    "output_dir": "transcriptomics_results",
}


def mapping_path(value: str | None, bundled: Path) -> Path | None:
    if value == "bundled":
        if not bundled.exists():
            raise FileNotFoundError(f"Bundled mapping is missing: {bundled}")
        return bundled
    return Path(value) if value else None


def load_config(path: str | Path) -> dict:
    supplied = json.loads(Path(path).read_text(encoding="utf-8"))
    config = DEFAULTS | supplied
    required = ["counts_path", "metadata_path"]
    missing = [x for x in required if not config.get(x)]
    if missing:
        raise ValueError(f"Configuration lacks required fields: {missing}")
    if not config["contrasts"] and config["run_de"]:
        raise ValueError("run_de=true requires at least one explicit contrast")
    config["pca_color_column"] = config["pca_color_column"] or config["condition_column"]
    return config


def read_lengths(config: dict, embedded: pd.Series | None) -> pd.Series | None:
    path = config.get("gene_length_path")
    if not path:
        return embedded
    return read_gene_annotation(
        path,
        config.get("gene_length_gene_id_column", "gene_id"),
        config.get("gene_length_column", "length_bp"),
    ).astype(float)


def methods_text(config: dict, design: str, retained: int) -> str:
    norm = (
        "DESeq2 variance-stabilizing transformation with use_design=" + str(config["vst_use_design"])
        if config["normalization"] == "vst"
        else "log2(TPM + 1), with raw counts retained for DESeq2"
    )
    top = "all retained genes" if config["top_variable_genes"] is None else f"the {config['top_variable_genes']} most variable genes"
    return (
        f"Genes were retained when they had at least {config['count_threshold']} raw reads in at least "
        f"{config['min_samples']} samples, followed by biotype filtering to {config['biotype']}. "
        f"This retained {retained} genes. QC and PCA used {norm}. PCA used {top}, with each gene "
        f"standardized across samples before fitting. Differential expression used filtered raw counts "
        f"with PyDESeq2 0.5.4 and design {design}; Cook's refitting was enabled, while Cook's filtering "
        f"and independent filtering were disabled for contrasts. The volcano y-axis used nominal p-values, "
        f"while significance and coloring used adjusted p-values. The FDR threshold was {config['alpha']}. "
        f"Gene-level tables used Ensembl IDs as the index and gene symbols as a separate column. "
        f"Separate up/down Hallmark enrichment used all retained mapped genes as background"
        f"{' and MSigDB ' + config['enrichment_dbver'] if config['run_enrichment'] else ' (disabled)'}. "
        f"The top {config['enrichment_top_terms']} pathways per direction were exported and plotted "
        f"regardless of significance, with the adjusted-p = {config['alpha']} threshold marked.\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="JSON configuration file")
    args = parser.parse_args()
    config = load_config(args.config)
    out = Path(config["output_dir"])
    figures = out / "figures"
    tables = out / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    counts, embedded_lengths = read_counts(
        config["counts_path"],
        config["counts_format"],
        config.get("count_id_column"),
        config["strip_ensembl_version"],
    )
    metadata = read_metadata(config["metadata_path"], config["sample_id_column"])
    counts, metadata, alignment = align_samples(counts, metadata, config["count_sample_id_column"])
    qc_counts = counts.copy()
    write_json(alignment, out / "sample_alignment.json")
    raw_gene_count = counts.shape[1]
    counts, keep = filter_low_counts(counts, config["count_threshold"], config["min_samples"])

    biotypes = None
    biotype_path = (
        mapping_path(config.get("biotype_path"), BUNDLED_BIOTYPES)
        if config["biotype"] != "all"
        else None
    )
    if biotype_path:
        biotypes = read_gene_annotation(
            biotype_path,
            config.get("biotype_gene_id_column", "gene_id"),
            config.get("biotype_column", "biotype"),
        )
    counts, biotype_report = filter_biotype(counts, biotypes, config["biotype"])
    if counts.shape[1] < 2:
        raise ValueError("Fewer than two genes remain after filtering")
    gene_name_path = mapping_path(config.get("gene_name_path"), BUNDLED_GENE_NAMES)
    gene_names = (
        read_gene_annotation(
            gene_name_path,
            config.get("gene_name_gene_id_column", "gene_id"),
            config.get("gene_name_column", "gene_name"),
        )
        if gene_name_path else None
    )
    annotate_gene_table(counts.T, gene_names).to_csv(tables / "filtered_raw_counts.csv")

    design = build_design(
        config["condition_column"], config["covariates"], config["paired"], config["subject_column"]
    )
    validate_design(metadata, design, config["contrasts"])
    pair_report = {}
    if config["paired"]:
        for factor, test, reference in config["contrasts"]:
            pair_report[f"{factor}__{test}_vs_{reference}"] = incomplete_pairs(
                metadata, config["subject_column"], factor, test, reference
            )

    if config["normalization"] == "vst":
        normalized, _ = normalize_vst(
            counts, metadata, design, config["n_cpus"], config["vst_use_design"]
        )
    elif config["normalization"] == "tpm":
        lengths = read_lengths(config, embedded_lengths)
        if lengths is None:
            raise ValueError("TPM normalization requires featureCounts Length or gene_length_path")
        normalized = np.log2(normalize_tpm(counts, lengths) + 1)
    else:
        raise ValueError("normalization must be 'vst' or 'tpm'")
    annotate_gene_table(normalized.T, gene_names).to_csv(tables / f"normalized_{config['normalization']}.csv")

    scores, pca, pca_genes = compute_pca(normalized, config["top_variable_genes"])
    scores.to_csv(tables / "pca_scores.csv")
    annotate_gene_table(pd.DataFrame(index=pca_genes), gene_names).to_csv(tables / "pca_genes.csv")
    plot_qc(
        qc_counts,
        normalized,
        metadata,
        figures,
        config["pca_color_column"],
        config["palette"],
        config["count_threshold"],
        config["figure_formats"],
    )
    plot_pca(
        scores,
        pca,
        metadata,
        figures / "pca",
        config["pca_color_column"],
        config["pca_shape_column"],
        config["sample_label_column"],
        config["palette"],
        config["figure_formats"],
    )

    de_summary = {}
    enrichment_summary = {}
    if config["run_de"]:
        _, de_results = fit_deseq2(
            counts, metadata, design, config["contrasts"], config["n_cpus"], config["alpha"]
        )
        gene_sets = None
        if config["run_enrichment"]:
            if gene_names is None:
                raise ValueError("Pathway enrichment requires gene_name_path")
            try:
                from gseapy import Msigdb
            except ImportError as exc:
                raise ImportError("Pathway enrichment requires gseapy") from exc
            gene_sets = Msigdb().get_gmt(
                category=config["enrichment_category"], dbver=config["enrichment_dbver"]
            )
            if not gene_sets:
                raise RuntimeError("MSigDB returned no pathway gene sets")
        for key, result in de_results.items():
            annotated = annotate_gene_table(result, gene_names)
            de_summary[key] = int(result["padj"].lt(config["alpha"]).sum())
            plot_volcano(
                result,
                figures / f"volcano_{key}",
                key.replace("__", ": ").replace("_vs_", " vs "),
                config["alpha"],
                config.get("volcano_label_lfc", 1.0),
                gene_names,
                config["figure_formats"],
            )
            if gene_sets is not None:
                pathway_maps = []
                for direction in ("up", "down"):
                    enrichment, pathway_map = run_pathway_enrichment(
                        result,
                        gene_names,
                        gene_sets,
                        counts.columns,
                        figures / f"hallmark_{direction}_{key}",
                        tables / f"hallmark_{direction}_{key}_all.csv",
                        tables / f"hallmark_{direction}_{key}_top{config['enrichment_top_terms']}.csv",
                        direction,
                        config["alpha"],
                        config["enrichment_top_terms"],
                        config["figure_formats"],
                    )
                    enrichment_summary[f"{key}__{direction}"] = {
                        "tested_pathways": len(enrichment),
                        "significant_pathways": int(enrichment.get("Adjusted P-value", pd.Series(dtype=float)).lt(config["alpha"]).sum()),
                    }
                    pathway_maps.append(pathway_map)
                combined = {symbol: sorted(set(pathway_maps[0].get(symbol, []) + pathway_maps[1].get(symbol, []))) for symbol in set(pathway_maps[0]) | set(pathway_maps[1])}
                annotated["Hallmark"] = annotated["gene_symbol"].map(combined).map(
                    lambda terms: "; ".join(terms) if isinstance(terms, list) else np.nan
                )
            annotated.to_csv(tables / f"de_{key}_annotated.csv")

    versions = {}
    for package in ["pydeseq2", "pandas", "numpy", "scikit-learn", "matplotlib", "seaborn", "gseapy"]:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    supplied = json.loads(Path(args.config).read_text(encoding="utf-8"))
    summary = {
        "samples": counts.shape[0],
        "genes_input": raw_gene_count,
        "genes_after_low_count_filter": int(keep.sum()),
        "genes_after_biotype_filter": counts.shape[1],
        "biotype_report": biotype_report,
        "design": design,
        "normalization": config["normalization"],
        "pca_genes": len(pca_genes),
        "pca_explained_variance_percent": (pca.explained_variance_ratio_ * 100).tolist(),
        "incomplete_pairs": pair_report,
        "significant_genes_by_contrast": de_summary,
        "pathway_enrichment_by_contrast_and_direction": enrichment_summary,
        "overrides_from_defaults": {k: v for k, v in supplied.items() if k in DEFAULTS and v != DEFAULTS[k]},
        "package_versions": versions,
    }
    write_json(summary, out / "run_summary.json")
    (out / "methods.txt").write_text(methods_text(config, design, counts.shape[1]), encoding="utf-8")


if __name__ == "__main__":
    main()
