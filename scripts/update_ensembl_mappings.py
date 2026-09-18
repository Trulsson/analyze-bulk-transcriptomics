#!/usr/bin/env python3
"""Refresh the bundled human Ensembl gene-name and biotype mappings."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path


def parse_biomart_tsv(data: str) -> tuple[dict[str, str], dict[str, str]]:
    gene_names: dict[str, str] = {}
    biotypes: dict[str, str] = {}
    for row in data.splitlines():
        fields = row.split("\t")
        if len(fields) != 3 or not fields[0]:
            continue
        gene_id, gene_name, biotype = fields
        gene_names[gene_id] = gene_name
        biotypes[gene_id] = biotype
    return gene_names, biotypes


def fetch_mappings(host: str) -> tuple[dict[str, str], dict[str, str]]:
    try:
        import biomart
    except ImportError as exc:
        raise SystemExit("Install the updater dependency with: python -m pip install biomart") from exc
    server = biomart.BiomartServer(host)
    mart = server.datasets["hsapiens_gene_ensembl"]
    response = mart.search(
        {"attributes": ["ensembl_gene_id", "external_gene_name", "gene_biotype"]}
    )
    return parse_biomart_tsv(response.raw.data.decode("utf-8"))


def write_pickle(mapping: dict[str, str], path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(mapping, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Directory for the two mapping pickle files (default: skill root)",
    )
    parser.add_argument("--host", default="https://www.ensembl.org/biomart")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    gene_names, biotypes = fetch_mappings(args.host)
    if not gene_names or not biotypes:
        raise SystemExit("BioMart returned empty mappings; existing files were not replaced")
    write_pickle(gene_names, args.output_dir / "ensemble_mappings_gene.pkl")
    write_pickle(biotypes, args.output_dir / "ensemble_mappings_biotype.pkl")
    print(f"Saved {len(gene_names):,} gene names and {len(biotypes):,} biotypes")


if __name__ == "__main__":
    main()
