"""Ingestão em lote, Parquet e carregamento de dados."""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import polars as pl
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
)

from repeatdb.vcf_parser import parse_vcf

logger = logging.getLogger(__name__)

_COLUMNS_SCHEMA: dict[str, pl.DataType] = {
    "sample_id": pl.Utf8,
    "locus_id": pl.Utf8,
    "allele1": pl.Int32,
    "allele2": pl.Int32,
    "filter": pl.Utf8,
    "genotype": pl.Utf8,
}


def _parse_vcf_rows(vcf_path: str) -> list[dict[str, Any]]:
    """Worker picklável: converte um VCF em linhas dict para ``DataFrame``."""
    return [c.model_dump() for c in parse_vcf(vcf_path)]


def _empty_calls_frame() -> pl.DataFrame:
    return pl.DataFrame(schema=_COLUMNS_SCHEMA)


def ingest_vcf_list(vcf_paths: list[str], n_workers: int = 4) -> pl.DataFrame:
    """
    Ingere vários VCFs em paralelo e devolve um único ``DataFrame`` polars.

    Colunas: ``sample_id``, ``locus_id``, ``allele1``, ``allele2``, ``filter``,
    ``genotype``. Tipos: ``allele1``/``allele2`` como ``Int32`` nullable;
    restantes ``Utf8``.
    """
    if not vcf_paths:
        return _empty_calls_frame()

    workers = max(1, n_workers)
    rows: list[dict[str, Any]] = []

    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeElapsedColumn(),
    ) as progress:
        task_id = progress.add_task("A ingestar VCFs", total=len(vcf_paths))
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_path = {
                executor.submit(_parse_vcf_rows, path): path for path in vcf_paths
            }
            for future in as_completed(future_to_path):
                path = future_to_path[future]
                try:
                    rows.extend(future.result())
                except Exception as exc:
                    logger.warning("Falha ao processar %s: %s", path, exc)
                progress.advance(task_id)

    if not rows:
        return _empty_calls_frame()

    return pl.DataFrame(rows, schema=_COLUMNS_SCHEMA)


def save_parquet(df: pl.DataFrame, output_path: str) -> None:
    """Grava ``df`` em Parquet com compressão Snappy."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(output_path, compression="snappy")


def load_parquet(input_path: str) -> pl.DataFrame:
    """Lê um Parquet e devolve o ``DataFrame``."""
    return pl.read_parquet(input_path)
