"""Agregação por locus para análise de outliers."""

from __future__ import annotations

import numpy as np
import polars as pl
from scipy import stats

_REQUIRED_COLS = frozenset({"sample_id", "locus_id", "allele1", "allele2"})


def _ensure_call_columns(df: pl.DataFrame) -> None:
    missing = _REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"DataFrame em falta colunas: {sorted(missing)}")


def _per_sample_locus_max(df: pl.DataFrame) -> pl.DataFrame:
    """
    Filtra no-calls, calcula max(allele1, allele2) por linha e agrega uma linha
    por (sample_id, locus_id) com o máximo entre duplicados.
    """
    _ensure_call_columns(df)
    return (
        df.filter(~(pl.col("allele1").is_null() & pl.col("allele2").is_null()))
        .with_columns(
            pl.max_horizontal(pl.col("allele1"), pl.col("allele2"))
            .cast(pl.Float64)
            .alias("max_allele")
        )
        .group_by(["sample_id", "locus_id"])
        .agg(pl.col("max_allele").max())
    )


def build_locus_matrix_with_sample_ids(
    df: pl.DataFrame,
) -> tuple[dict[str, list[float]], dict[str, list[str]]]:
    """
    Como :func:`build_locus_matrix`, mas devolve também, por locus, a lista de
    ``sample_id`` alinhada com cada valor (ordem por ``sample_id``).
    """
    prep = _per_sample_locus_max(df).sort(["locus_id", "sample_id"])
    if prep.is_empty():
        return {}, {}

    values: dict[str, list[float]] = {}
    samples: dict[str, list[str]] = {}
    grouped = prep.group_by("locus_id", maintain_order=True).agg(
        pl.col("sample_id"),
        pl.col("max_allele"),
    )
    for row in grouped.iter_rows(named=True):
        lid = str(row["locus_id"])
        values[lid] = [float(x) for x in row["max_allele"]]
        samples[lid] = [str(s) for s in row["sample_id"]]
    return values, samples


def build_locus_matrix(df: pl.DataFrame) -> dict[str, list[float]]:
    """
    Por cada ``locus_id``, devolve os máximos alelos por amostra (potencialmente
    patogénico), como lista de ``float``.

    Exclui no-calls (``allele1`` e ``allele2`` ambos nulos). Valores por locus
    estão ordenados por ``sample_id`` para estabilidade.
    """
    vals, _ = build_locus_matrix_with_sample_ids(df)
    return vals


def build_sample_locus_table(df: pl.DataFrame) -> pl.DataFrame:
    """
    Tabela larga: linhas = amostras, colunas = loci, valores = max allele.
    Combinações sem chamada ficam nulas.
    """
    prep = _per_sample_locus_max(df)
    if prep.is_empty():
        return pl.DataFrame({"sample_id": pl.Series([], dtype=pl.Utf8)})

    wide = prep.pivot(
        on="locus_id",
        index="sample_id",
        values="max_allele",
    )
    return wide.sort("sample_id")


def locus_stats(values: list[float]) -> dict[str, float | int]:
    """
    Estatísticas descritivas para um locus (lista de máximos alelos por amostra).

    Inclui ``n``, ``mean``, ``median``, ``std`` (desvio-padrão amostral, ``ddof=1``),
    ``q1``, ``q3``, ``iqr`` (``scipy.stats.iqr``) e ``mad`` (desvio absoluto mediano
    em torno da mediana, ``scipy.stats.median_abs_deviation`` com ``scale=1``).
    """
    arr = np.asarray(values, dtype=np.float64)
    n = int(arr.size)
    if n == 0:
        nan = float("nan")
        return {
            "n": 0,
            "mean": nan,
            "median": nan,
            "std": nan,
            "q1": nan,
            "q3": nan,
            "iqr": nan,
            "mad": nan,
        }

    mean = float(np.mean(arr))
    median = float(np.median(arr))
    std = float(np.std(arr, ddof=1)) if n > 1 else float("nan")
    q1 = float(np.percentile(arr, 25.0))
    q3 = float(np.percentile(arr, 75.0))
    iqr = float(stats.iqr(arr))
    mad = float(stats.median_abs_deviation(arr, scale=1, nan_policy="omit"))

    return {
        "n": n,
        "mean": mean,
        "median": median,
        "std": std,
        "q1": q1,
        "q3": q3,
        "iqr": iqr,
        "mad": mad,
    }
