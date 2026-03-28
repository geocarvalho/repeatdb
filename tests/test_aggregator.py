"""Testes da agregação por locus."""

from __future__ import annotations

import math

import polars as pl
import pytest

from repeatdb.aggregator import (
    build_locus_matrix,
    build_sample_locus_table,
    locus_stats,
)


def _synthetic_calls_df() -> pl.DataFrame:
    """DataFrame no formato do store (allele1/allele2 Int32 nullable)."""
    return pl.DataFrame(
        {
            "sample_id": ["S1", "S1", "S2", "S2", "S3", "S3"],
            "locus_id": ["HTT", "FMR1", "HTT", "FMR1", "HTT", "FMR1"],
            "allele1": [19, 30, 40, 55, None, 20],
            "allele2": [22, None, 35, None, None, 25],
            "filter": ["PASS"] * 6,
            "genotype": ["0/1"] * 6,
        },
        schema={
            "sample_id": pl.Utf8,
            "locus_id": pl.Utf8,
            "allele1": pl.Int32,
            "allele2": pl.Int32,
            "filter": pl.Utf8,
            "genotype": pl.Utf8,
        },
    )


def test_build_locus_matrix_excludes_no_calls_and_uses_max_allele() -> None:
    df = _synthetic_calls_df()
    matrix = build_locus_matrix(df)

    # S3 HTT: both null -> excluded from HTT list; max(19,22)=22 para S1
    assert "HTT" in matrix
    assert matrix["HTT"] == [22.0, 40.0]

    # FMR1: S1 max(30,null)=30, S2 55, S3 max(20,25)=25; sorted by sample_id
    assert matrix["FMR1"] == [30.0, 55.0, 25.0]


def test_build_locus_matrix_empty_after_all_no_calls() -> None:
    df = pl.DataFrame(
        {
            "sample_id": ["A"],
            "locus_id": ["X"],
            "allele1": [None],
            "allele2": [None],
            "filter": ["PASS"],
            "genotype": ["./."],
        },
        schema={
            "sample_id": pl.Utf8,
            "locus_id": pl.Utf8,
            "allele1": pl.Int32,
            "allele2": pl.Int32,
            "filter": pl.Utf8,
            "genotype": pl.Utf8,
        },
    )
    assert build_locus_matrix(df) == {}


def test_build_sample_locus_table_pivot_and_nulls() -> None:
    df = _synthetic_calls_df()
    wide = build_sample_locus_table(df)

    assert wide["sample_id"].to_list() == ["S1", "S2", "S3"]
    # S3 sem HTT (no-call) -> null na coluna HTT
    row_s3 = wide.filter(pl.col("sample_id") == "S3").row(0, named=True)
    htt = row_s3["HTT"]
    assert htt is None or (isinstance(htt, float) and math.isnan(htt))

    r1 = wide.filter(pl.col("sample_id") == "S1").row(0, named=True)
    assert r1["HTT"] == 22.0
    assert r1["FMR1"] == 30.0


def test_build_sample_locus_table_empty_input() -> None:
    df = pl.DataFrame(
        {
            "sample_id": pl.Series([], dtype=pl.Utf8),
            "locus_id": pl.Series([], dtype=pl.Utf8),
            "allele1": pl.Series([], dtype=pl.Int32),
            "allele2": pl.Series([], dtype=pl.Int32),
            "filter": pl.Series([], dtype=pl.Utf8),
            "genotype": pl.Series([], dtype=pl.Utf8),
        }
    )
    out = build_sample_locus_table(df)
    assert out.columns == ["sample_id"]
    assert out.height == 0


def test_locus_stats_known_values() -> None:
    vals = [1.0, 2.0, 3.0, 4.0, 5.0]
    s = locus_stats(vals)
    assert s["n"] == 5
    assert s["mean"] == pytest.approx(3.0)
    assert s["median"] == pytest.approx(3.0)
    assert s["q1"] == pytest.approx(2.0)
    assert s["q3"] == pytest.approx(4.0)
    assert s["iqr"] == pytest.approx(2.0)
    assert s["mad"] == pytest.approx(1.0)
    assert s["std"] == pytest.approx(math.sqrt(2.5))


def test_locus_stats_empty_and_single() -> None:
    empty = locus_stats([])
    assert empty["n"] == 0
    assert math.isnan(empty["mean"])

    one = locus_stats([42.0])
    assert one["n"] == 1
    assert one["mean"] == pytest.approx(42.0)
    assert math.isnan(one["std"])


def test_missing_columns_raises() -> None:
    df = pl.DataFrame({"sample_id": ["a"]})
    with pytest.raises(ValueError, match="em falta colunas"):
        build_locus_matrix(df)
