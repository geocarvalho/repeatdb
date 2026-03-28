"""Testes de detecção de outliers por locus."""

from __future__ import annotations

import logging

import polars as pl
import pytest

from repeatdb.outlier import OutlierDetector

# Dez inliers ligeiramente variados + um outlier claro (n = 10, limite mínimo)
_BASE_VALS = [38.0, 39.0, 40.0, 41.0, 40.0, 39.0, 40.0, 41.0, 40.0, 100.0]
_BASE_IDS = [f"S{i}" for i in range(len(_BASE_VALS))]


def _matrix_l1() -> dict[str, list[float]]:
    return {"L1": list(_BASE_VALS)}


def _ids_l1() -> dict[str, list[str]]:
    return {"L1": list(_BASE_IDS)}


def test_mad_flags_known_outlier() -> None:
    det = OutlierDetector(sample_ids_by_locus=_ids_l1())
    rows = det.detect(_matrix_l1(), method="mad", threshold=3.5)
    outliers = [r for r in rows if r["is_outlier"]]
    assert len(outliers) == 1
    assert outliers[0]["sample_id"] == "S9"
    assert outliers[0]["locus_id"] == "L1"
    assert outliers[0]["allele_value"] == 100.0
    assert abs(outliers[0]["score"]) > 3.5


def test_iqr_flags_known_outlier() -> None:
    det = OutlierDetector(sample_ids_by_locus=_ids_l1())
    rows = det.detect(_matrix_l1(), method="iqr", threshold=1.5)
    outliers = [r for r in rows if r["is_outlier"]]
    assert len(outliers) == 1
    assert outliers[0]["sample_id"] == "S9"


def test_zscore_flags_known_outlier() -> None:
    det = OutlierDetector(sample_ids_by_locus=_ids_l1())
    # z ≈ 2,84 para o valor 100 com estes inliers; limiar < 3
    rows = det.detect(_matrix_l1(), method="zscore", threshold=2.5)
    outliers = [r for r in rows if r["is_outlier"]]
    assert len(outliers) == 1
    assert outliers[0]["sample_id"] == "S9"


def test_default_method_is_mad() -> None:
    det = OutlierDetector(sample_ids_by_locus=_ids_l1())
    rows_default = det.detect(_matrix_l1())
    rows_mad = det.detect(_matrix_l1(), method="mad", threshold=3.5)
    assert rows_default == rows_mad


def test_skips_locus_with_few_samples(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.WARNING)
    det = OutlierDetector(sample_ids_by_locus={"L2": [f"A{i}" for i in range(5)]})
    rows = det.detect({"L2": [1.0, 2.0, 3.0, 4.0, 5.0]}, method="mad", threshold=3.5)
    assert rows == []
    assert any("apenas 5 amostras" in r.message for r in caplog.records)


def test_invalid_method_raises() -> None:
    det = OutlierDetector()
    with pytest.raises(ValueError, match="método inválido"):
        det.detect({"L": [1.0] * 10}, method="robust", threshold=1.0)


def test_results_roundtrip_polars() -> None:
    det = OutlierDetector(sample_ids_by_locus=_ids_l1())
    rows = det.detect(_matrix_l1(), method="mad", threshold=3.5)
    df = pl.DataFrame(rows)
    assert df.columns == [
        "sample_id",
        "locus_id",
        "allele_value",
        "score",
        "is_outlier",
    ]
    assert df.height == 10
