"""Detecção de outliers em repeat counts por locus."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from scipy import stats

logger = logging.getLogger(__name__)

_METHODS = frozenset({"zscore", "iqr", "mad"})


class OutlierDetector:
    """
    Deteta outliers por locus a partir de uma matriz ``locus_id → [valores]``.

    Os valores devem estar na mesma ordem que ``sample_ids_by_locus`` para cada
    locus (por exemplo, saída de :func:`repeatdb.aggregator.build_locus_matrix`
    alinhada com os ``sample_id`` ordenados). Se ``sample_ids_by_locus`` for
    ``None``, usam-se identificadores ``"0"``, ``"1"``, … por locus.
    """

    def __init__(self, sample_ids_by_locus: dict[str, list[str]] | None = None) -> None:
        self._sample_ids_by_locus = sample_ids_by_locus

    def detect(
        self,
        locus_matrix: dict[str, list[float]],
        method: str = "mad",
        threshold: float = 3.5,
    ) -> list[dict[str, Any]]:
        """
        Devolve uma linha por (amostra, locus) com ``allele_value``, ``score`` e
        ``is_outlier``.

        Métodos:

        - ``"zscore"``: :math:`z=(x-\\mu)/\\sigma` (``\\sigma`` amostral); outlier
          se :math:`|z| > \\text{threshold}`. Limiar típico: 3.0.
        - ``"iqr"``: cercas de Tukey :math:`Q_1 - k\\cdot\\mathrm{IQR}` e
          :math:`Q_3 + k\\cdot\\mathrm{IQR}` com ``k = threshold``; típico
          ``threshold=1.5``.
        - ``"mad"``: :math:`z_{\\mathrm{mad}} = 0{,}6745\\,(x-\\mathrm{med})/\\mathrm{MAD}`
          com ``MAD`` = :func:`scipy.stats.median_abs_deviation` e
          ``scale="normal"``; outlier se :math:`|z_{\\mathrm{mad}}| > \\text{threshold}`.
          Limiar típico: 3.5 (Leys et al., 2013).

        Loci com menos de 10 valores são ignorados (aviso em log).
        """
        if method not in _METHODS:
            raise ValueError(
                f"método inválido {method!r}; use um de {sorted(_METHODS)}"
            )

        out: list[dict[str, Any]] = []

        for locus_id, values in locus_matrix.items():
            vals = list(values)
            n = len(vals)
            if n < 10:
                logger.warning(
                    "Locus %r tem apenas %s amostras com call (< 10); ignorado",
                    locus_id,
                    n,
                )
                continue

            sample_ids, vals = self._align_samples(locus_id, vals)
            arr = np.asarray(vals, dtype=np.float64)
            if arr.size < 10:
                logger.warning(
                    "Locus %r tem apenas %s amostras alinhadas (< 10); ignorado",
                    locus_id,
                    arr.size,
                )
                continue

            if method == "zscore":
                out.extend(self._detect_zscore(locus_id, sample_ids, arr, threshold))
            elif method == "iqr":
                out.extend(self._detect_iqr(locus_id, sample_ids, arr, threshold))
            else:
                out.extend(self._detect_mad(locus_id, sample_ids, arr, threshold))

        return out

    def _align_samples(self, locus_id: str, values: list[float]) -> tuple[list[str], list[float]]:
        n = len(values)
        if self._sample_ids_by_locus is None:
            return [str(i) for i in range(n)], values

        ids = self._sample_ids_by_locus.get(locus_id)
        if ids is None:
            logger.warning(
                "Sem sample_ids para locus %r; a usar índices 0..%s",
                locus_id,
                n - 1,
            )
            return [str(i) for i in range(n)], values

        if len(ids) != n:
            logger.warning(
                "locus %r: %s sample_ids para %s valores; truncando ao mínimo",
                locus_id,
                len(ids),
                n,
            )
            m = min(len(ids), n)
            return [str(ids[i]) for i in range(m)], values[:m]

        return [str(s) for s in ids], values

    def _detect_zscore(
        self,
        locus_id: str,
        sample_ids: list[str],
        arr: np.ndarray,
        threshold: float,
    ) -> list[dict[str, Any]]:
        mean = float(np.mean(arr))
        std = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
        rows: list[dict[str, Any]] = []
        if std == 0.0 or not np.isfinite(std):
            logger.warning(
                "Locus %r (zscore): desvio-padrão nulo ou inválido; scores a 0",
                locus_id,
            )
            for sid, x in zip(sample_ids, arr.tolist(), strict=True):
                rows.append(
                    {
                        "sample_id": sid,
                        "locus_id": locus_id,
                        "allele_value": float(x),
                        "score": 0.0,
                        "is_outlier": False,
                    }
                )
            return rows

        for sid, x in zip(sample_ids, arr.tolist(), strict=True):
            z = (float(x) - mean) / std
            rows.append(
                {
                    "sample_id": sid,
                    "locus_id": locus_id,
                    "allele_value": float(x),
                    "score": float(z),
                    "is_outlier": abs(z) > threshold,
                }
            )
        return rows

    def _detect_iqr(
        self,
        locus_id: str,
        sample_ids: list[str],
        arr: np.ndarray,
        threshold: float,
    ) -> list[dict[str, Any]]:
        q1 = float(np.percentile(arr, 25.0))
        q3 = float(np.percentile(arr, 75.0))
        iqr = float(q3 - q1)
        lower = q1 - threshold * iqr
        upper = q3 + threshold * iqr
        rows: list[dict[str, Any]] = []

        if iqr == 0.0 or not np.isfinite(iqr):
            logger.warning(
                "Locus %r (iqr): IQR nulo ou inválido; nenhum outlier por IQR",
                locus_id,
            )
            for sid, x in zip(sample_ids, arr.tolist(), strict=True):
                rows.append(
                    {
                        "sample_id": sid,
                        "locus_id": locus_id,
                        "allele_value": float(x),
                        "score": 0.0,
                        "is_outlier": False,
                    }
                )
            return rows

        for sid, x in zip(sample_ids, arr.tolist(), strict=True):
            xf = float(x)
            if xf < lower:
                score = (lower - xf) / iqr
                is_out = True
            elif xf > upper:
                score = (xf - upper) / iqr
                is_out = True
            else:
                score = 0.0
                is_out = False
            rows.append(
                {
                    "sample_id": sid,
                    "locus_id": locus_id,
                    "allele_value": xf,
                    "score": float(score),
                    "is_outlier": is_out,
                }
            )
        return rows

    def _detect_mad(
        self,
        locus_id: str,
        sample_ids: list[str],
        arr: np.ndarray,
        threshold: float,
    ) -> list[dict[str, Any]]:
        med = float(np.median(arr))
        mad = float(
            stats.median_abs_deviation(arr, scale="normal", nan_policy="omit")
        )
        rows: list[dict[str, Any]] = []
        const = 0.6745

        if mad == 0.0 or not np.isfinite(mad):
            logger.warning(
                "Locus %r (mad): MAD nulo ou inválido; scores a 0",
                locus_id,
            )
            for sid, x in zip(sample_ids, arr.tolist(), strict=True):
                rows.append(
                    {
                        "sample_id": sid,
                        "locus_id": locus_id,
                        "allele_value": float(x),
                        "score": 0.0,
                        "is_outlier": False,
                    }
                )
            return rows

        for sid, x in zip(sample_ids, arr.tolist(), strict=True):
            xf = float(x)
            z_mad = const * (xf - med) / mad
            rows.append(
                {
                    "sample_id": sid,
                    "locus_id": locus_id,
                    "allele_value": xf,
                    "score": float(z_mad),
                    "is_outlier": abs(z_mad) > threshold,
                }
            )
        return rows
