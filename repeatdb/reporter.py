"""Relatórios de outliers (TSV, JSON, HTML)."""

from __future__ import annotations

import html
import json
from pathlib import Path

import polars as pl

_HTML_COLUMNS = ("sample_id", "locus_id", "allele_value", "score", "is_outlier")


def generate_report(
    outliers_df: pl.DataFrame,
    output_path: str,
    format: str,
) -> None:
    """
    Grava ``outliers_df`` em ``output_path``.

    * ``tsv`` — todas as colunas, separador TAB.
    * ``json`` — todas as colunas, lista de objetos com ``indent=2``.
    * ``html`` — página autónoma com DataTables (CDN); tabela só com as colunas
      de outlier; linhas com ``is_outlier=True`` com fundo âmbar suave.
    """
    fmt = format.strip().lower()
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "tsv":
        outliers_df.write_csv(path, separator="\t")
    elif fmt == "json":
        records = outliers_df.to_dicts()
        path.write_text(json.dumps(records, indent=2), encoding="utf-8")
    elif fmt == "html":
        path.write_text(_build_html_datatables(outliers_df), encoding="utf-8")
    else:
        raise ValueError(
            f"formato inválido {format!r}; use 'tsv', 'json' ou 'html'"
        )


def _cell_str(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if v != v:  # NaN
            return ""
        return repr(v) if v.is_integer() else f"{v:.6g}"
    return str(v)


def _build_html_datatables(df: pl.DataFrame) -> str:
    missing = [c for c in _HTML_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"DataFrame em falta colunas para HTML: {missing}")

    view = df.select([pl.col(c) for c in _HTML_COLUMNS])
    rows_html: list[str] = []
    for row in view.iter_rows(named=True):
        is_out = bool(row["is_outlier"])
        tr_cls = ' class="outlier-row"' if is_out else ""
        tds = []
        for col in _HTML_COLUMNS:
            raw = row[col]
            if col == "is_outlier":
                text = "true" if raw else "false"
            else:
                text = _cell_str(raw)
            tds.append(f"<td>{html.escape(text)}</td>")
        rows_html.append(f"<tr{tr_cls}>{''.join(tds)}</tr>")

    tbody = "\n".join(rows_html) if rows_html else ""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>repeatdb — outliers</title>
  <link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css"/>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 1.5rem; }}
    h1 {{ font-size: 1.25rem; }}
    tr.outlier-row td {{ background-color: #fff3cd !important; }}
  </style>
</head>
<body>
  <h1>Outliers</h1>
  <table id="outliers" class="display compact nowrap" style="width:100%">
    <thead>
      <tr>
        <th>sample_id</th>
        <th>locus_id</th>
        <th>allele_value</th>
        <th>score</th>
        <th>is_outlier</th>
      </tr>
    </thead>
    <tbody>
{tbody}
    </tbody>
  </table>
  <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
  <script>
    $(function () {{
      $("#outliers").DataTable({{ pageLength: 25, order: [[3, "desc"]] }});
    }});
  </script>
</body>
</html>
"""
