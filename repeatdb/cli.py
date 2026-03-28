"""Interface de linha de comando (Typer)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import polars as pl
import typer
from rich.console import Console
from rich.table import Table

from repeatdb.aggregator import build_locus_matrix_with_sample_ids
from repeatdb.outlier import OutlierDetector
from repeatdb.reporter import generate_report
from repeatdb.store import ingest_vcf_list, load_parquet, save_parquet

app = typer.Typer(no_args_is_help=True, help="repeatdb — repeat expansions (ExpansionHunter)")
console = Console()


def _read_vcf_paths(list_file: Path) -> list[str]:
    text = list_file.read_text(encoding="utf-8")
    return [line.strip() for line in text.splitlines() if line.strip()]


@app.command("ingest")
def cmd_ingest(
    vcf_list: Path = typer.Option(
        ...,
        "--vcf-list",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Ficheiro com um path de VCF por linha",
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        help="Path do ficheiro Parquet de saída",
    ),
    workers: int = typer.Option(4, "--workers", min=1, help="Processos paralelos"),
) -> None:
    paths = _read_vcf_paths(vcf_list)
    if not paths:
        console.print("[yellow]Lista de VCFs vazia; nada a fazer.[/yellow]")
        raise typer.Exit(code=1)
    df = ingest_vcf_list(paths, n_workers=workers)
    save_parquet(df, str(output))
    console.print(f"[green]Ingestão concluída:[/green] {df.height} linhas → {output}")


@app.command("detect")
def cmd_detect(
    input_parquet: Path = typer.Option(
        ...,
        "--input",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Parquet de chamadas (formato store)",
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        help="Path do relatório de outliers",
    ),
    method: str = typer.Option(
        "mad",
        "--method",
        help="mad, iqr ou zscore",
    ),
    threshold: float = typer.Option(3.5, "--threshold", help="Limiar do método"),
    out_format: str = typer.Option(
        "tsv",
        "--format",
        help="tsv, json ou html",
    ),
) -> None:
    m = method.strip().lower()
    if m not in {"mad", "iqr", "zscore"}:
        console.print(f"[red]Método inválido:[/red] {method!r}")
        raise typer.Exit(code=1)

    fmt = out_format.strip().lower()
    if fmt not in {"tsv", "json", "html"}:
        console.print(f"[red]Formato inválido:[/red] {out_format!r}")
        raise typer.Exit(code=1)

    df = load_parquet(str(input_parquet))
    matrix, sample_ids = build_locus_matrix_with_sample_ids(df)
    detector = OutlierDetector(sample_ids_by_locus=sample_ids)
    rows = detector.detect(matrix, method=m, threshold=threshold)

    if rows:
        outliers_df = pl.DataFrame(rows)
    else:
        outliers_df = pl.DataFrame(
            schema={
                "sample_id": pl.Utf8,
                "locus_id": pl.Utf8,
                "allele_value": pl.Float64,
                "score": pl.Float64,
                "is_outlier": pl.Boolean,
            }
        )

    generate_report(outliers_df, str(output), fmt)
    n_out = int(outliers_df.filter(pl.col("is_outlier")).height) if outliers_df.height else 0
    console.print(
        f"[green]Deteção concluída:[/green] {outliers_df.height} linhas "
        f"({n_out} outliers) → {output}"
    )


@app.command("query")
def cmd_query(
    input_parquet: Path = typer.Option(
        ...,
        "--input",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Parquet de chamadas",
    ),
    gene: str = typer.Option(
        ...,
        "--gene",
        help="Identificador do locus (coluna locus_id / REPID)",
    ),
    sample: Optional[str] = typer.Option(
        None,
        "--sample",
        help="Filtrar por sample_id (opcional)",
    ),
) -> None:
    df = load_parquet(str(input_parquet))
    q = df.filter(pl.col("locus_id") == gene)
    if sample is not None:
        q = q.filter(pl.col("sample_id") == sample)

    if q.is_empty():
        console.print("[yellow]Sem linhas para os filtros indicados.[/yellow]")
        raise typer.Exit(code=0)

    table = Table(show_header=True, header_style="bold")
    for col in q.columns:
        table.add_column(col)

    for row in q.iter_rows():
        table.add_row(*(str(x) if x is not None else "" for x in row))

    console.print(table)


def main() -> None:
    app()


if __name__ == "__main__":
    main()
