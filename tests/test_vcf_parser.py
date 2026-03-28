"""Testes do parser de VCF (ExpansionHunter)."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import pytest

from repeatdb.schema import RepeatCall
from repeatdb.vcf_parser import parse_vcf

# VCF mínimo estilo ExpansionHunter (cyvcf2 exige ficheiro no disco)
_EH_VCF_HEADER = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=REPCN,Number=.,Type=String,Description="Repeat counts per allele">
##FORMAT=<ID=REPCI,Number=.,Type=String,Description="Confidence interval">
##FORMAT=<ID=LC,Number=1,Type=Float,Description="Local coverage">
##INFO=<ID=REPID,Number=1,Type=String,Description="Repeat ID">
##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference length in units">
##contig=<ID=chr4,length=191154276>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A
"""


def _write_vcf(content: str) -> Path:
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".vcf",
        delete=False,
        encoding="utf-8",
        newline="\n",
    ) as f:
        f.write(content)
        return Path(f.name)


def test_parse_vcf_diploid_repcn_and_pass() -> None:
    body = "chr4\t3076604\t.\tCAG\tCAG,C\t.\tPASS\tREPID=HTT;RU=CAG;RL=19\tGT:REPCN:REPCI:LC\t0/1:19/22:18-20/21-23:150.5\n"
    path = _write_vcf(_EH_VCF_HEADER + body)
    try:
        got = parse_vcf(str(path))
    finally:
        path.unlink(missing_ok=True)

    assert got == [
        RepeatCall(
            sample_id="SAMPLE_A",
            locus_id="HTT",
            allele1=19,
            allele2=22,
            filter="PASS",
            genotype="0/1",
        )
    ]


def test_parse_vcf_haploid_repcn_single_sample() -> None:
    header = _EH_VCF_HEADER.replace("SAMPLE_A", "MALE_X")
    body = "chrX\t123\t.\tC\tC\t.\tPASS\tREPID=AR\tGT:REPCN:LC\t0:42:90.0\n"
    path = _write_vcf(header + body)
    try:
        got = parse_vcf(str(path))
    finally:
        path.unlink(missing_ok=True)

    assert len(got) == 1
    assert got[0].allele1 == 42
    assert got[0].allele2 is None
    assert got[0].genotype == "0"


def test_parse_vcf_repcn_dot_no_call() -> None:
    body = "chr4\t1\t.\tC\tC\t.\tPASS\tREPID=FOO\tGT:REPCN\t./.:.\n"
    path = _write_vcf(_EH_VCF_HEADER + body)
    try:
        got = parse_vcf(str(path))
    finally:
        path.unlink(missing_ok=True)

    assert got[0].allele1 is None
    assert got[0].allele2 is None
    assert got[0].genotype == "./."


def test_parse_vcf_filter_lowqual() -> None:
    body = "chr4\t1\t.\tC\tC\t.\tLowQual\tREPID=Z\tGT:REPCN\t0/1:5/5\n"
    path = _write_vcf(_EH_VCF_HEADER + body)
    try:
        got = parse_vcf(str(path))
    finally:
        path.unlink(missing_ok=True)

    assert got[0].filter == "LowQual"


def test_parse_vcf_malformed_path_returns_empty(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.WARNING)
    got = parse_vcf("/nonexistent/repeatdb_fake.vcf")
    assert got == []
    assert any("Não foi possível abrir" in r.message for r in caplog.records)


def test_parse_vcf_skips_missing_repid(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.WARNING)
    body = "chr4\t1\t.\tC\tC\t.\tPASS\t.\tGT:REPCN\t0/1:1/1\n"
    path = _write_vcf(_EH_VCF_HEADER + body)
    try:
        got = parse_vcf(str(path))
    finally:
        path.unlink(missing_ok=True)

    assert got == []
    assert any("sem INFO/REPID" in r.message for r in caplog.records)
