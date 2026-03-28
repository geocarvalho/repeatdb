"""Leitura e interpretação de VCFs do ExpansionHunter (cyvcf2)."""

from __future__ import annotations

import logging
from typing import Any

import cyvcf2

from repeatdb.schema import RepeatCall

logger = logging.getLogger(__name__)


def _filter_to_string(vcf_filter: str | None) -> str:
    if vcf_filter is None:
        return "PASS"
    return str(vcf_filter)


def _parse_repcn(raw: Any) -> tuple[int | None, int | None]:
    """Parse REPCN: ``X/Y``, ``X``, ``.``, or empty → (allele1, allele2)."""
    if raw is None:
        return None, None
    s = str(raw).strip()
    if not s or s == ".":
        return None, None
    if "/" in s:
        parts = s.split("/")
        if len(parts) != 2:
            raise ValueError(f"REPCN diplóide inválido: {s!r}")

        def one(x: str) -> int | None:
            t = x.strip()
            if not t or t == ".":
                return None
            return int(t)

        return one(parts[0]), one(parts[1])
    return int(s), None


def _genotype_tuple_to_string(gt: list[Any]) -> str:
    """Converte o tuple de genótipo do cyvcf2 em string estilo VCF (ex.: ``0/1``, ``.|.``)."""
    if not gt:
        return "."

    def al(x: Any) -> str:
        try:
            xi = int(x)
        except (TypeError, ValueError):
            return "."
        return "." if xi < 0 else str(xi)

    if len(gt) == 2:
        return al(gt[0])

    if len(gt) >= 3:
        a1, a2, phased = gt[0], gt[1], gt[2]
        sep = "|" if phased else "/"
        return f"{al(a1)}{sep}{al(a2)}"

    return "."


def _variant_to_call(
    rec: cyvcf2.Variant,
    sample_id: str,
    sample_idx: int,
) -> RepeatCall | None:
    repid = rec.INFO.get("REPID")
    if repid is None or repid == "":
        logger.warning(
            "Variante %s:%s sem INFO/REPID; registo ignorado",
            rec.CHROM,
            rec.POS,
        )
        return None
    locus_id = str(repid)

    try:
        repcn_arr = rec.format("REPCN")
    except KeyError:
        logger.warning(
            "Variante %s:%s sem FORMAT/REPCN; usando no-call",
            rec.CHROM,
            rec.POS,
        )
        allele1, allele2 = None, None
    else:
        raw_repcn = repcn_arr[sample_idx]
        if raw_repcn is None:
            s = ""
        else:
            s = str(raw_repcn).strip()
        try:
            allele1, allele2 = _parse_repcn(s if s else None)
        except ValueError as e:
            logger.warning(
                "Variante %s:%s REPCN %r inválido: %s",
                rec.CHROM,
                rec.POS,
                s,
                e,
            )
            allele1, allele2 = None, None

    filt = _filter_to_string(rec.FILTER)

    gts = rec.genotypes
    if sample_idx >= len(gts):
        logger.warning(
            "Variante %s:%s sem genótipo para amostra índice %s",
            rec.CHROM,
            rec.POS,
            sample_idx,
        )
        genotype = "."
    else:
        genotype = _genotype_tuple_to_string(list(gts[sample_idx]))

    return RepeatCall(
        sample_id=sample_id,
        locus_id=locus_id,
        allele1=allele1,
        allele2=allele2,
        filter=filt,
        genotype=genotype,
    )


def parse_vcf(vcf_path: str) -> list[RepeatCall]:
    """
    Lê um VCF do ExpansionHunter e devolve uma lista de :class:`RepeatCall`.

    Assume uma única amostra por ficheiro (padrão ExpansionHunter). Se existirem
    várias, usa a primeira e regista um aviso.
    """
    calls: list[RepeatCall] = []

    try:
        vcf = cyvcf2.VCF(vcf_path)
    except Exception as e:
        logger.warning("Não foi possível abrir o VCF %s: %s", vcf_path, e)
        return calls

    with vcf:
        samples = vcf.samples
        if not samples:
            logger.warning("VCF %s não tem colunas de amostra", vcf_path)
            return calls

        if len(samples) > 1:
            logger.warning(
                "VCF %s tem %d amostras; a usar apenas a primeira (%r)",
                vcf_path,
                len(samples),
                samples[0],
            )

        sample_id = str(samples[0])
        sample_idx = 0

        try:
            for rec in vcf:
                try:
                    call = _variant_to_call(rec, sample_id, sample_idx)
                    if call is not None:
                        calls.append(call)
                except Exception as e:
                    logger.warning(
                        "Variante %s:%s ignorada (erro ao processar): %s",
                        getattr(rec, "CHROM", "?"),
                        getattr(rec, "POS", "?"),
                        e,
                    )
        except Exception as e:
            logger.warning("Erro ao iterar o VCF %s: %s", vcf_path, e)

    return calls
