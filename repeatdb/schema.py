"""Modelos de dados (Pydantic) para loci, chamadas e amostras."""

from pydantic import BaseModel


class Locus(BaseModel):
    locus_id: str
    gene: str
    motif: str
    chrom: str
    start: int
    end: int


class RepeatCall(BaseModel):
    sample_id: str
    locus_id: str
    allele1: int | None
    allele2: int | None
    filter: str
    genotype: str


class Sample(BaseModel):
    sample_id: str
    vcf_path: str
    metadata: dict
