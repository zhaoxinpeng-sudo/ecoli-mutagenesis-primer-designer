from __future__ import annotations

import re


CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def normalize_cds(raw: str) -> str:
    lines = [line.strip() for line in (raw or "").splitlines()]
    if lines and lines[0].startswith(">"):
        lines = [line for line in lines if not line.startswith(">")]
    seq = re.sub(r"\s+", "", "".join(lines)).upper().replace("U", "T")
    if not seq:
        raise ValueError("CDS 不能为空")
    bad = sorted(set(seq) - set("ACGT"))
    if bad:
        raise ValueError(f"CDS 含有非法碱基：{', '.join(bad)}")
    if len(seq) % 3:
        raise ValueError(f"CDS 长度为 {len(seq)} bp，不是 3 的倍数")
    return seq


def parse_mutation(text: str) -> tuple[str, int, str]:
    cleaned = re.sub(r"\s+", "", str(text or "")).upper()
    match = re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])(\d+)(?:->|-|>)*([ACDEFGHIKLMNPQRSTVWY])", cleaned)
    if not match:
        raise ValueError("突变格式无效，请使用 A123V 或 A123->V")
    old, position, new = match.group(1), int(match.group(2)), match.group(3)
    if position < 1:
        raise ValueError("氨基酸位置必须从 1 开始")
    return old, position, new


def translate_codon(codon: str) -> str:
    return CODON_TABLE[codon]


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

