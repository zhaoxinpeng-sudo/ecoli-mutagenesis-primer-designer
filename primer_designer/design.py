from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache

from .codons import ECOLI_CODONS
from .sequence import normalize_cds, parse_mutation, reverse_complement, translate_codon

RESTRICTION_SITES = ("GAATTC", "GGATCC", "AAGCTT", "GGTACC", "CTCGAG", "GTCGAC", "GCGGCCGC")


@dataclass(frozen=True)
class DesignParameters:
    min_flank: int = 20
    min_overlap: int = 20
    max_flank: int = 30
    target_tm: float = 68.0
    min_gc: float = 40.0
    max_gc: float = 60.0
    max_tm_difference: float = 5.0
    max_hairpin_run: int = 7
    max_dimer_run: int = 7


@dataclass
class DesignResult:
    sample_id: str
    mutation: str
    status: str
    message: str = ""
    original_codon: str = ""
    optimized_codon: str = ""
    forward_name: str = ""
    reverse_name: str = ""
    forward_primer: str = ""
    reverse_primer: str = ""
    forward_length: int | None = None
    reverse_length: int | None = None
    forward_tm: float | None = None
    reverse_tm: float | None = None
    forward_gc: float | None = None
    reverse_gc: float | None = None
    overlap_length: int | None = None
    left_flank: int | None = None
    right_flank: int | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def gc_percent(seq: str) -> float:
    return 100.0 * sum(base in "GC" for base in seq) / len(seq)


def melting_temperature(seq: str) -> float:
    if len(seq) < 14:
        return float(2 * sum(base in "AT" for base in seq) + 4 * sum(base in "GC" for base in seq))
    return 64.9 + 41.0 * (sum(base in "GC" for base in seq) - 16.4) / len(seq)


@lru_cache(maxsize=32768)
def longest_complementary_run(a: str, b: str) -> int:
    target = reverse_complement(b)
    best = 0
    for offset in range(-len(target) + 1, len(a)):
        run = 0
        for i, base in enumerate(a):
            j = i - offset
            if 0 <= j < len(target) and base == target[j]:
                run += 1
                best = max(best, run)
            else:
                run = 0
    return best


@lru_cache(maxsize=8192)
def hairpin_run(seq: str) -> int:
    best = 0
    for split in range(4, len(seq) - 7):
        left, right = seq[:split], seq[split + 3:]
        best = max(best, longest_complementary_run(left, right))
    return best


@lru_cache(maxsize=32768)
def three_prime_dimer_run(a: str, b: str) -> int:
    """Return complementarity runs that touch a 3' end.

    The intentional mutagenic overlap is located at the 5' ends of both
    primers and must not be treated as an unwanted primer dimer.
    """
    target = reverse_complement(b)
    best = 0
    for offset in range(-len(target) + 1, len(a)):
        run_start_a = None
        run_start_target = None
        run = 0
        for i, base in enumerate(a):
            j = i - offset
            if 0 <= j < len(target) and base == target[j]:
                if run == 0:
                    run_start_a, run_start_target = i, j
                run += 1
                # In this orientation, a's 3' end is its right edge and
                # b's 3' end is the left edge of reverse_complement(b).
                if i == len(a) - 1 or run_start_target == 0:
                    best = max(best, run)
            else:
                run = 0
                run_start_a = run_start_target = None
    return best


def newly_created_sites(original: str, mutated: str) -> list[str]:
    return [site for site in RESTRICTION_SITES if site in mutated and site not in original]


def _candidate_score(forward: str, reverse: str, codon_rank: int, params: DesignParameters) -> tuple[float, dict]:
    ftm, rtm = melting_temperature(forward), melting_temperature(reverse)
    fgc, rgc = gc_percent(forward), gc_percent(reverse)
    hp = max(hairpin_run(forward), hairpin_run(reverse))
    # Cross-primer complementarity is intentional for overlapping mutagenic
    # primers. Screen only dangerous 3'-anchored self-dimers here.
    dimer = max(three_prime_dimer_run(forward, forward), three_prime_dimer_run(reverse, reverse))
    gc_penalty = max(0, params.min_gc - fgc) + max(0, fgc - params.max_gc)
    gc_penalty += max(0, params.min_gc - rgc) + max(0, rgc - params.max_gc)
    score = abs(ftm - params.target_tm) + abs(rtm - params.target_tm) + 2 * abs(ftm - rtm)
    score += 2 * gc_penalty + 3 * max(0, hp - 4) + 3 * max(0, dimer - 4) + codon_rank * 1.5
    return score, {"ftm": ftm, "rtm": rtm, "fgc": fgc, "rgc": rgc, "hp": hp, "dimer": dimer}


def design_primers(cds: str, mutation: str, sample_id: str = "", params: DesignParameters | None = None) -> DesignResult:
    params = params or DesignParameters()
    result = DesignResult(sample_id=str(sample_id), mutation=str(mutation), status="失败")
    try:
        seq = normalize_cds(cds)
        old_aa, position, new_aa = parse_mutation(mutation)
        canonical_mutation = f"{old_aa}{position}{new_aa}"
        result.mutation = canonical_mutation
        result.forward_name = f"{canonical_mutation}-F"
        result.reverse_name = f"{canonical_mutation}-R"
        if position > len(seq) // 3:
            raise ValueError(f"位点 {position} 超出 CDS 编码范围（共 {len(seq)//3} 个密码子）")
        codon_start = (position - 1) * 3
        original_codon = seq[codon_start:codon_start + 3]
        actual = translate_codon(original_codon)
        result.original_codon = original_codon
        if actual != old_aa:
            raise ValueError(f"原氨基酸不匹配：位点 {position} 实际为 {actual}（{original_codon}），不是 {old_aa}")
        if codon_start < params.min_flank or len(seq) - (codon_start + 3) < params.min_flank:
            raise ValueError("突变位点距离 CDS 边界不足，无法在左右两侧各保留至少 20 bp")

        best: tuple[float, str, str, str, int, int, int, dict] | None = None
        original_window = seq[max(0, codon_start - 10):min(len(seq), codon_start + 13)]
        for rank, codon in enumerate(ECOLI_CODONS[new_aa]):
            mutated = seq[:codon_start] + codon + seq[codon_start + 3:]
            mutated_window = mutated[max(0, codon_start - 10):min(len(mutated), codon_start + 13)]
            if newly_created_sites(original_window, mutated_window) or any(base * 5 in mutated_window for base in "ACGT"):
                continue
            # The overlap contains the complete mutated codon and is shared by both primers.
            for overlap_left in range(0, params.min_overlap - 2):
                overlap_right = params.min_overlap - 3 - overlap_left
                overlap_start = codon_start - overlap_left
                overlap_end = codon_start + 3 + overlap_right
                if overlap_start < 0 or overlap_end > len(mutated):
                    continue
                for left_flank in range(params.min_flank, params.max_flank + 1):
                    left_start = codon_start - left_flank
                    if left_start < 0:
                        continue
                    reverse_template = mutated[left_start:overlap_end]
                    reverse = reverse_complement(reverse_template)
                    for right_flank in range(params.min_flank, params.max_flank + 1):
                        right_end = codon_start + 3 + right_flank
                        if right_end > len(mutated):
                            continue
                        forward = mutated[overlap_start:right_end]
                        score, metrics = _candidate_score(forward, reverse, rank, params)
                        if metrics["hp"] > params.max_hairpin_run or metrics["dimer"] > params.max_dimer_run:
                            continue
                        if abs(metrics["ftm"] - metrics["rtm"]) > params.max_tm_difference:
                            continue
                        candidate = (score, codon, forward, reverse, left_flank, right_flank, overlap_end-overlap_start, metrics)
                        if best is None or candidate[0] < best[0]:
                            best = candidate
        if best is None:
            raise ValueError("没有候选方案同时通过同源长度、Tm、GC及二级结构筛选")
        _, codon, forward, reverse, left_flank, right_flank, overlap, metrics = best
        result.status = "成功"
        result.message = "设计成功"
        result.optimized_codon = codon
        result.forward_primer = forward
        result.reverse_primer = reverse
        result.forward_length = len(forward)
        result.reverse_length = len(reverse)
        result.forward_tm = round(metrics["ftm"], 1)
        result.reverse_tm = round(metrics["rtm"], 1)
        result.forward_gc = round(metrics["fgc"], 1)
        result.reverse_gc = round(metrics["rgc"], 1)
        result.overlap_length = overlap
        result.left_flank = left_flank
        result.right_flank = right_flank
        return result
    except (ValueError, KeyError) as exc:
        result.message = str(exc)
        return result
