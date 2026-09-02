from __future__ import annotations

from dataclasses import dataclass, field
from io import BytesIO
from pathlib import PurePath
import re
from zipfile import ZipFile, BadZipFile

from Bio import Align, SeqIO

from .design import DesignParameters, design_primers
from .sequence import normalize_cds, parse_mutation, reverse_complement, translate_codon


@dataclass
class RawRead:
    filename: str
    sample_id: str
    sequence: str
    qualities: list[int] | None = None
    peaks: list[int] = field(default_factory=list)
    traces: dict[str, list[int]] = field(default_factory=dict)
    source: str = "SEQ"


@dataclass
class AlignedRead:
    raw: RawRead
    sequence: str
    qualities: list[int] | None
    direction: str
    score: float
    ref_to_read: dict[int, int]
    insertions: list[str]
    deletions: list[str]
    coverage_start: int | None
    coverage_end: int | None


@dataclass
class VerificationResult:
    sample_id: str
    mutation: str
    expected_codon: str = ""
    actual_codon: str = ""
    verdict: str = "需复核"
    reason: str = ""
    valid_reads: int = 0
    directions: str = ""
    target_min_quality: int | None = None
    coverage_percent: float = 0.0
    coverage_range: str = "—"
    extra_variants: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    consensus: str = ""
    aligned_reads: list[AlignedRead] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "样品编号": self.sample_id, "目标突变": self.mutation,
            "预期密码子": self.expected_codon, "实际密码子": self.actual_codon,
            "判定": self.verdict, "有效读段数": self.valid_reads,
            "测序方向": self.directions, "参考覆盖率（%）": self.coverage_percent,
            "目标最低Q值": self.target_min_quality, "覆盖范围": self.coverage_range,
            "额外变异": "；".join(self.extra_variants),
            "警告/原因": "；".join(filter(None, [self.reason, *self.warnings])),
        }


def extract_sample_id(filename: str) -> str:
    stem = PurePath(filename).name.rsplit(".", 1)[0]
    sample = re.split(r"[_-]", stem, maxsplit=1)[0].strip()
    return sample or stem


def parse_seq(data: bytes, filename: str) -> RawRead:
    text = None
    for encoding in ("utf-8-sig", "gb18030", "latin1"):
        try:
            text = data.decode(encoding)
            break
        except UnicodeDecodeError:
            pass
    lines = [line for line in (text or "").splitlines() if not line.lstrip().startswith(">")]
    sequence = re.sub(r"[^A-Za-z]", "", "".join(lines)).upper().replace("U", "T")
    bad = sorted(set(sequence) - set("ACGTN"))
    if not sequence or bad:
        raise ValueError(f"{filename} 不是有效的 DNA 序列文件")
    return RawRead(filename, extract_sample_id(filename), sequence, source="SEQ")


def parse_ab1(data: bytes, filename: str) -> RawRead:
    record = SeqIO.read(BytesIO(data), "abi")
    sequence = str(record.seq).upper()
    qualities = list(record.letter_annotations.get("phred_quality", [])) or None
    raw = record.annotations.get("abif_raw", {})
    peaks = list(raw.get("PLOC2", raw.get("PLOC1", [])))
    order = raw.get("FWO_1", b"GATC")
    if isinstance(order, bytes):
        order = order.decode("ascii", errors="ignore")
    traces = {}
    for index, base in enumerate(str(order)[:4]):
        values = raw.get(f"DATA{9 + index}")
        if values is not None:
            traces[base] = list(values)
    return RawRead(filename, extract_sample_id(filename), sequence, qualities, peaks, traces, "AB1")


def unpack_files(files) -> tuple[list[RawRead], list[str]]:
    payloads: list[tuple[str, bytes]] = []
    warnings: list[str] = []
    for uploaded in files or []:
        name = PurePath(uploaded.name).name
        data = uploaded.getvalue()
        if name.lower().endswith(".zip"):
            try:
                with ZipFile(BytesIO(data)) as archive:
                    for member in archive.infolist():
                        leaf = PurePath(member.filename).name
                        if not member.is_dir() and leaf.lower().endswith((".ab1", ".seq")):
                            payloads.append((leaf, archive.read(member)))
            except BadZipFile:
                warnings.append(f"{name} 不是有效的 ZIP 文件")
        elif name.lower().endswith((".ab1", ".seq")):
            payloads.append((name, data))

    # 同名文件以 AB1 为主。
    chosen: dict[str, tuple[str, bytes]] = {}
    for name, data in payloads:
        key = PurePath(name).stem.lower()
        if key not in chosen or name.lower().endswith(".ab1"):
            chosen[key] = (name, data)
    reads = []
    for name, data in chosen.values():
        try:
            reads.append(parse_ab1(data, name) if name.lower().endswith(".ab1") else parse_seq(data, name))
        except Exception as exc:
            warnings.append(f"{name}：{exc}")
    return reads, warnings


def q20_trim(read: RawRead, threshold: int = 20, min_length: int = 50, window: int = 10) -> RawRead | None:
    if read.qualities is None:
        return read if len(read.sequence) >= min_length else None
    qualities = read.qualities[:len(read.sequence)]
    left, right = 0, len(qualities)
    while left + window <= right and sum(qualities[left:left + window]) / window < threshold:
        left += 1
    while right - window >= left and sum(qualities[right - window:right]) / window < threshold:
        right -= 1
    if right - left < min_length:
        return None
    return RawRead(read.filename, read.sample_id, read.sequence[left:right], qualities[left:right],
                   read.peaks[left:right] if read.peaks else [], read.traces, read.source)


def _align_one(reference: str, read: RawRead) -> AlignedRead:
    aligner = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-2,
                                    open_gap_score=-5, extend_gap_score=-1)
    candidates = []
    for direction, sequence, qualities in (
        ("正向", read.sequence, read.qualities),
        ("反向", reverse_complement(read.sequence.replace("N", "A")), list(reversed(read.qualities)) if read.qualities else None),
    ):
        alignment = aligner.align(reference, sequence)[0]
        candidates.append((alignment.score, direction, sequence, qualities, alignment))
    score, direction, sequence, qualities, alignment = max(candidates, key=lambda item: item[0])
    coords = alignment.coordinates
    mapping: dict[int, int] = {}
    insertions: list[str] = []
    deletions: list[str] = []
    for i in range(coords.shape[1] - 1):
        r0, r1 = int(coords[0, i]), int(coords[0, i + 1])
        q0, q1 = int(coords[1, i]), int(coords[1, i + 1])
        if r1 > r0 and q1 > q0:
            for offset in range(min(r1 - r0, q1 - q0)):
                mapping[r0 + offset] = q0 + offset
        elif r1 == r0 and q1 > q0:
            insertions.append(f"c.{r0}_{r0 + 1}ins{sequence[q0:q1]}")
        elif q1 == q0 and r1 > r0:
            deletions.append(f"c.{r0 + 1}_{r1}del")
    covered = sorted(mapping)
    return AlignedRead(read, sequence, qualities, direction, score, mapping, insertions, deletions,
                       covered[0] if covered else None, covered[-1] if covered else None)


def verify_sample(reference: str, sample_id: str, mutation: str, reads: list[RawRead],
                  q_threshold: int = 20, min_length: int = 50) -> VerificationResult:
    result = VerificationResult(str(sample_id), str(mutation))
    try:
        reference = normalize_cds(reference)
        old_aa, position, new_aa = parse_mutation(mutation)
        codon_start = (position - 1) * 3
        if codon_start + 3 > len(reference):
            raise ValueError("目标位点超出参考 CDS")
        if translate_codon(reference[codon_start:codon_start + 3]) != old_aa:
            raise ValueError("参考 CDS 的原氨基酸与目标突变不一致")
        design = design_primers(reference, mutation, params=DesignParameters())
        if design.status != "成功":
            raise ValueError(f"无法按引物设计规则得到预期密码子：{design.message}")
        result.mutation = design.mutation
        result.expected_codon = design.optimized_codon
        trimmed = [item for read in reads if (item := q20_trim(read, q_threshold, min_length))]
        result.valid_reads = len(trimmed)
        if not trimmed:
            result.reason = "没有达到长度和质量要求的有效读段"
            return result
        aligned = [_align_one(reference, read) for read in trimmed]
        result.aligned_reads = aligned
        result.directions = "、".join(sorted({item.direction for item in aligned}))

        calls: dict[int, list[tuple[str, int | None, str]]] = {}
        for item in aligned:
            for ref_pos, read_pos in item.ref_to_read.items():
                base = item.sequence[read_pos]
                quality = item.qualities[read_pos] if item.qualities else None
                calls.setdefault(ref_pos, []).append((base, quality, item.raw.source))
        consensus = list("N" * len(reference))
        consensus_q: dict[int, int | None] = {}
        conflicts = []
        for pos, values in calls.items():
            ranked = sorted(values, key=lambda x: -1 if x[1] is None else x[1], reverse=True)
            best = ranked[0]
            disagree = [x for x in ranked[1:] if x[0] != best[0]]
            if disagree and best[1] is not None and disagree[0][1] is not None and abs(best[1] - disagree[0][1]) < 5:
                consensus[pos] = "N"
                conflicts.append(pos + 1)
            else:
                consensus[pos] = best[0]
                consensus_q[pos] = best[1]
        result.consensus = "".join(consensus)
        covered = sorted(calls)
        result.coverage_percent = round(100 * len(covered) / len(reference), 1)
        if covered:
            result.coverage_range = f"c.{covered[0] + 1}–{covered[-1] + 1}"
        result.actual_codon = result.consensus[codon_start:codon_start + 3]
        target_q = [consensus_q.get(i) for i in range(codon_start, codon_start + 3)]
        numeric_q = [q for q in target_q if q is not None]
        result.target_min_quality = min(numeric_q) if len(numeric_q) == 3 else None

        extras = []
        for pos in covered:
            if codon_start <= pos < codon_start + 3:
                continue
            base = result.consensus[pos]
            quality = consensus_q.get(pos)
            if base in "ACGT" and base != reference[pos] and quality is not None and quality >= q_threshold:
                extras.append(f"c.{pos + 1}{reference[pos]}>{base}")
        for item in aligned:
            extras.extend(item.insertions + item.deletions)
        result.extra_variants = sorted(set(extras))
        if conflicts:
            result.warnings.append("读段冲突：" + "、".join(f"c.{x}" for x in conflicts[:12]))
        seq_only = all(read.source == "SEQ" for read in trimmed)
        target_complete = all(base in "ACGT" for base in result.actual_codon)
        target_high_quality = result.target_min_quality is not None and result.target_min_quality >= q_threshold

        if target_complete and target_high_quality and result.actual_codon != result.expected_codon:
            result.verdict, result.reason = "失败", "目标密码子与引物设计预期密码子不一致"
        elif result.extra_variants:
            result.verdict, result.reason = "失败", "检测到目标密码子之外的高质量变异"
        elif result.actual_codon == result.expected_codon and target_high_quality and not conflicts and not seq_only:
            result.verdict, result.reason = "通过", "目标密码子正确，未检出额外高质量变异"
        else:
            result.verdict = "需复核"
            if not target_complete:
                result.reason = "目标密码子未完整覆盖或存在读段冲突"
            elif seq_only:
                result.reason = "仅有无质量值的 SEQ 文件支持"
            elif not target_high_quality:
                result.reason = "目标密码子质量不足"
            else:
                result.reason = "目标密码子需要人工复核"
        return result
    except (ValueError, KeyError) as exc:
        result.reason = str(exc)
        return result


def local_alignment(result: VerificationResult, flank: int = 24) -> tuple[str, str, str]:
    try:
        _, position, _ = parse_mutation(result.mutation)
        start = max(0, (position - 1) * 3 - flank)
        end = min(len(result.consensus), (position - 1) * 3 + 3 + flank)
        ref = getattr(result, "_reference", "")
        con = result.consensus[start:end]
        return f"c.{start + 1}–{end}", ref[start:end] if ref else "", con
    except ValueError:
        return "", "", ""
