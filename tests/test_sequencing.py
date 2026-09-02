import pandas as pd

from primer_designer.design import design_primers
from primer_designer.sequence import reverse_complement
from primer_designer.sequencing import RawRead, extract_sample_id, parse_seq, verify_sample
from primer_designer.verification_excel import read_verification_map, verification_results_bytes, verification_template_bytes

CDS = "ATG" + "GCT" * 18 + "GAA" + "GCT" * 25 + "TAA"

def mutated_sequence():
    design = design_primers(CDS, "E20V"); start = 19 * 3
    return CDS[:start] + design.optimized_codon + CDS[start + 3:], design.optimized_codon

def test_sample_id_and_seq_cleaning():
    assert extract_sample_id("1_T7_TSS.seq") == "1"
    assert extract_sample_id("ABC-T7TER.ab1") == "ABC"
    read = parse_seq(b">sample\nATG GCT\n123\n", "1_T7.seq")
    assert read.sequence == "ATGGCT" and read.sample_id == "1"

def test_exact_expected_codon_passes_and_reverse_is_detected():
    sequence, codon = mutated_sequence()
    result = verify_sample(CDS, "1", "E20V", [RawRead("1_T7.ab1", "1", sequence, [35] * len(sequence), source="AB1")])
    assert result.verdict == "通过" and result.expected_codon == result.actual_codon == codon
    reverse = RawRead("1_T7TER.ab1", "1", reverse_complement(sequence), [35] * len(sequence), source="AB1")
    result = verify_sample(CDS, "1", "E20V", [reverse])
    assert result.verdict == "通过" and result.directions == "反向"

def test_synonymous_target_codon_fails_exact_verification():
    sequence, expected = mutated_sequence(); start = 19 * 3
    alternative = next(x for x in ("GTG", "GTT", "GTC", "GTA") if x != expected)
    sequence = sequence[:start] + alternative + sequence[start + 3:]
    result = verify_sample(CDS, "1", "E20V", [RawRead("1.ab1", "1", sequence, [35]*len(sequence), source="AB1")])
    assert result.verdict == "失败" and "预期密码子" in result.reason

def test_extra_variant_fails_and_seq_only_reviews():
    sequence, _ = mutated_sequence(); changed = sequence[:30] + ("C" if sequence[30] != "C" else "A") + sequence[31:]
    failed = verify_sample(CDS, "1", "E20V", [RawRead("1.ab1", "1", changed, [35]*len(changed), source="AB1")])
    assert failed.verdict == "失败" and failed.extra_variants
    review = verify_sample(CDS, "1", "E20V", [RawRead("1.seq", "1", sequence, source="SEQ")])
    assert review.verdict == "需复核" and "无质量值" in review.reason

def test_chinese_mapping_and_report(tmp_path):
    template = tmp_path / "map.xlsx"; template.write_bytes(verification_template_bytes())
    assert list(read_verification_map(template).columns) == ["样品编号", "目标突变"]
    sequence, _ = mutated_sequence()
    result = verify_sample(CDS, "1", "E20V", [RawRead("1.ab1", "1", sequence, [35]*len(sequence), source="AB1")])
    output = tmp_path / "report.xlsx"; output.write_bytes(verification_results_bytes([result]))
    report = pd.read_excel(output)
    assert "预期密码子" in report.columns and report.iloc[0]["判定"] == "通过"
