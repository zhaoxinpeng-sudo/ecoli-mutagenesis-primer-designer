import pandas as pd
import pytest

from primer_designer.design import DesignParameters, design_primers
from primer_designer.excel_io import flatten_results, read_mutations, results_bytes, template_bytes
from primer_designer.sequence import normalize_cds, parse_mutation, reverse_complement, translate_codon


CDS = "ATG" + "GCT" * 18 + "GAA" + "GCT" * 25 + "TAA"


def test_normalize_fasta_and_validate():
    assert normalize_cds(">x\natg gct\n") == "ATGGCT"
    with pytest.raises(ValueError): normalize_cds("ATGN")
    with pytest.raises(ValueError): normalize_cds("ATGG")


@pytest.mark.parametrize("text,expected", [("A2V", ("A", 2, "V")), ("A2->V", ("A", 2, "V")), ("A2-V", ("A", 2, "V"))])
def test_mutation_formats(text, expected):
    assert parse_mutation(text) == expected


def test_success_and_hard_constraints():
    result = design_primers(CDS, "E20V", "s1")
    assert result.status == "成功", result.message
    assert result.left_flank >= 20 and result.right_flank >= 20 and result.overlap_length >= 20
    assert result.optimized_codon in ("GTG", "GTT", "GTC", "GTA")
    assert translate_codon(result.optimized_codon) == "V"
    assert result.optimized_codon in result.forward_primer
    assert reverse_complement(result.optimized_codon) in result.reverse_primer
    assert result.forward_name == "E20V-F"
    assert result.reverse_name == "E20V-R"


def test_primer_names_use_canonical_mutation():
    result = design_primers(CDS, "E20->V")
    assert result.mutation == "E20V"
    assert result.forward_name == "E20V-F"
    assert result.reverse_name == "E20V-R"


def test_success_result_is_two_primer_rows_without_codon_or_status():
    record = design_primers(CDS, "E20V", "x").to_dict()
    rows = flatten_results([record])
    assert [row["primer_name"] for row in rows] == ["E20V-F", "E20V-R"]
    assert [row["direction"] for row in rows] == ["上游引物", "下游引物"]
    assert all("status" not in row and "original_codon" not in row and "optimized_codon" not in row for row in rows)


def test_original_amino_acid_mismatch():
    result = design_primers(CDS, "A20V")
    assert result.status == "失败"
    assert "实际为 E" in result.message
    assert not result.forward_primer


def test_boundary_and_out_of_range_fail():
    assert design_primers(CDS, "M1V").status == "失败"
    assert "超出" in design_primers(CDS, "A999V").message


def test_synonymous_substitution_supported():
    result = design_primers(CDS, "A10A")
    assert result.status == "成功"
    assert translate_codon(result.optimized_codon) == "A"


def test_excel_round_trip(tmp_path):
    template = tmp_path / "template.xlsx"
    template.write_bytes(template_bytes())
    frame = read_mutations(template)
    assert list(frame.columns) == ["mutation"]
    assert frame.iloc[0]["mutation"] == "A123V"
    output = tmp_path / "results.xlsx"
    output.write_bytes(results_bytes([design_primers(CDS, "E20V", "x").to_dict()]))
    check = pd.read_excel(output)
    assert list(check["引物名称"]) == ["E20V-F", "E20V-R"]
    assert list(check.columns) == ["突变", "引物名称", "方向", "引物序列（5′→3′）", "长度（nt）", "Tm（°C）", "GC（%）", "共享区长度（bp）", "失败原因"]
    from openpyxl import load_workbook
    workbook = load_workbook(output)
    assert "A2:A3" in {str(cell_range) for cell_range in workbook["设计结果"].merged_cells.ranges}
