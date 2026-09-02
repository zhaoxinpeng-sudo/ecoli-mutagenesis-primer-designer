from __future__ import annotations

import pandas as pd

from .excel_io import _write_styled_excel

MAP_COLUMNS = ["样品编号", "目标突变"]


def verification_template_bytes() -> bytes:
    return _write_styled_excel(pd.DataFrame([["A01", "A123V"]], columns=MAP_COLUMNS), "测序映射")


def read_verification_map(file) -> pd.DataFrame:
    frame = pd.read_excel(file, dtype=str)
    aliases = {"sample_id": "样品编号", "mutation": "目标突变"}
    frame = frame.rename(columns=aliases)
    missing = [column for column in MAP_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError("Excel 缺少必需列：" + "、".join(missing))
    frame = frame[MAP_COLUMNS].fillna("")
    for column in MAP_COLUMNS:
        frame[column] = frame[column].str.strip()
    frame = frame[(frame["样品编号"] != "") & (frame["目标突变"] != "")].drop_duplicates().reset_index(drop=True)
    if frame.empty:
        raise ValueError("Excel 中没有可用的样品映射")
    return frame


def verification_results_bytes(results) -> bytes:
    columns = ["样品编号", "目标突变", "目标氨基酸", "实际氨基酸", "预期密码子", "实际密码子", "判定", "有效读段数", "测序方向",
               "参考覆盖率（%）", "目标最低Q值", "覆盖范围", "额外变异", "警告/原因"]
    frame = pd.DataFrame([result.to_dict() for result in results], columns=columns)
    return _write_styled_excel(frame, "测序验证结果")
