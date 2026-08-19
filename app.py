from __future__ import annotations

from html import escape
import re

import pandas as pd
import streamlit as st

from primer_designer.design import DesignParameters, design_primers
from primer_designer.excel_io import read_mutations, results_bytes, template_bytes
from primer_designer.sequence import normalize_cds


st.set_page_config(page_title="单点突变引物设计", page_icon="🧬", layout="wide")
st.markdown("""
<style>
:root{--ink:#0b302d;--muted:#637875;--line:#d6e3df;--paper:#ffffff;--wash:#f4f8f5;--mint:#e5f5ef;--green:#08735f;--lime:#dff186}
.stApp{background:radial-gradient(circle at 88% 3%,#f2f6df 0,transparent 29rem),#f4f8f5;color:var(--ink)}
.block-container{max-width:1220px;padding:2.2rem 2.2rem 3rem}
#MainMenu,footer{visibility:hidden}.stDeployButton{display:none}
.hero{display:grid;grid-template-columns:1fr auto;align-items:end;gap:3rem;margin-bottom:1.9rem}
.eyebrow{color:var(--green);font-size:.76rem;font-weight:750;letter-spacing:.08em;margin-bottom:.6rem}
.hero h1{font-size:3.45rem;line-height:1.02;letter-spacing:-.065em;margin:0;color:var(--ink);font-weight:850}
.hero h1 span{background:linear-gradient(transparent 62%,var(--lime) 62%);padding:0 .08em}
.hero p{color:var(--muted);font-size:1rem;line-height:1.7;max-width:670px;margin:.85rem 0 0}
.steps{display:flex;align-items:flex-start;gap:0;min-width:310px;padding-bottom:.45rem}.step{min-width:95px;color:var(--muted);font-size:.72rem}.step b{display:block;color:var(--green);font-size:.68rem;margin-bottom:.36rem}.step:not(:last-child):after{content:"";display:block;position:relative;width:42px;height:1px;background:#bfd2cc;left:47px;top:-27px}
[data-testid="stVerticalBlockBorderWrapper"]{border:1px solid var(--line)!important;border-radius:18px!important;background:rgba(255,255,255,.94)!important;box-shadow:0 20px 55px rgba(20,55,47,.08);overflow:hidden}
.panel-head{display:flex;align-items:flex-start;gap:.75rem;margin:.15rem 0 1.1rem}.panel-num{display:inline-flex;width:28px;height:28px;align-items:center;justify-content:center;border-radius:7px;background:var(--mint);color:var(--green);font-size:.68rem;font-weight:800}.panel-title{font-size:1.04rem;font-weight:750;color:var(--ink);line-height:1.25}.panel-sub{font-size:.76rem;color:var(--muted);margin-top:.2rem}
.field-help{font-size:.73rem;color:var(--muted);margin:-.25rem 0 .55rem}.required{color:var(--green);font-weight:700}
.rule-strip{display:flex;flex-wrap:wrap;gap:.42rem;margin:1rem 0 .25rem}.rule-strip span{border:1px solid var(--line);background:#fff;border-radius:7px;padding:.42rem .56rem;color:var(--muted);font-size:.68rem}.rule-strip b{color:var(--ink)}
.empty{min-height:410px;display:flex;flex-direction:column;align-items:center;justify-content:center;text-align:center;color:var(--muted)}.dna{display:flex;align-items:center;margin-bottom:1.1rem}.dna i{display:flex;align-items:center;justify-content:center;width:38px;height:38px;border:1px solid #b9d6cd;border-radius:50%;background:#fff;color:var(--green);font-style:normal;font-size:.72rem;font-weight:800;box-shadow:0 4px 12px rgba(8,115,95,.06)}.dna em{width:24px;height:1px;background:#bfd2cc}.empty strong{color:var(--ink);font-size:1rem;margin-bottom:.45rem}.empty p{max-width:470px;font-size:.78rem;line-height:1.65;margin:0}
.result-meta{font-size:.76rem;color:var(--muted);margin:-.5rem 0 .7rem}
.primer-wrap{overflow-x:auto;border-top:1px solid var(--ink);border-bottom:1px solid var(--line)}.primer-table{width:100%;border-collapse:collapse;font-size:12px}.primer-table th{background:#f5f9f7;color:var(--muted);font-size:11px;font-weight:700}.primer-table th,.primer-table td{padding:9px 7px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}.primer-table tr:last-child td{border-bottom:0}.primer-table tbody tr:hover td{background:#f8fbf9}.primer-table .mutation-cell{text-align:center;vertical-align:middle;font-weight:800;color:var(--green);background:#eaf7f2!important}.primer-table .sequence{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;letter-spacing:.02em}
.stTextArea textarea,.stTextInput input,.stNumberInput input{border-color:var(--line)!important;background:#fbfdfc!important;border-radius:9px!important}.stTextArea textarea:focus,.stTextInput input:focus,.stNumberInput input:focus{border-color:var(--green)!important;box-shadow:0 0 0 1px var(--green)!important}
[data-testid="stFileUploader"]{border-color:var(--line);border-radius:10px;background:#f8fbf9}.stButton>button,.stDownloadButton>button{border-radius:9px!important;min-height:2.7rem;font-weight:700!important}.stButton>button[kind="primary"]{background:var(--ink)!important;border-color:var(--ink)!important;color:#fff!important}.stButton>button[kind="primary"]:hover{background:var(--green)!important;border-color:var(--green)!important}.stDownloadButton>button{border-color:var(--line)!important;color:var(--green)!important;background:#fff!important}
[data-testid="stExpander"]{border:0!important;border-top:1px solid var(--line)!important;border-bottom:1px solid var(--line)!important;border-radius:0!important;background:#fff}div[role="radiogroup"]{gap:.25rem}
.fineprint{text-align:center;color:#7b8f8a;font-size:.67rem;margin-top:1rem}
@media(max-width:850px){.block-container{padding:1.25rem 1rem 2.5rem}.hero{grid-template-columns:1fr;gap:1.2rem}.hero h1{font-size:2.45rem}.steps{min-width:0}.empty{min-height:260px}.primer-table th,.primer-table td{padding:8px 6px}}
@media(prefers-reduced-motion:reduce){*{scroll-behavior:auto!important;transition:none!important}}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<header class="hero">
  <div><div class="eyebrow">一键完成定点突变引物设计</div><h1>从突变位点到 <span>可订购引物</span></h1><p>粘贴 CDS、填写氨基酸突变，自动完成位点核验、大肠杆菌密码子优化、同源臂设计与 Excel 整理。</p></div>
  <nav class="steps" aria-label="设计流程"><div class="step"><b>01</b>输入序列</div><div class="step"><b>02</b>生成引物</div><div class="step"><b>03</b>导出 Excel</div></nav>
</header>
""", unsafe_allow_html=True)

with st.container(border=True):
    input_col, result_col = st.columns([0.42, 0.58], gap="large")
    with input_col:
        st.markdown('<div class="panel-head"><span class="panel-num">01</span><div><div class="panel-title">输入设计信息</div><div class="panel-sub">序列须为完整 CDS，氨基酸编号从 1 开始</div></div></div>', unsafe_allow_html=True)
        cds_raw = st.text_area("CDS DNA 序列", height=190, placeholder=">gene_name\nATGGAACAGCTG...", help="仅接受 A、C、G、T；CDS 长度必须为 3 的倍数。")
        st.markdown('<div class="field-help"><span class="required">必填</span> · 支持纯序列或 FASTA 格式</div>', unsafe_allow_html=True)

        mode = st.radio("突变位点", ["手动填写", "Excel 批量导入"], horizontal=True)
        rows = None
        if mode == "手动填写":
            mutation_text = st.text_area("突变列表", height=95, placeholder="每行一个，例如：\nA123V\nG156D", label_visibility="collapsed")
            mutations = [item.strip() for item in re.split(r"[\n,;，；]+", mutation_text) if item.strip()]
            if mutations:
                rows = pd.DataFrame(mutations, columns=["mutation"])
            st.markdown('<div class="field-help">支持 A123V、A123-&gt;V；多个位点请换行或用逗号分隔。</div>', unsafe_allow_html=True)
        else:
            upload_col, template_col = st.columns(2)
            template_col.download_button("下载模板", template_bytes(), "突变批量导入模板.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")
            uploaded = st.file_uploader("上传 Excel", type=["xlsx"], label_visibility="collapsed")
            if uploaded:
                try:
                    rows = read_mutations(uploaded)
                    upload_col.success(f"已读取 {len(rows)} 条")
                except Exception as exc:
                    st.error(str(exc))

        with st.expander("高级参数"):
            c1, c2 = st.columns(2)
            min_flank = c1.number_input("左右同源 (bp)", 20, 40, 20)
            min_overlap = c2.number_input("共享区 (bp)", 20, 40, 20)
            c3, c4 = st.columns(2)
            target_tm = c3.number_input("目标 Tm (°C)", 55.0, 80.0, 68.0, 0.5)
            max_tm_difference = c4.number_input("最大 Tm 差", 1.0, 10.0, 5.0, 0.5)
            c5, c6 = st.columns(2)
            min_gc = c5.number_input("最低 GC (%)", 20.0, 60.0, 40.0, 1.0)
            max_gc = c6.number_input("最高 GC (%)", 40.0, 80.0, 60.0, 1.0)

        design_clicked = st.button("生成引物　→", type="primary", width="stretch")
        if design_clicked:
            if rows is None or rows.empty:
                st.error("请先填写或上传至少一个突变位点")
            else:
                try:
                    cds = normalize_cds(cds_raw)
                    params = DesignParameters(min_flank=int(min_flank), min_overlap=int(min_overlap), max_flank=max(30, int(min_flank) + 10), target_tm=float(target_tm), min_gc=float(min_gc), max_gc=float(max_gc), max_tm_difference=float(max_tm_difference))
                    st.session_state["design_records"] = [design_primers(cds, row.mutation, params=params).to_dict() for row in rows.itertuples(index=False)]
                except ValueError as exc:
                    st.error(str(exc))

    with result_col:
        st.markdown('<div class="panel-head"><span class="panel-num">02</span><div><div class="panel-title">设计结果</div><div class="panel-sub">引物名称、序列与关键理化参数</div></div></div>', unsafe_allow_html=True)
        if "design_records" not in st.session_state:
            st.markdown("""
            <div class="empty"><div class="dna"><i>A</i><em></em><i>T</i><em></em><i>G</i></div><strong>准备开始设计</strong><p>系统将核验原氨基酸，选择大肠杆菌偏好密码子，并生成满足左右同源与共享区规则的引物。</p><div class="rule-strip"><span><b>上游</b> 右侧 ≥20 nt</span><span><b>下游</b> 左侧 ≥20 nt</span><span><b>共享区</b> ≥20 nt</span></div></div>
            """, unsafe_allow_html=True)
        else:
            records = st.session_state["design_records"]
            successful_designs = sum(record["status"] == "成功" for record in records)
            st.markdown(f'<div class="result-meta">{len(records)} 个突变 · {successful_designs * 2} 条引物</div>', unsafe_allow_html=True)
            table_rows = []
            for record in records:
                if record["status"] == "成功":
                    primers = [(record["forward_name"], "上游", record["forward_primer"], record["forward_length"], record["forward_tm"], record["forward_gc"]), (record["reverse_name"], "下游", record["reverse_primer"], record["reverse_length"], record["reverse_tm"], record["reverse_gc"])]
                    for index, (name, direction, sequence, length, tm, gc) in enumerate(primers):
                        mutation_cell = f'<td rowspan="2" class="mutation-cell">{escape(record["mutation"])}</td>' if index == 0 else ""
                        table_rows.append(f'<tr>{mutation_cell}<td>{escape(name)}</td><td>{direction}</td><td class="sequence">{escape(sequence)}</td><td>{length}</td><td>{tm}</td><td>{gc}</td><td>{record["overlap_length"]}</td></tr>')
                else:
                    table_rows.append(f'<tr><td class="mutation-cell">{escape(record["mutation"])}</td><td colspan="7">{escape(record["message"])}</td></tr>')
            st.markdown('<div class="primer-wrap"><table class="primer-table"><thead><tr><th>突变</th><th>引物名称</th><th>方向</th><th>序列 (5′→3′)</th><th>nt</th><th>Tm</th><th>GC%</th><th>共享</th></tr></thead><tbody>' + "".join(table_rows) + '</tbody></table></div>', unsafe_allow_html=True)
            st.markdown('<div class="rule-strip"><span><b>上游</b> 右侧 ≥20 nt</span><span><b>下游</b> 左侧 ≥20 nt</span><span><b>共享区</b> ≥20 nt</span></div>', unsafe_allow_html=True)
            download_col, _ = st.columns([1, 1])
            download_col.download_button("导出设计结果", results_bytes(records), "引物设计结果.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")

st.markdown('<div class="fineprint">v2.1 · 设计结果用于实验前辅助判断；订购前请结合完整质粒与具体试剂盒要求复核。</div>', unsafe_allow_html=True)

