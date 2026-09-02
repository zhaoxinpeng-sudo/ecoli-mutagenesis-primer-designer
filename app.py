from __future__ import annotations

from html import escape
import re
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from primer_designer.design import DesignParameters, design_primers
from primer_designer.excel_io import read_mutations, results_bytes, template_bytes
from primer_designer.sequence import normalize_cds, parse_mutation
from primer_designer.sequencing import sample_id_key, unpack_files, verify_sample
from primer_designer.verification_excel import read_verification_map, verification_results_bytes, verification_template_bytes

st.set_page_config(page_title="大肠杆菌突变设计与验证", page_icon="🧬", layout="wide")
st.markdown("""<style>
:root{--ink:#0b302d;--muted:#637875;--line:#d6e3df;--mint:#e5f5ef;--green:#08735f;--lime:#dff186}.stApp{background:radial-gradient(circle at 88% 3%,#f2f6df 0,transparent 29rem),#f4f8f5;color:var(--ink)}.block-container{max-width:1220px;padding:2.1rem 2.2rem 3rem}#MainMenu,footer{visibility:hidden}.stDeployButton{display:none}.hero{display:grid;grid-template-columns:1fr auto;align-items:end;gap:3rem;margin-bottom:1.35rem}.eyebrow{color:var(--green);font-size:.76rem;font-weight:750;letter-spacing:.08em;margin-bottom:.6rem}.hero h1{font-size:3.25rem;line-height:1.02;letter-spacing:-.06em;margin:0;color:var(--ink);font-weight:850}.hero h1 span{background:linear-gradient(transparent 62%,var(--lime) 62%);padding:0 .08em}.hero p{color:var(--muted);font-size:1rem;line-height:1.7;max-width:720px;margin:.8rem 0 0}.workspace-note{font-size:.72rem;color:var(--muted)}[data-testid="stVerticalBlockBorderWrapper"]{border:1px solid var(--line)!important;border-radius:18px!important;background:rgba(255,255,255,.95)!important;box-shadow:0 20px 55px rgba(20,55,47,.08)}.panel-head{display:flex;gap:.75rem;margin:.1rem 0 1rem}.panel-num{display:inline-flex;width:28px;height:28px;align-items:center;justify-content:center;border-radius:7px;background:var(--mint);color:var(--green);font-size:.68rem;font-weight:800}.panel-title{font-size:1.04rem;font-weight:750;color:var(--ink)}.panel-sub{font-size:.76rem;color:var(--muted);margin-top:.2rem}.empty{min-height:390px;display:flex;flex-direction:column;align-items:center;justify-content:center;text-align:center;color:var(--muted)}.dna{display:flex;align-items:center;margin-bottom:1.1rem}.dna i{display:flex;align-items:center;justify-content:center;width:38px;height:38px;border:1px solid #b9d6cd;border-radius:50%;background:#fff;color:var(--green);font-style:normal;font-weight:800}.dna em{width:24px;height:1px;background:#bfd2cc}.empty strong{color:var(--ink)}.empty p{max-width:470px;font-size:.78rem;line-height:1.65}.primer-wrap{overflow-x:auto;border-top:1px solid var(--ink);border-bottom:1px solid var(--line)}.primer-table{width:100%;border-collapse:collapse;font-size:12px}.primer-table th{background:#f5f9f7;color:var(--muted)}.primer-table th,.primer-table td{padding:9px 7px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}.mutation-cell{text-align:center!important;vertical-align:middle!important;font-weight:800;color:var(--green);background:#eaf7f2!important}.sequence,.codons{font-family:ui-monospace,Consolas,monospace}.status{display:inline-flex;border-radius:999px;padding:.25rem .62rem;font-size:.72rem;font-weight:800}.pass{background:#ddf5e9;color:#08735f}.review{background:#fff3d6;color:#946200}.fail{background:#fee7e7;color:#b42318}.sample-head{display:flex;justify-content:space-between;align-items:center}.sample-title{font-weight:800;color:var(--ink)}.codons{color:var(--green);font-weight:800}.alignment{overflow-x:auto;background:#0d2725;color:#d9f5eb;padding:1rem;border-radius:10px;font:12px/1.65 ui-monospace,Consolas,monospace;white-space:pre}.stButton>button,.stDownloadButton>button{border-radius:9px!important;min-height:2.65rem;font-weight:700!important}.stButton>button[kind="primary"]{background:var(--ink)!important;border-color:var(--ink)!important;color:white!important}.fineprint{text-align:center;color:#7b8f8a;font-size:.67rem;margin-top:1rem}@media(max-width:850px){.block-container{padding:1.2rem 1rem}.hero{grid-template-columns:1fr}.hero h1{font-size:2.35rem}.empty{min-height:250px}}
</style>""", unsafe_allow_html=True)
st.markdown("""<style>
.block-container{max-width:1500px;padding:3.8rem 1.4rem 2rem}.hero{display:none}
.appbar{display:flex;align-items:center;gap:.75rem;padding:.15rem 0 .8rem;border-bottom:1px solid var(--line);margin-bottom:.75rem}.appmark{width:30px;height:30px;border:1.5px solid var(--green);border-radius:9px;display:grid;place-items:center;color:var(--green);font-weight:900}.apptitle{font-size:1.25rem;font-weight:850;letter-spacing:-.025em}.appsub{margin-left:auto;color:var(--muted);font-size:.72rem}
.workspace-shell{border:1px solid var(--line);border-radius:14px;background:#fff;box-shadow:0 12px 38px rgba(20,55,47,.06);overflow:hidden}.section-kicker{font-size:1.12rem;font-weight:800;margin:.1rem 0 .9rem}.step-label{display:flex;align-items:center;gap:.55rem;font-weight:800;margin:.2rem 0 .65rem}.step-dot{display:inline-grid;place-items:center;width:25px;height:25px;border-radius:50%;background:var(--ink);color:#fff;font-size:.75rem}.sequence-meta{display:flex;gap:.55rem;flex-wrap:wrap;color:var(--muted);font-size:.74rem;margin:-.35rem 0 .75rem}.sequence-meta b{color:var(--green)}
.summary-bar{display:flex;gap:1.6rem;align-items:center;border:1px solid var(--line);background:#f8fbf9;border-radius:10px;padding:.7rem .9rem;margin-bottom:.7rem}.summary-item{font-size:.72rem;color:var(--muted)}.summary-item b{display:block;color:var(--ink);font-size:1.05rem;margin-top:.08rem}.summary-item.good b{color:var(--green)}.summary-item.warn b{color:#a86700}
.result-toolbar{display:flex;justify-content:space-between;align-items:center;margin:.15rem 0 .65rem}.result-toolbar small{color:var(--muted)}.primer-wrap{max-height:690px;overflow:auto;border:1px solid var(--line);border-radius:9px}.primer-table{min-width:800px}.primer-table thead{position:sticky;top:0;z-index:2}.primer-table th{padding:10px 8px}.primer-table td{padding:9px 8px}.primer-table tbody tr:nth-child(4n+1) td,.primer-table tbody tr:nth-child(4n+2) td{background:#fbfdfc}.primer-table .failure-row td{background:#fff8e8!important;color:#8a5b00;white-space:normal}.mutation-cell{border-left:3px solid var(--green)}
[data-testid="stVerticalBlockBorderWrapper"]{box-shadow:none!important;border-radius:11px!important}.stButton>button[kind="primary"]{min-width:190px}.fineprint{margin-top:.75rem}
@media(max-width:850px){.block-container{padding:3.6rem .75rem 1.5rem}.appsub{display:none}.summary-bar{gap:.8rem;overflow-x:auto}.primer-wrap{max-height:520px}.apptitle{font-size:1.05rem}}
</style><header class="appbar"><div class="appmark">∿</div><div class="apptitle">突变引物设计与测序验证</div><div class="appsub">E. coli · DESIGN & VERIFY</div></header>""", unsafe_allow_html=True)
workspace = st.segmented_control("工作区", ["引物设计", "测序验证"], default="引物设计", label_visibility="collapsed")

def head(n, title, sub):
    st.markdown(f'<div class="panel-head"><span class="panel-num">{n}</span><div><div class="panel-title">{title}</div><div class="panel-sub">{sub}</div></div></div>', unsafe_allow_html=True)

def render_design():
    st.markdown('<div class="section-kicker">输入并生成引物</div>',unsafe_allow_html=True)
    with st.container(border=True):
        left, right = st.columns([.36,.64], gap="large")
        with left:
            st.markdown('<div class="step-label"><span class="step-dot">1</span>参考 CDS</div>',unsafe_allow_html=True)
            cds_raw=st.text_area("CDS 序列（5′→3′）",value=st.session_state.get("reference_cds_raw",""),height=220,placeholder=">ecoli_target_cds\nATGGAACAGCTG...",key="design_cds")
            meta='等待输入完整 CDS'
            try:
                preview=normalize_cds(cds_raw); meta=f'长度 {len(preview):,} nt　·　{len(preview)//3:,} aa　·　<b>阅读框正确</b>'
            except ValueError: pass
            st.markdown(f'<div class="sequence-meta">{meta}</div>',unsafe_allow_html=True)
            st.markdown('<div class="step-label"><span class="step-dot">2</span>填写目标突变</div>',unsafe_allow_html=True)
            mode=st.radio("录入方式",["手动填写","Excel 批量导入"],horizontal=True,label_visibility="collapsed"); rows=None
            if mode=="手动填写":
                text=st.text_area("目标突变",height=135,placeholder="每行一个，例如：\nE176C\nA123V",label_visibility="collapsed")
                values=[x.strip() for x in re.split(r"[\n,;，；]+",text) if x.strip()]
                if values: rows=pd.DataFrame(values,columns=["mutation"])
                x,y=st.columns(2); x.download_button("下载填写模板",template_bytes(),"突变批量导入模板.xlsx",width="stretch"); upload=y.file_uploader("上传突变位点表",type=["xlsx"],key="design_quick_excel",label_visibility="collapsed")
                if upload:
                    try: rows=read_mutations(upload); y.success(f"已读取 {len(rows)} 条")
                    except Exception as exc: st.error(str(exc))
            else:
                st.download_button("下载填写模板",template_bytes(),"突变批量导入模板.xlsx",width="stretch")
                upload=st.file_uploader("上传突变位点表（Excel）",type=["xlsx"],key="design_excel")
                if upload:
                    try: rows=read_mutations(upload); st.success(f"已读取 {len(rows)} 条突变")
                    except Exception as exc: st.error(str(exc))
            st.markdown('<div class="step-label"><span class="step-dot">3</span>设计参数</div>',unsafe_allow_html=True)
            with st.expander("左右同源 ≥20 bp · 共享区 ≥20 bp · Tm 68°C"):
                c1,c2=st.columns(2); flank=c1.number_input("左右同源 (bp)",20,40,20); overlap=c2.number_input("共享区 (bp)",20,40,20)
                c3,c4=st.columns(2); tm=c3.number_input("目标 Tm (°C)",55.,80.,68.,.5); diff=c4.number_input("最大 Tm 差",1.,10.,5.,.5)
                c5,c6=st.columns(2); mingc=c5.number_input("最低 GC (%)",20.,60.,40.); maxgc=c6.number_input("最高 GC (%)",40.,80.,60.)
            _,action,_=st.columns([.15,.7,.15])
            if action.button("生成引物",type="primary",width="stretch"):
                if rows is None or rows.empty: st.error("请先填写或上传至少一个突变位点")
                else:
                    try:
                        cds=normalize_cds(cds_raw); params=DesignParameters(int(flank),int(overlap),max(30,int(flank)+10),float(tm),float(mingc),float(maxgc),float(diff))
                        records=[design_primers(cds,r.mutation,params=params).to_dict() for r in rows.itertuples(index=False)]
                        st.session_state.update(design_records=records,reference_cds=cds,reference_cds_raw=cds_raw)
                    except ValueError as exc: st.error(str(exc))
        with right:
            st.markdown('<div class="section-kicker">引物设计结果</div>',unsafe_allow_html=True); records=st.session_state.get("design_records")
            if not records: st.markdown('<div class="empty"><div class="dna"><i>A</i><em></em><i>T</i><em></em><i>G</i></div><strong>结果将在这里集中显示</strong><p>输入 CDS 和目标突变后，系统会生成上下游引物并整理为可订购表格。</p></div>',unsafe_allow_html=True)
            else:
                success=sum(r["status"]=="成功" for r in records); failed=len(records)-success
                st.markdown(f'<div class="summary-bar"><div class="summary-item">突变<b>{len(records)}</b></div><div class="summary-item good">可订购引物<b>{success*2}</b></div><div class="summary-item warn">需调整<b>{failed}</b></div></div>',unsafe_allow_html=True)
                filter_col,export_col=st.columns([.72,.28]); selected=filter_col.segmented_control("筛选",["全部","设计成功","需调整"],default="全部",label_visibility="collapsed"); export_col.download_button("导出订购表",results_bytes(records),"引物设计结果.xlsx",width="stretch")
                shown=[r for r in records if selected=="全部" or (selected=="设计成功" and r["status"]=="成功") or (selected=="需调整" and r["status"]!="成功")]
                table=[]
                for r in shown:
                    if r["status"]=="成功":
                        for i,p in enumerate(((r["forward_name"],"上游",r["forward_primer"],r["forward_length"],r["forward_tm"],r["forward_gc"]),(r["reverse_name"],"下游",r["reverse_primer"],r["reverse_length"],r["reverse_tm"],r["reverse_gc"]))):
                            cell=f'<td rowspan="2" class="mutation-cell">{escape(r["mutation"])}</td>' if i==0 else ""; table.append(f'<tr>{cell}<td>{p[0]}</td><td>{p[1]}</td><td class="sequence">{p[2]}</td><td>{p[3]}</td><td>{p[4]}</td><td>{p[5]}</td></tr>')
                    else: table.append(f'<tr class="failure-row"><td class="mutation-cell">{escape(r["mutation"])}</td><td colspan="6">需调整：{escape(r["message"])}</td></tr>')
                st.markdown('<div class="primer-wrap"><table class="primer-table"><thead><tr><th>突变</th><th>引物名称</th><th>方向</th><th>序列 (5′→3′)</th><th>nt</th><th>Tm</th><th>GC%</th></tr></thead><tbody>'+''.join(table)+'</tbody></table></div>',unsafe_allow_html=True)

def trace_plot(aligned,codon_start):
    raw=aligned.raw
    if not raw.traces or not raw.peaks:return None
    indices=[aligned.ref_to_read.get(i) for i in range(codon_start,codon_start+3)]; indices=[i for i in indices if i is not None]
    if not indices:return None
    idx=indices[len(indices)//2]
    if aligned.direction=="反向":idx=len(raw.sequence)-1-idx
    idx=max(0,min(idx,len(raw.peaks)-1)); lo=max(0,idx-20); hi=min(len(raw.peaks),idx+21); x0=raw.peaks[lo]; x1=raw.peaks[hi-1]
    fig=go.Figure(); colors={"A":"#18a36b","C":"#2878d0","G":"#26332f","T":"#e24b4b"}
    for base in "ACGT":
        if base in raw.traces: fig.add_trace(go.Scatter(x=list(range(x0,min(x1+1,len(raw.traces[base])))),y=raw.traces[base][x0:x1+1],mode="lines",name=base,line=dict(color=colors[base],width=1.4)))
    fig.add_vrect(x0=raw.peaks[idx]-8,x1=raw.peaks[idx]+8,fillcolor="#dff186",opacity=.25,line_width=0); fig.update_layout(height=260,margin=dict(l=10,r=10,t=35,b=15),title=f"{raw.filename} · {aligned.direction} · 目标附近 ±20 bp",legend=dict(orientation="h"))
    return fig

def render_verification():
    with st.container(border=True):
        left,right=st.columns([.40,.60],gap="large")
        with left:
            head("01","导入测序数据","自动识别方向、合并读段并精确核对密码子")
            source=st.radio("参考 CDS 来源",["沿用引物设计 CDS","单独输入参考 CDS"],horizontal=True)
            if source=="沿用引物设计 CDS":
                reference_raw=st.session_state.get("reference_cds","")
                if reference_raw:
                    st.success(f"已载入 {len(reference_raw)} bp 参考 CDS")
                else:
                    st.info("尚无引物设计 CDS，请切换为单独输入。")
            else:
                fasta=st.file_uploader("上传 FASTA（可选）",type=["fa","fasta","fna"],key="ref_fasta"); default=fasta.getvalue().decode("utf-8-sig") if fasta else ""
                reference_raw=st.text_area("参考 CDS",value=default,height=135,placeholder=">reference\nATG...")
            st.markdown("**填写每个样品要验证的突变**")
            st.caption("先下载模板，填写样品编号和目标突变，再把表格上传回来。")
            d1,d2=st.columns(2); d2.download_button("下载填写模板",verification_template_bytes(),"测序验证映射模板.xlsx",width="stretch")
            mapfile=st.file_uploader("上传突变位点表（Excel）",type=["xlsx"],key="verify_map"); mapping=None
            if mapfile:
                try: mapping=read_verification_map(mapfile); d1.success(f"{len(mapping)} 个样品")
                except Exception as exc: st.error(str(exc))
            files=st.file_uploader("上传测序结果（AB1 / SEQ / ZIP）",type=["ab1","seq","zip"],accept_multiple_files=True,key="sequence_files")
            with st.expander("质量参数"):
                q=st.slider("AB1 质量阈值",10,40,20); minimum=st.number_input("最短有效读段 (nt)",20,300,50)
            if st.button("一键验证　→",type="primary",width="stretch"):
                try:
                    reference=normalize_cds(reference_raw)
                    if mapping is None:raise ValueError("请上传样品映射 Excel")
                    reads,warnings=unpack_files(files); grouped={}
                    for read in reads:grouped.setdefault(sample_id_key(read.sample_id),[]).append(read)
                    results=[]
                    for row in mapping.itertuples(index=False):
                        result=verify_sample(reference,str(row[0]),str(row[1]),grouped.get(sample_id_key(row[0]),[]),int(q),int(minimum)); result.warnings.extend(warnings); results.append(result)
                    st.session_state.update(verification_results=results,verification_reference=reference)
                except ValueError as exc:st.error(str(exc))
        with right:
            head("02","验证结果","预期密码子由引物设计规则自动生成"); results=st.session_state.get("verification_results"); reference=st.session_state.get("verification_reference","")
            if not results:st.markdown('<div class="empty"><div class="dna"><i>A</i><em></em><i>T</i><em></em><i>G</i></div><strong>等待测序文件</strong><p>同名 AB1/SEQ 自动配对并以 AB1 为主，测序方向由序列比对自动判断。</p></div>',unsafe_allow_html=True)
            else:
                frame=pd.DataFrame([r.to_dict() for r in results]); st.dataframe(frame[["样品编号","目标突变","目标氨基酸","实际氨基酸","预期密码子","实际密码子","判定","有效读段数","测序方向","目标最低Q值","额外变异"]],hide_index=True,width="stretch")
                st.download_button("下载中文验证报告",verification_results_bytes(results),"测序验证结果.xlsx",width="stretch")
                for r in results:
                    cls={"通过":"pass","需复核":"review","失败":"fail"}[r.verdict]
                    with st.expander(f"样品 {r.sample_id} · {r.mutation} · {r.verdict}"):
                        st.markdown(f'<div class="sample-head"><div><div class="sample-title">样品 {escape(r.sample_id)}</div><div class="codons">预期 {r.expected_codon or "—"} · 实际 {r.actual_codon or "—"}</div></div><span class="status {cls}">{r.verdict}</span></div>',unsafe_allow_html=True); st.caption(r.reason+(" · "+"；".join(r.warnings) if r.warnings else ""))
                        if r.consensus and reference:
                            _,pos,_=parse_mutation(r.mutation); start=max(0,(pos-1)*3-24); end=min(len(reference),(pos-1)*3+27); ref=reference[start:end]; con=r.consensus[start:end]; marks=''.join('|' if a==b else '^' for a,b in zip(ref,con)); st.markdown(f'<div class="alignment">坐标 c.{start+1}–{end}\n参考  {ref}\n      {marks}\n共识  {con}</div>',unsafe_allow_html=True)
                            for aligned in r.aligned_reads:
                                quality="无质量值" if aligned.qualities is None else f"平均 Q{sum(aligned.qualities)/len(aligned.qualities):.1f}"; st.markdown(f"**{aligned.raw.filename}** · {aligned.direction} · {quality} · {len(aligned.ref_to_read)} bp 比对"); fig=trace_plot(aligned,(pos-1)*3)
                                if fig:st.plotly_chart(fig,width="stretch",config={"displayModeBar":False})

if workspace=="引物设计":
    render_design()
else:
    render_verification()
st.markdown('<div class="fineprint">v3.0 · 测序判定用于实验辅助；低质量、冲突峰及关键样品请结合原始峰图人工复核。</div>',unsafe_allow_html=True)
