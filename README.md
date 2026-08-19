# 大肠杆菌单点突变引物设计工具

本地 Streamlit 应用：输入 CDS 和 `A123V` 格式的氨基酸替换，生成经过 E. coli K-12 密码子优化的定点突变引物，并支持 Excel 批量处理。

## 启动

```powershell
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

浏览器打开终端显示的本地地址。所有序列计算均在本机进行。

## Excel 格式

上传文件只需包含 `mutation` 列。可直接在页面下载模板。

## 说明

该工具使用简化的局部 Tm、GC、发卡和二聚体筛选，适合方案初筛。实验前仍应结合完整质粒序列、聚合酶/试剂盒要求和实验条件复核。
