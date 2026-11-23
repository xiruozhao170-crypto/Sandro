# 报告数据与图像来源核对

核对日期：2026-11-23。对应[报告正文](Report_2026-11-23.md)及[图像附录](Report_figure_appendix_2026-11-23.md)。本文件区分原始资料、下载副本和分析方法版本；文件时间戳仅用于辅助追溯。

## 1　ACCESS-CM2 两批 zg 下载的结论

移动硬盘根目录为 `E:\Sandro Dachuang`。能够确认的两次 9 月下载分别是：

|批次|本地位置|保存内容|完成时间，北京时间|
|---|---|---|---|
|Z-sub|`ACCESS-CM2_zg_script_20260911/500hPa`|14 个 500 hPa 子集文件|2026-09-12 16:14|
|Z-full|硬盘根目录的 `zg_Amon_ACCESS-CM2_*.nc`|14 个原始全层文件，合计 18,921,284,485 字节|2026-09-13 04:06|

早期 `ACCESS-CM2_500hPa` 目录只包含 r1 historical 的先行提取/验证，后来被 Z-sub 重用，不是第三套独立气候试验。

证据来自 Z-sub 的 `completion_manifest.json`、`requested_files.json`，Z-full 的 `download_job_20260912/original_download_status.json`，本次读取的 NetCDF 文件头、抽样数值，以及仓库的 [Task 6 输入清单](../results/task6_wpsh/input_manifest.json)。逐文件证据存于 [report_data_provenance_2026-11-23.json](../results/report_data_provenance_2026-11-23.json)。

- 两批下载记录中，同名原始文件的 SHA256 一致：14/14。
- 当前全层文件字节数与下载记录一致：14/14；tracking ID 和字节数与已提交的早期 Task 6 输入清单一致：14/14。
- 每文件读取 500 hPa 的首、中、末月份，42 个完整空间切片的值和缺测位置均相同，最大绝对差为 0。
- 经纬度及解码后的全部时间坐标一致。部分子集改变了时间编码单位，例如从 `days since 1850-01-01` 改为另一基准的小时数；原始时间数字不同不代表实际日期不同。
- 子集文件的 SHA256 与全层文件自然不同，不能用这种不同误判为不同气候数据版本。

本次没有重新对 18.9 GB 原始文件计算整文件哈希，也没有逐月比较所有场。这里“SHA256 一致”特指两次下载清单记录的原始文件校验码相同；下载当时的整文件校验与本次抽样验证分开记录。

由此可以确定，这两次下载属于同源资料的不同保存形式，不能作为独立集合成员、独立模式或独立试验重复计数。如果记忆中的“两组 ACCESS”还包括 8 月 tos/pr 或后来的风场补充，则它们是不同变量资料，见第 3 节。

### 原始 zg 版本与时间范围

|试验|成员|ESGF 发布版本|文件时间范围|
|---|---|---|---|
|historical|r1i1p1f1|v20191108|185001–194912；195001–201412|
|historical|r2i1p1f1|v20191125|185001–194912；195001–201412|
|historical|r5i1p1f1|v20210607|185001–194912；195001–201412|
|ssp245|r1i1p1f1|v20191108|201501–210012|
|ssp245|r2i1p1f1|v20200303|201501–210012|
|ssp245|r5i1p1f1|v20210802|201501–210012|
|ssp585|r1i1p1f1|v20210317|201501–210012；210101–220012；220101–230012|
|ssp585|r2i1p1f1|v20200303|201501–210012|
|ssp585|r5i1p1f1|v20210802|201501–210012|

SSP585 r1 的 2015—2100 文件内部 `version` 可仍为 v20191108，而下载目录发布版本为 v20210317；两者应分别保存，不能仅以文件头日期替代发布版本。14 个文件均列入核对，其中 12 个提供至 2100 年的分析输入，两个 r1 的 2101—2300 文件未用于当前报告。

## 2　Task 6 图像究竟对应哪一批

|图像/结果|代码输入及来源|本次核对结论|
|---|---|---|
|`figures/task6_wpsh/` 六张早期 PNG|`task6_wpsh_codex.py` 读取硬盘根目录全层 zg，即 Z-full|六张图与 `task6_wpsh_outputs` 下同名 PNG 的 SHA256 全部相同；原始文件 tracking ID 与输入清单全部匹配|
|`Task6_Figure1_WPSH_5870_JJA.png`|现行 `task6_wpsh.py` 读取 `data/processed/zg500_*`|固定 5870 gpm 诊断，保留历史缺线和未来阈值退化|
|`Task6_Figure2_WPSH_5870_by_member.png`|同上|文件名含 5870，但当前模式面板实际采用等面积阈值 Tm|
|`Task6_Figure3_WPSH_equivalent_threshold.png`|同上|现行等面积阈值总览；不是旧版的统一高度偏差订正图|

当前处理子集的仓库说明指向从原始 zg 提取 500 hPa 的流程，但没有提交完整的原始字节级输入清单。因此不凭文件名宣称处理子集与某份原始文件逐字节相同。本次改为从早期 Z-full 流程的已提交年度高度场重建现行诊断：观测面积比例为 0.3176599145，九组模式 Tm 和正负位相西界均与当前 `task6_summary.json` 一致到其保存精度。见 [report_audit_2026-11-23.json](../results/report_audit_2026-11-23.json)。这为两种流程的关键数值一致性提供了独立检查，但不等于重新生成并比较现行图的每个像素。

Git 历史也表明这些是方法版本变更：`a74c7cd` 将早期平均偏差订正改为等面积阈值。此前 README 中的 `5870+Δm` 和旧西界数值不应再用于当前 Figure 2/3 的解释。报告附录保留早期固定阈值图作为方法对照，不把它们算成另一组气候样本。

## 3　其他 ACCESS 资料与方法版本

### 海温和降水：批次 A

`ACCESS-CM2_data` 含 30 个 tos/pr 文件，即 historical、SSP245、SSP585 × r1—r5 × 两变量。文件系统修改时间集中在 2026-08-25；这一时间不等于 CMIP 发布日期，也不一定等于文件最初下载日期。本次逐个读取文件头，保存 experiment、variant、grid、tracking ID、creation_date 和维度。

仓库明确记录 SSP245 r3 的 pr 原网盘副本损坏，并在计算端从 ESGF NCI v20200428 重下替换。本次只检查本地文件头，没有完整解压检验 30 个文件的数据块，也没有将本地同名 pr 与计算端替换文件做全文件比对。因此图注使用“批次 A 对应的资料体系/现行结果”，不承诺硬盘每个字节就是生成现行图的计算端副本。

最早 ACCESS 实现与后续实现曾使用不同的模态选择和去趋势方式；Git 中 `d210f0f` 对齐旧 MIROC 流程后，现行 Task 1—4 统一取 EOF1 和线性去趋势。不能把方法改变造成的图形变化误归因为重新下载了不同成员。当前报告只使用现行图，并说明方法与原任务文字之间的实际差别。

### 850 hPa 风场：批次 W

Task 7 为单独补充的全球 ua/va 风场资料，不能由 500 hPa zg 或 tos/pr 替代。版本、来源 URL、校验记录和处理参数见 `results/task7_wpsh_eof/*_manifest.json`，下载审计见 [access_download_validation.json](../results/task7_wpsh_eof/access_download_validation.json)。当前模型分析已补齐 2100 年；历史基准是 1900—2014 年，未来是 2015—2100 年。本文使用现行六组结果，不沿用较早时段版本。

## 4　MIROC6 来源与避免旧图混入的措施

指定笔记本为：

- `task1-2(MICRO6+ERSST).ipynb`：执行记录指向 MIROC6 r10 的 tos 和 pr，正/负 PDO 夏季为 62/45。
- `task3-4_regression_covariance_pattern_T4-6_T4-7_added (1).ipynb`：历史 r1—r10；未来执行记录载入 r3/r4/r6，实际情景为 SSP126/245/370/585。

文件名中的 `MICRO6` 是原笔记本命名，实际模式为 MIROC6。文件筛选的字典序会令 `r10` 排在 `r1` 前，故不能把早期单成员图自动认作 r1。

硬盘 `task1-2_output` 中旧模式图显示 44/63 个正/负位相夏季，且降水分布与指定笔记本输出不同。报告没有用这些高清旧图替换指定 notebook 图。当前 Task 1—2 的 12 张有效图均从指定 notebook 的内嵌 PNG 提取；Task 3—4 的九张有效图采用硬盘对应 subplot 目录的高清输出，已与 notebook 的面板、成员、情景和分布逐项目视核对。源文件、嵌入输出名称及 SHA256 均列入[图像清单](../results/report_figure_manifest_2026-11-23.json)。高清图与 notebook 显示图分辨率不同，并非字节相同。

## 5　选图范围与排除项

正文 17 张、图像附录 52 张，共 69 张研究图；另有导师论文 Figure 1 截图一张。组合图和其中单图同时保留是为了选版，不增加独立证据数量。

排除所有 Task 1 regioncheck/region validation 图；排除 MIROC6 T4-5 三成员平均降水合成差、MIROC6 T4-7 三成员平均的 SSP585−SSP126，以及 ACCESS T4-5 五成员平均降水合成差。Task 3 Spread 面板是成员间方差，因此保留。早期 Task 6 的平均场图是位势高度，不属于本次要求排除的降水成员平均图。

现行 ACCESS T4-7、T4-8 分别是各成员 SSP245−historical、SSP585−historical 的合成差变化。该更改见 Git `7fb16af`，旧 README 中 SSP585−SSP245、SSP585−obs 的说明已同步纠正。

## 6　核查范围与复核入口

已读取的项目材料包括四份 Task 文档、`Materials to Read.docx`、开题 `Review.docx` 和 `Proposal.docx`、两份指定 notebook、所提供的七篇独立 PDF，以及现行仓库代码/结果。`Materials to Read` 指定的顺序和章节已用于组织背景与方法；Yang 等论文的两个路径重复，只作为一篇文献。开题计划中的 ±1σ 位相筛选、部分时间窗和计划性归因任务未直接当作已经执行的方法。

本次复核完成：降水正减负恒等式、位相样本数、观测空间相关、未来与自身历史降水型相关、Task 5 共同成员峰带统计、Task 6 等面积阈值和西界、Task 7 保存结果一致性，以及图像路径和校验码。它不是全部原始模型资料的完整重算；未新增正式的趋势归因或显著性试验。

```powershell
python scripts/audit_report_results.py
python scripts/audit_report_data_provenance.py --root "E:/Sandro Dachuang"
python scripts/validate_task7_wpsh.py
python scripts/validate_report_package.py
```

原始大文件未复制进 GitHub。报告与图像可以直接在仓库浏览；需要重跑原始资料审计时，应连接移动硬盘并调整 `--root`。本次文档及笔记本来源校验记录见 [report_source_documents_2026-11-23.json](../results/report_source_documents_2026-11-23.json)。
