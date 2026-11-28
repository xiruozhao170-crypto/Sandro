# Sandro: PDO–EASM 分析（ACCESS-CM2，按旧 MIROC6 流程重做）

## 完整研究报告（2026-11-23）

[报告正文](docs/Report_2026-11-23.md) · [图像附录](docs/Report_figure_appendix_2026-11-23.md) · [数据批次及图像来源](docs/Report_data_provenance_2026-11-23.md)

新版按 MIROC6 Task 1—4、ACCESS 补充分析及 Task 5—7 重写，正文与附录收录 65 张研究图和 1 张方法参考图，列出 22 条参考文献。已按要求排除区域检查图及 Task 3—4 成员平均降水图。旧初稿由本版替代。

两批 9 月 ACCESS zg 下载已核对为同源的 500 hPa 子集与全层文件；具体图像对应关系、验证范围及原始文件版本见来源附录。当前报告依据代码与保存结果，修正了旧说明中的 T4-7/T4-8 差分对象、Task 5 显著性解释和 Task 6 等面积阈值方法。

2026-11-25 清理：已从仓库删除 ACCESS T4-5 成员平均降水图及其生成逻辑，移除 Task 3 降水均值字段和 Task 3—4 依赖降水成员平均的统计。降水分析保留各成员正/负位相、合成差、同一成员的未来−历史比较及成员间方差；Task 3 图中相关系数改为各成员与观测的相关。MIROC6 T4-5/T4-7 平均降水图此前已排除。Task 6 位势高度资料不属于此次清理范围。

## 新增 Task 7：850 hPa 风场 → 流函数 → WPSH EOF 指数

按 [Yang et al. (2022), PNAS](https://www.pnas.org/doi/10.1073/pnas.2120335119)
计算 **1979–2020 年 ERA5 全球 0.25°、850 hPa u/v 风场**的流函数，取
**10–40°N、110–180°E** 内的 JJA 流函数距平，线性去趋势后做面积加权 EOF。
标准化 PC1 为 WPSH 指数，正值表示副高增强。另提供 NCEP/NCAR 再分析对照。

| 资料 | EOF1 方差解释率 | EOF2 | PC1 与区域流函数相关 | PC1 > 1.5 SD 年份 |
|---|---:|---:|---:|---|
| ERA5（主结果） | **77.86%** | 9.12% | **0.953** | **1995、1998、2010、2020** |
| NCEP/NCAR R1 | 75.85% | 9.09% | 0.957 | 1988、1998、2010、2020 |

[ERA5 图](figures/task7_wpsh_eof/era5_1979_2020.png) ·
[逐年指数 CSV](results/task7_wpsh_eof/era5_1979_2020_index.csv) ·
[完整 NetCDF](results/task7_wpsh_eof/era5_1979_2020.nc) ·
[方法、敏感性与复现步骤](docs/Task7_WPSH_streamfunction_EOF.md) ·
[验证记录](results/task7_wpsh_eof/validation_raw.json)

已通过解析风场、独立 SVD 和原始 ERA5 月风场复核。主结果采用 JJA 去趋势；
同时保存“先月距平去趋势、再取 JJA”及不做纬度加权的敏感性结果。
原 Task 1–6 保留。另完成 **ACCESS-CM2 的 r1、r2、r5 × SSP245、SSP585 六组
1900–2100 年分析**：全球风场反演流函数，逐月二次去趋势，JJA EOF，基准期为
1900–2014 年（115 年），未来比较期为 2015–2100 年（86 年），与历史试验边界一致。
2100 年风场已从同一版本的原始资料补齐，并纳入完整时段的去趋势和 EOF 重算。
强事件采用固定历史期均值和标准差定义；由于时段长度不同，图中比较发生率（强事件夏季占比），
同时在表中保留实际次数。高地下方缺测采用球面调和延拓，每组均有最近邻延拓对照。
ERA5 人工遮挡试验中，调和延拓 PC1 与完整风场结果相关 **0.99986**，四个强事件
年份保持一致；这检验了补值对指数的影响，并不意味着恢复了真实地下风场。

[六组结果汇总](docs/Task7_ACCESS_results.md) ·
[模式指数总览](figures/task7_wpsh_eof/access_overview.png) ·
[历史与未来比较](figures/task7_wpsh_eof/access_period_comparison.png) ·
[模式汇总 CSV](results/task7_wpsh_eof/access_summary.csv) ·
[补值检验](results/task7_wpsh_eof/wind_fill_benchmark.json)

```bash
python -m pip install -r requirements-wpsh.txt
python scripts/fetch_wpsh_winds.py era5       # 需已配置 CDS API
python scripts/task7_wpsh_eof.py --source era5
python scripts/fetch_access_wpsh_winds.py --workers 4
python scripts/task7_wpsh_eof.py --source access --wind-fill harmonic
python scripts/task7_wpsh_eof.py --source access --wind-fill nearest
python scripts/summarize_access_wpsh.py
python scripts/task7_wpsh_eof.py --redraw     # 仅用提交的结果重画，无需原始风场
python -m pytest tests/test_wpsh_streamfunction.py tests/test_wpsh_wind_fill.py -q
python scripts/validate_task7_wpsh.py
```

用导师提供的新数据 **ACCESS-CM2**（5 成员 r1–r5，替代旧 MIROC6 10 成员）重跑 Task 1–4。
方法严格沿用旧数据版本的代码思路（任务文档 Task1–4.docx 的规定流程），
不做任何额外的"改进"（上一版曾引入 EOF 模态自动选择和二次去趋势，本版已全部去掉）：

- **PDO 一律取 EOF1**（Task 1.docx：PDO = 北太平洋 SST 异常的第一 EOF 模态）；
- **SST 一律逐格点线性去趋势**（历史与未来情景相同处理）；
- 位相判定：PC1 做 **9 年（108 个月）中心滑动平均**，每年用**夏季（JJA）均值**的正负号分正/负位相（Task 2.docx）；
- 降水合成：华东 15–45°N / 100–130°E，每月总量 ÷ 当月天数得 mm/day，JJA 三个月平均成季节值，
  减全时段季节平均得异常（与旧 notebook 完全相同的定义）；正位相平均 − 负位相平均 = 合成差；
  模式降水统一插值到 **GPCC 0.5° 网格**（旧版 npy 的 "ec" 网格）。
- Task 4 对**全部 5 个成员 r1–r5** 跑未来情景（不再选 best members）；
  Task 5 的未来谱使用指定成员 r2/r3/r5。各成员降水合成差与观测的空间相关作为诊断输出，
  不再计算降水成员平均及与其相关的统计。
- **降水只保留陆地**：模式降水插值到 GPCC 网格后套用 GPCC 陆地掩膜（GPCC 为站点资料，
  海洋格点恒为 NaN），所有降水图和空间相关都只含陆地格点，海洋在图上涂浅灰。

## 数据

| 位置 | 内容 |
|---|---|
| `data/ACCESS-CM2_data/` | 原始 CMIP6 月数据：`tos`（海洋原生网格）、`pr`（1.25°×1.875°）；historical(1850–2014) + ssp245 + ssp585(2015–2100)，各 5 成员 r1–r5。⚠ `pr_..._ssp245_r3` 原网盘版本比特损坏，已从 ESGF（NCI 节点，v20200428）重下替换 |
| `data/ersst.v4.sst.mnmean.nc` | ERSSTv4 观测 SST（NOAA PSL，1854–2020） |
| `data/precip.mon.total.v2018.nc` | GPCC v2018 观测月降水（0.5°，1891–2016，mm/月） |
| `data/processed/tos_2deg_*.nc` | tos 重网格化到 ERSST 同款 2°×2° 网格（`regrid_tos.py`，球面 KDTree IDW k=4） |
| `data/processed/zg500_*.nc` | 500 hPa 位势高度子集（`extract_zg500.py` 从原始 zg 抽取，60–240°E / 15°S–70°N）：ACCESS-CM2 historical(1850–2014)+ssp245+ssp585(2015–2100) × r1/r2/r5（导师 Google Drive `1JNrGqpUA6vcpK0Zv5AOU-ORThdWAJSZj`，原始文件每个 1–1.6 GB，抽取后即删）；`zg500_ncep_...` = NCEP/NCAR R1 月平均 hgt 500 hPa（NOAA PSL，1948–今） |

数据源：南大网盘 `box.nju.edu.cn/d/bb53a9fd82f84563a4c9`；任务文档与旧 MIROC6 示例输出在导师
Google Drive `1hFPx8SuuVn399gsVioY4fCWsB1HON10F`。

## 流水线

```bash
venv/bin/python scripts/regrid_tos.py            # tos 原生网格 -> 2°（15 个文件，幂等）
venv/bin/python scripts/task1_pdo.py ersst       # Task 1 观测 PDO（access_r1..r5 同理）
venv/bin/python scripts/task1_compare.py         # 成员 vs 观测：方差比例、空间型相关
venv/bin/python scripts/task2_composite.py       # Task 2 观测 + 5 成员位相合成
venv/bin/python scripts/task3_spread.py          # Task 3 spread（r2/r3/r5 供 Task 5 未来谱）
venv/bin/python scripts/task4_future.py          # Task 4 未来情景（ssp245/ssp585 × r1–r5）
venv/bin/python scripts/task5_spectrum.py        # Task 5 PC1 傅立叶功率谱 / cycle+energy
venv/bin/python scripts/extract_zg500.py <raw_zg.nc>  # zg 原始文件 -> 500 hPa 子集（逐个）
venv/bin/python scripts/task6_wpsh.py            # Task 6 WPSH 5870 gpm 位置（JJA，按 PDO 位相）
```

共享模块：`pdo_phase.py`（EOF1 PDO + 9 年平滑 JJA 定位相）、`easm_rain.py`（华东 JJA 降水异常、
GPCC 网格插值、合成 + Welch t 检验）、`easm_plots.py`（三联图/时间序列图）。

## Task 1（1900–2014，20–60°N / 110–250°E，√cos(lat) 加权，线性去趋势+去年循环，EOF1）

| 数据 | EOF1 方差 | 与观测 PDO 型相关（centered，面积加权） |
|---|---|---|
| ERSSTv4 | 20.5% | — |
| r1 | 17.0% | 0.84 |
| r2 | 16.1% | 0.76 |
| r3 | 21.2% | 0.76 |
| r4 | 17.3% | **0.21**（EOF1 为黑潮延伸体主导型，EOF1/2 方差 17.3% vs 14.0% 近简并） |
| r5 | 19.1% | 0.84 |

图：`figures/task1-4/report_access/task1_pdo_*.png`（EOF1 + PC1）、
[观测与五成员区域检查](figures/task1-4/report_access/task1_regioncheck_obs_members.png)（3 行 × 2 列，
依次为 ERSSTv4 观测和 ACCESS-CM2 R1–R5，六个子图共用 SST 色标）、
`task1_pdo_members.png`（成员总览）；
结果：`results/task1_pdo_*.nc`、`task1_summary.json`。
区域检查仍取 1900 年 8 月；该月五成员 SST 保存在 `results/task1_regioncheck_access.nc`，
观测 SST 保存在 `results/task1_regioncheck_ersst.nc`，可用 `python scripts/task1_regioncheck.py` 直接重画。
若要从原始资料重建，运行 `python scripts/task1_regioncheck.py --data-root <数据目录> --ersst-file <ERSSTv4月海温文件>`
（数据目录含 `processed/` 或 `ACCESS-CM2_data/`；两个参数也可单独使用以更新对应缓存）。

## Task 2（PDO 位相合成华东 JJA 降水，9 年滑动平均 + JJA 定位相）

- **观测**（ERSSTv4 + GPCC）：**pos = 53 年 / neg = 54 年，与导师示例完全一致**。
  合成型：PDO+ 时华北/黄淮显著偏旱、江南-华南偏湿差异显著（白点 p<0.1），与 Yu 2013 图 2d 一致。
- 模式成员（全部用 EOF1）：r1 62/45、r2 57/50、r3 54/53、r4 55/52、r5 50/57（pos/neg 年数）。
- 图：`task2_pdo_timeseries.png`（观测+5 成员，红蓝柱 + 9 年滑动黑线）、
  `task2_obs_three_panels.png`（观测三联图）、
  [五成员位相降水总图](figures/task1-4/report_access/task2_members_three_panels.png)（5 行 × 3 列；
  各行依次为 R1–R5，列为正位相 / 负位相 / 合成差；15 个子图共用色标，保留 p<0.1 白点；
  `python scripts/task2_member_panels.py` 从已存结果重画）、
  `task2_composite_obs_members.png`（obs + r1–r5 合成差 2×3 总览，共用色标，
  白点 p<0.1；`task2_composite_grid.py` 从已存 results 出图）。

## Task 3（5 成员 spread，图式沿用旧版 Task3_Figure1–4）

- `Task3_Figure1_PDO_EOF1_Spread.png`：r1–r5 EOF1 covariance 型 + 跨成员方差
  （方差极大值在黑潮延伸体锋区——成员间 PDO 型差异主要来自 KOE）。
- `Task3_Figure2_PC1.png`：r1–r5 PC1 + 9 年滑动平均。
- `Task3_Figure3_Rainfall_Composite_Spread.png`：各成员合成差 + 跨成员方差，标注各成员与观测
  合成型的空间相关 r(obs)，仅使用陆地有效格点：
  **r1 = −0.17、r2 = −0.06、r3 = 0.20、r4 = −0.23、r5 = 0.12**
  （这是按各自 PDO 位相合成后的空间型相关，不能仅用自由耦合模式与观测时间不同步解释其偏低）。
- `Task3_Figure4_Positive_Negative_Rainfall_Spread.png`：正/负位相降水异常的跨成员方差。
- **Task 4 使用全部成员 r1–r5；Task 5 未来谱成员 = r2, r3, r5（指定）**。

## Task 4（未来情景 ssp245 / ssp585，全部成员 r1–r5，2015–2100，方法与历史完全相同）

Task 1–4 的图像统一位于 [`figures/task1-4/`](figures/task1-4/)；
其中 [`report_access/`](figures/task1-4/report_access/) 存放现行 ACCESS 分析及配套观测图，
[`report_miroc6/`](figures/task1-4/report_miroc6/) 保存 MIROC6 Task 1–4 图。2026-11-29 已为 Task 3 的 EOF1、PC1、降水合成图加入观测面板，将 Task 4 四情景的五张图统一改为 r1、r5、r9，并将该目录所有降水图改为仅显示陆地。

生成代码为 [`scripts/miroc6_report.py`](scripts/miroc6_report.py)，数值及输入清单位于 [`results/miroc6_report/`](results/miroc6_report/)。历史十成员沿用已核对的缓存，ERSSTv4/GPCC 观测与未来新成员从原始数据计算。观测不参与十成员方差；保留原 PNG 文件名以维持报告链接，其中历史 PC1 图实际增加了第六行观测。

```powershell
# 安装重画与重算所需依赖
python -m pip install -r requirements-miroc6.txt
# 从已提交的数值结果重画
python scripts/miroc6_report.py
# 从本地原始资料和已核对的历史缓存重新计算
python scripts/miroc6_report.py --data-root "E:/Sandro Dachuang" --compute
python scripts/validate_miroc6_report.py
```

| 情景 | 成员 | EOF1 方差 | 与观测历史 PDO 型相关 | pos/neg 年数 |
|---|---|---|---|---|
| ssp245 | r1 | 15.8% | 0.70 | 38/40 |
| ssp245 | r2 | 15.1% | 0.42 | 41/37 |
| ssp245 | r3 | 18.5% | 0.61 | 36/42 |
| ssp245 | r4 | 16.5% | 0.51 | 39/39 |
| ssp245 | r5 | 15.5% | 0.56 | 40/38 |
| ssp585 | r1 | 16.3% | 0.42 | 30/48 |
| ssp585 | r2 | 18.0% | 0.73 | 45/33 |
| ssp585 | r3 | 18.1% | **−0.40** | 50/28 |
| ssp585 | r4 | 19.0% | **−0.30** | 45/33 |
| ssp585 | r5 | 19.3% | 0.45 | 30/48 |

- `Task4_T4-1`—`T4-4` 和 `T4-6` 包含未来 EOF、PC1、逐成员正负位相降水及合成差；原 `T4-5` 成员平均降水图已删除并停用生成，保留其余图号。
- `Task4_T4-7_Scenario_Difference_SSP245_minus_Historical.png` 为每个成员的 SSP245 合成差减自身历史合成差。
- `Task4_T4-8_Composite_Difference_SSP585_minus_Historical.png` 为每个成员的 SSP585 合成差减自身历史合成差。
- [T4-7 / T4-8 合并图](figures/task1-4/report_access/Task4_T4-7-8_Composite_Difference_SSP_minus_Historical.png)：
  2 行 × 5 列，上行为 SSP245−historical，下行为 SSP585−historical，列为 R1–R5，
  十个子图共用色标；原 T4-7、T4-8 图同时保留在新目录中。
  运行 `python scripts/task4_historical_comparison.py` 可从既有结果重画，完整 Task 4 流程也会生成此图。
- 各成员未来合成差与自身历史合成差的空间相关（陆地格点）：SSP245 为 −0.126—0.234，
  SSP585 为 −0.346—0.333；逐成员数值及同一成员两情景间的相关见 `task4_summary.json`。
  未来降水变化应逐成员比较，以展示成员间的差异。
- 注意：ssp585 下强迫增暖强烈，线性去趋势后 EOF1 不一定是经典 PDO 型
  （如 ssp585 r3 = −0.40、r4 = −0.30，与观测型相关为负）；这是任务文档规定方法
  （EOF1 + 线性去趋势）下的真实结果，解释未来合成图时需留意。

## Task 5（PDO PC1 功率谱）

借鉴 [Veiga 等（2020）](https://doi.org/10.3390/atmos11010003) Figure 1 的谱诊断思路，对未经九年平滑的标准化月度 PC1 计算周期图、五点 Daniell 平滑、AR(1) 红噪声和近似逐频点 90%/95% 阈值。原论文研究 AMM，本项目未复现其海温—风场 MCA。

历史使用 r1–r5，未来使用指定 r2/r3/r5。主峰搜索范围为 2 年至记录长度的三分之一；没有超过阈值的峰会回退到超过背景的峰，因此并非所有报告的周期都“显著”。所谓 energy 为所选峰带的统计方差，不是物理能量。

[历史谱](figures/Task5_Figure1_Historical_PC1_Spectrum.png) · [未来谱](figures/Task5_Figure2_Future_PC1_Spectrum.png) · [汇总](figures/Task5_Figure3_Cycle_Energy_Summary.png) · [数值](results/task5_summary.json)

观测主峰为 28.75 年。共同成员 r2/r3/r5 的所选峰带占比平均为历史 35.55%、SSP245 12.14%、SSP585 12.77%，但频带和记录长度不同，且 SSP245 r5 主峰为 2.21 年。这些数值不能作为 PDO 年代际变率“显著减弱”的检验；接近一年周期的尖峰成因也尚未单独验证。完整解释见报告第 6 节。

## Task 6（PDO 位相下的 500 hPa 副高轮廓）

沿用 Task 2/4 位相年份，模式成员 r1/r2/r5；历史 1900–2014、未来 2015–2100，观测为 NCEP/NCAR 1948–2014（实际位相样本 27/36）。参考 Dong（2016）绘制固定 5870 gpm，但保留全部 ENSO 状态。

现行模式主诊断使用**等面积阈值 Tm**：在 10–40°N、110–180°E 内，使各自气候态达到阈值的面积比例等于观测高于 5870 gpm 的比例（0.31766）。这替代了早期的统一偏差订正 `5870+Δm`。

|试验|r1 Tm|r2 Tm|r5 Tm|
|---|---:|---:|---:|
|historical|5856.6|5855.9|5860.1|
|SSP245|5912.2|5912.4|5918.2|
|SSP585|5932.1|5933.1|5938.9|

单位 gpm。固定阈值下历史有缺线、未来最西格点落到约 100.3°E 的诊断边界；不能把它作为精确副高西脊点。等面积阈值下，各成员的 PDO 位相西界变化仍不一致，未进行西界差的显著性检验。

- [固定 5870 gpm](figures/Task6_Figure1_WPSH_5870_JJA.png)
- [逐成员等面积阈值](figures/Task6_Figure2_WPSH_5870_by_member.png)：保留旧文件名，但实际使用 Tm。
- [等面积阈值总览](figures/Task6_Figure3_WPSH_equivalent_threshold.png)
- [当前数值](results/task6_summary.json)

### Task 6 早期本地实现

`figures/task6_wpsh/` 的六张图及对应 SVG、`results/task6_wpsh/` 和 [早期方法](docs/Task6_WPSH_method.md)保留。代码为 `scripts/task6_wpsh_codex.py`，输入是移动硬盘根目录的原始全层 zg；14 个文件列入清单，12 个用于至 2100 年的分析。这套图使用固定 5870 gpm。它与现行方法的关系及两次下载核对见[来源附录](docs/Report_data_provenance_2026-11-23.md)。

## 备注与历史变更记录

- Task 1–4 文档在 `docs/`；旧 MIROC6 版示例输出（Drive `task1-2_output`、`task3/4_*_outputs`）为图式参照。
  导师论文 `atmosphere-11-00003-v2.pdf`（Task 5 方法来源）与其重点图（`微信图片_*.png`，
  即论文 Figure 1，图 c 为功率谱模板）在仓库根目录。
- 早期改动记录（后续已改为等面积阈值）：新增 **Task 6（WPSH 5870 gpm 位置）**：导师在 Google Drive 提供 zg
  （historical/ssp245/ssp585 × r1/r2/r5，共 12 个原始文件 ~15 GB），
  逐个下载→`extract_zg500.py` 抽 500 hPa 子集（每个 ~13 MB）→删原始文件；
  观测另下 NCEP/NCAR R1 月平均 hgt（320 MB，同样只保留 500 hPa 子集）。
- 前轮改动：**原 Task 6（PC1 功率谱）更名为 Task 5**（脚本/图/结果文件与图题全部改名）——
  因原 Task 5（SSP−观测气候态差）已整体删除（脚本、图、结果均已移除）。
- 早期改动记录（后续已改为 SSP245−hist / SSP585−hist）：T4-7 曾改为逐成员 SSP585−SSP245（1×5）；新增 **T4-8** 逐成员 SSP585−obs
  合成差之差（1×5）；Task 4 改为**全部 5 个成员 r1–r5**（2×5 图版式），不再选 best members；
  功率谱任务的未来谱暂仍用 r2/r3/r5。
- 前轮改动：①同一幅图内的 colorbar 统一（`task1_pdo_members` 六面板共用一套色标、
  task2 三联图 (a)(b)(c) 共用一套色标）；②Task 4 与功率谱任务的成员曾改为指定 r2/r3/r5
  （此前按"vs 全体成员平均合成差相关"排名取 r5/r2/r4）。
- 更早改动：所有降水图改为只保留陆地（GPCC 掩膜）并重跑 Task 2–4；
  新增功率谱任务（PC1 cycle/energy，时称 Task 6）与 SSP−观测气候态差任务（时称 Task 5，已删除）。
- 相对上一版 GitHub 结果的修正：去掉了 PDO 模态自动选择（一律 EOF1）、去掉了未来情景 SST 的
  二次去趋势（一律线性）、模式降水改为插值到 GPCC 0.5° 网格再合成（与旧 MIROC6 流程一致）、
  图集版式改回旧版 Task3_Figure1–4 / Task4_T4-1–7。
- 本版按导师上传的旧 notebook（`task1-2(MICRO6+ERSST).ipynb`、`task3-4_..._added (1).ipynb`）
  重新核对后修正的三处出入：
  1. **Task 3 best3 选择标准**：改回旧标准"各成员合成差 vs 全体成员平均合成差的空间相关"
     （此前误用了与观测合成差的相关）→ best3 由 r3/r5/r2 改为 r5/r2/r4，Task 4 全部重算
     （本轮已再次改为指定 r2/r3/r5，见"本轮改动"）；
  2. **JJA 降水季节值**：改回旧定义"逐月总量÷当月天数→三个月平均"（此前是三个月总量÷92 天，
     数值差异极小但定义不同）；
  3. **位相零值**：负位相判据由 ≤0 改回旧代码的严格 <0（对结果无实际影响）。
