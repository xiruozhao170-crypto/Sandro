# Sandro: PDO–EASM 分析（ACCESS-CM2，按旧 MIROC6 流程重做）

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
  Task 6 的未来谱使用指定成员 r2/r3/r5。各成员合成差与全体成员平均合成差的空间相关
  （旧 notebook `select_best_members` 的排名标准）以及与观测的相关仍作诊断输出。
- **降水只保留陆地**：模式降水插值到 GPCC 网格后套用 GPCC 陆地掩膜（GPCC 为站点资料，
  海洋格点恒为 NaN），所有降水图和空间相关都只含陆地格点，海洋在图上涂浅灰。

## 数据

| 位置 | 内容 |
|---|---|
| `data/ACCESS-CM2_data/` | 原始 CMIP6 月数据：`tos`（海洋原生网格）、`pr`（1.25°×1.875°）；historical(1850–2014) + ssp245 + ssp585(2015–2100)，各 5 成员 r1–r5。⚠ `pr_..._ssp245_r3` 原网盘版本比特损坏，已从 ESGF（NCI 节点，v20200428）重下替换 |
| `data/ersst.v4.sst.mnmean.nc` | ERSSTv4 观测 SST（NOAA PSL，1854–2020） |
| `data/precip.mon.total.v2018.nc` | GPCC v2018 观测月降水（0.5°，1891–2016，mm/月） |
| `data/processed/tos_2deg_*.nc` | tos 重网格化到 ERSST 同款 2°×2° 网格（`regrid_tos.py`，球面 KDTree IDW k=4） |

数据源：南大网盘 `box.nju.edu.cn/d/bb53a9fd82f84563a4c9`；任务文档与旧 MIROC6 示例输出在导师
Google Drive `1hFPx8SuuVn399gsVioY4fCWsB1HON10F`。

## 流水线

```bash
venv/bin/python scripts/regrid_tos.py            # tos 原生网格 -> 2°（15 个文件，幂等）
venv/bin/python scripts/task1_pdo.py ersst       # Task 1 观测 PDO（access_r1..r5 同理）
venv/bin/python scripts/task1_compare.py         # 成员 vs 观测：方差比例、空间型相关
venv/bin/python scripts/task2_composite.py       # Task 2 观测 + 5 成员位相合成
venv/bin/python scripts/task3_spread.py          # Task 3 spread（r2/r3/r5 供 Task 6 未来谱）
venv/bin/python scripts/task4_future.py          # Task 4 未来情景（ssp245/ssp585 × r1–r5）
venv/bin/python scripts/task5_ssp_minus_obs.py   # Task 5 SSP 降水气候态 − GPCC 观测
venv/bin/python scripts/task6_spectrum.py        # Task 6 PC1 傅立叶功率谱 / cycle+energy
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

图：`figures/task1_pdo_*.png`（EOF1 + PC1）、`task1_regioncheck_*.png`（区域抽查）、
`task1_pdo_members.png`（成员总览）；结果：`results/task1_pdo_*.nc`、`task1_summary.json`。

## Task 2（PDO 位相合成华东 JJA 降水，9 年滑动平均 + JJA 定位相）

- **观测**（ERSSTv4 + GPCC）：**pos = 53 年 / neg = 54 年，与导师示例完全一致**。
  合成型：PDO+ 时华北/黄淮显著偏旱、江南-华南偏湿差异显著（白点 p<0.1），与 Yu 2013 图 2d 一致。
- 模式成员（全部用 EOF1）：r1 62/45、r2 57/50、r3 54/53、r4 55/52、r5 50/57（pos/neg 年数）。
- 图：`task2_pdo_timeseries.png`（观测+5 成员，红蓝柱 + 9 年滑动黑线）、
  `task2_<tag>_three_panels.png`（正位相 / 负位相 / 合成差三联图）、
  `task2_composite_obs_members.png`（obs + r1–r5 合成差 2×3 总览，共用色标，
  白点 p<0.1；`task2_composite_grid.py` 从已存 results 出图）。

## Task 3（5 成员 spread，图式沿用旧版 Task3_Figure1–4）

- `Task3_Figure1_PDO_EOF1_Spread.png`：r1–r5 EOF1 covariance 型 + 跨成员方差
  （方差极大值在黑潮延伸体锋区——成员间 PDO 型差异主要来自 KOE）。
- `Task3_Figure2_PC1.png`：r1–r5 PC1 + 9 年滑动平均。
- `Task3_Figure3_Rainfall_Composite_Spread.png`：各成员合成差 + 跨成员方差，标注与 5 成员平均
  合成差的空间相关（旧 notebook 的选择标准，改陆地掩膜后只算陆地格点）：
  **r1 = 0.45、r2 = 0.58、r3 = 0.25、r4 = 0.47、r5 = 0.77**。
  诊断参考——与观测合成型的相关：r1 = −0.17、r2 = −0.06、r3 = 0.20、r4 = −0.23、r5 = 0.12
  （自由耦合模式内部变率与观测不同步，低相关属预期）。
- `Task3_Figure4_Positive_Negative_Rainfall_Spread.png`：正/负位相降水异常的跨成员方差。
- **Task 4 使用全部成员 r1–r5；Task 6 未来谱成员 = r2, r3, r5（指定）**。

## Task 4（未来情景 ssp245 / ssp585，全部成员 r1–r5，2015–2100，方法与历史完全相同）

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

- 图 `Task4_T4-1`–`T4-7`（成员图均为 2×5：行 = 情景，列 = r1–r5）：未来 PDO EOF1、PC1、
  正/负位相降水、逐成员合成差、5 成员平均合成差（historical vs ssp245 vs ssp585，
  `Task4_T4-5_Composite_Difference_MemberMean.png`）、SSP585−SSP245 情景差。
- 5 成员平均合成差与历史型的空间相关（陆地格点）：ssp245 = −0.24、ssp585 = −0.00；
  两情景之间 −0.07（`task4_summary.json`）——**PDO–华东夏季降水关系在未来情景下
  与历史型差异明显，且情景间不一致**。
- 注意：ssp585 下强迫增暖强烈，线性去趋势后 EOF1 不一定是经典 PDO 型
  （如 ssp585 r3 = −0.40、r4 = −0.30，与观测型相关为负）；这是任务文档规定方法
  （EOF1 + 线性去趋势）下的真实结果，解释未来合成图时需留意。

## Task 5（SSP 情景降水 − 历史观测，全部 5 成员，陆地）

JJA 降水气候态（mm/day，GPCC 0.5° 网格、陆地掩膜）：各情景 2015–2100 减 GPCC 观测 1900–2014。

- `Task5_Figure1_SSP_minus_Obs_MemberMean.png`：观测气候态 | ssp245−obs | ssp585−obs（5 成员平均）。
- `Task5_Figure2_SSP_minus_Obs_by_member.png`：2×5 逐成员差值。
- 区域平均（陆地）：**ssp245 = +1.35 mm/day，ssp585 = +1.45 mm/day**（`task5_summary.json`）；
  空间上江南−华南显著偏湿（>+5 mm/day），朝鲜半岛与华北略偏干。
- 注意该差值同时包含**模式气候态偏差**与**未来增暖增湿信号**，不能全部解释为气候变化响应。

## Task 6（PC1 傅立叶功率谱：cycle 周期与 energy，方法按导师论文 Fig 1c）

方法 = 导师论文（Atmosphere 2020, 11, 3, doi:10.3390/atmos11010003）图 1c：
对**月分辨率归一化 PC1** 做 periodogram（傅立叶变换的功率谱），叠加由 lag-1 自相关拟合的
**AR1 红噪声理论谱**（Gilman et al. 1963）及其 **90%/95% χ² 置信线**；谱做 5 点 Daniell 平滑
（dof≈10，原始 periodogram 每 bin 仅 2 dof、单点噪声尖峰会假超线）。

- **cycle**：搜索窗口 [2 年, 记录长度/3] 内超过置信线的最强谱峰对应的周期
  （更长周期在记录内不足 3 个循环、与趋势不可分；1 年整峰为残余年循环，均排除）。
- **energy**：谱峰附近连续超过红噪声背景的频段上对 PSD 积分（Parseval：全谱积分 = PC1 总方差
  = 1 σ²），同时给出其占总方差的比例。

| 序列 | cycle（年） | 显著性 | energy（σ²，占方差） |
|---|---|---|---|
| obs (ERSSTv4) | 28.8 | 95% | 0.23（23%） |
| r1 / r2 / r3 / r4 / r5 | 10.5 / 38.3 / 38.3 / 38.3 / 28.8 | 90–95% | 0.18–0.42（18–41%） |
| ssp245 r2 / r3 / r5 | 28.7 / 28.7 / 2.2 | 不显著 / 90% / 95% | 0.07–0.19（7–19%） |
| ssp585 r2 / r3 / r5 | 28.7 / 28.7 / 28.7 | 95% / 不显著 / 不显著 | 0.07–0.15（8–17%） |

- 图：`Task6_Figure1_Historical_PC1_Spectrum.png`（obs + r1–r5）、
  `Task6_Figure2_Future_PC1_Spectrum.png`（ssp245/ssp585 × r2/r3/r5）、
  `Task6_Figure3_Cycle_Energy_Summary.png`（周期 + 能量汇总条形图）；
  数值：`results/task6_summary.json`。
- **主要结论**：历史时期（观测与 5 成员）PDO PC1 的年代际–多年代际 cycle（约 10–38 年）
  能量占总方差 18–41%（观测 23%）；未来情景下该 cycle 能量降到 7–19%
  （ssp245 平均 12%、ssp585 平均 13%），**PDO 低频循环显著减弱、变率向高频移动**——
  这正是 Task 3/4 中未来 PC1 子图彼此差异大、位相结构变乱的量化解释。
- ssp585 谱中 1 年整的尖峰为强迫增暖下残余年循环，已从 cycle 搜索中排除。

## 备注

- Task 1–4 文档在 `docs/`；旧 MIROC6 版示例输出（Drive `task1-2_output`、`task3/4_*_outputs`）为图式参照。
  导师论文 `atmosphere-11-00003-v2.pdf`（Task 6 方法来源）与其重点图（`微信图片_*.png`，
  即论文 Figure 1，图 c 为功率谱模板）在仓库根目录。
- 本轮改动：Task 4 改为**全部 5 个成员 r1–r5**（2×5 图版式），不再选 best members；
  Task 6 未来谱暂仍用 r2/r3/r5。
- 前轮改动：①同一幅图内的 colorbar 统一（`task1_pdo_members` 六面板共用一套色标、
  task2 三联图 (a)(b)(c) 共用一套色标）；②Task 4/6 成员曾改为指定 r2/r3/r5
  （此前按"vs 全体成员平均合成差相关"排名取 r5/r2/r4）。
- 上一轮改动：所有降水图改为只保留陆地（GPCC 掩膜）并重跑 Task 2–4；
  新增 Task 5（SSP−观测降水差）与 Task 6（PC1 功率谱 cycle/energy）。
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
