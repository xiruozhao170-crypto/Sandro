# Sandro: PDO–EASM 分析与预测（ACCESS-CM2 重做版）

用导师提供的新数据 **ACCESS-CM2**（替代旧 MIROC6）重跑 task1–4。
服务器磁盘有限（~38 GB 可用），旧 MIROC6 数据（~56 GB）不再下载。

## 数据

| 位置 | 内容 |
|---|---|
| `data/ACCESS-CM2_data/` | 原始 CMIP6 月数据：`tos`（海洋原生网格）、`pr`（1.25°×1.875°）；historical(1850–2014) + ssp245 + ssp585(2015–2100)，各 5 成员 r1–r5。⚠ `pr_..._ssp245_r3` 原网盘版本比特损坏，已从 ESGF（NCI 节点，v20200428）重下替换 |
| `data/ersst.v4.sst.mnmean.nc` | ERSSTv4 观测 SST（NOAA PSL，1854–2020） |
| `data/precip.mon.total.v2018.nc` | GPCC v2018 观测月降水（0.5°，1891–2016，mm/月；导师 Google Drive 共享） |
| `data/processed/tos_2deg_*.nc` | tos 重网格化到 ERSST 同款 2°×2° 规则网格（`regrid_tos.py`，球面 KDTree IDW k=4） |

数据源：南大网盘 `box.nju.edu.cn/d/bb53a9fd82f84563a4c9`（Seafile；列目录用
`/api/v2.1/share-links/<token>/dirents/`，下载用 `/d/<token>/files/?p=%2F<file>&dl=1`）。

## 环境

```bash
venv/bin/python  # xarray, netCDF4, eofs, cartopy, scipy, matplotlib, dask
```

## 流水线

```bash
venv/bin/python scripts/regrid_tos.py            # tos 原生网格 -> 2°（15 个文件，幂等）
venv/bin/python scripts/task1_pdo.py ersst       # Task 1 观测 PDO
venv/bin/python scripts/task1_pdo.py access_r1   # Task 1 模式 PDO（r1..r5 同理）
venv/bin/python scripts/task1_compare.py         # 成员 vs 观测：方差比例、空间型相关
venv/bin/python scripts/task1_r4_modecheck.py    # r4 模态互换检查（附录）
venv/bin/python scripts/task2_composite.py       # Task 2 观测+5 成员 位相合成（依赖 task1 无）
venv/bin/python scripts/task3_spread.py          # Task 3 spread（依赖 task2 输出）
venv/bin/python scripts/task4_future.py          # Task 4 未来情景（依赖 task3_summary.json）
venv/bin/python scripts/task4_detrend_sensitivity.py  # Task 4 降水去趋势敏感性（依赖 task2/task4 输出）
```

共享模块：`pdo_phase.py`（EOF+PDO 模态选择+9 年平滑定位相）、`easm_rain.py`（华东 JJA 降水异常与合成）、`easm_plots.py`（三联图/时间序列图；注意 cartopy `gridlines(draw_labels=True)` 与 constrained_layout+colorbar 不兼容，用手动 set_xticks）。

## Task 1 结果（1900–2014，20°N–60°N / 110°E–250°E，√cos(lat) 加权，去趋势+去年循环，EOF1）

- **ERSSTv4**：EOF1 方差 20.5%，空间型（中西部冷芯 + 北美沿岸暖马蹄）与 NOAA 官方 PDO 一致。
- **ACCESS-CM2 各成员**：EOF1 方差 16.1–21.2%；与观测空间型相关 r1=0.84、r2=0.76、r3=0.76、r5=0.84。
- **⚠ r4 模态互换**：r4 的 EOF1 是黑潮延伸体（KOE）主导模态（与观测 PDO 相关仅 0.21），
  经典 PDO 出现在其 **EOF2**（|r|=0.85）；EOF1/2 方差接近（17.3% vs 14.0%，North 准则近简并）。
  下游任务若用"PC1=PDO 指数"须注意此成员。
- 模式普遍在 KOE 锋区变率偏强（图上深蓝小尺度中心），属 CMIP6 常见偏差，非插值伪影。

结果：`results/task1_pdo_*.nc`、`results/task1_summary.json`；图：`figures/task1_*.png`。

## Task 2 结果（PDO 位相合成华东 15–45N/100–130E JJA 降水，9 年滑动平均定位相）

- 方法：PC 108 月中心滑动平均 → 每年 JJA 均值定正/负位相 → 正负位相 JJA 降水异常（mm/day）合成差值，Welch t 检验 p<0.1 打点。
- **观测**（ERSSTv4 + GPCC）：pos=53 年 / neg=54 年（与导师示例完全一致）。合成型：PDO+ 华北/黄淮显著偏旱、长江中下游偏涝，与 Yu 2013 图 2d 一致。
- 模式成员 PDO 模态自动对观测型选择：r1–r3、r5 用 EOF1，**r4 自动选 EOF2**（模态互换成员）。
- 文档按 MIROC6（10 成员）写，实际按导师新数据 ACCESS-CM2（5 成员）执行。

## Task 3 结果（5 成员 spread）

- 输出跨成员 ensemble mean + 方差图（PDO 型、pos/neg/diff 合成型）。
- 各成员合成型 vs 观测合成型（GPCC 插值到模式格点，陆地点）空间相关：
  r1=-0.10, r2=-0.08, **r3=0.16, r5=0.11, r4=0.07**（低相关属预期：自由耦合模式内部变率位相与观测不同步）。
- **best3 = r3, r5, r4**（按符号排序，负相关不入选）→ 供 Task 4。

## Task 4 结果（未来情景，best3 成员，2015–2100）

- SST 用二次去趋势（ssp 增暖非线性），PDO 模态仍对观测历史型选择：
  ssp245: r3=EOF1(0.64), r5=EOF3(0.60), r4=EOF1(0.67)；
  ssp585: r3=EOF1(0.64), r5=EOF2(0.84), r4=EOF2(0.85)——增暖越强，PDO 越偏离 EOF1。
- 每情景输出：PDO+/PDO− 各一图（3 成员+ensemble mean）、PDO 时间序列、
  `task4_scenario_comparison.png`（historical vs ssp245 vs ssp585 3×3 对比）。
- ensemble 合成差值与 historical 型的空间相关（task4_summary.json）：情景间型态差异明显，
  PDO–降水关系在强迫增强下不稳定（降水异常未去趋势，遵循任务文档方法——ssp585 下残余湿趋势会混入合成，见下节敏感性检验）。

### 降水去趋势敏感性（`task4_detrend_sensitivity.py`）

- 基线方法只减全时段 JJA 平均；此检验额外对 JJA 降水逐格点多项式去趋势
  （historical 一次、ssp 二次，与 SST 处理一致），PDO 位相年份直接复用 task2/task4 已存结果。
- ensemble 合成差值 原版 vs 去趋势版 空间相关：**historical 1.00、ssp245 0.98、ssp585 0.87**
  （逐成员最低 ssp585 r4=0.81）——历史与中等排放下结论对去趋势不敏感;
  ssp585 下强迫湿趋势确实混入合成（华南/台湾附近幅度被放大），但型态大体不变。
- 去趋势版情景 vs historical 相关（-0.02 / 0.05）与原版（-0.07 / -0.08）同样接近 0：
  "未来 PDO–降水关系与历史不同"这一结论不是降水趋势造成的伪影。
- 输出：`results/task4_detrend_*.nc`、`results/task4_detrend_summary.json`、
  `figures/task4_detrend_sensitivity.png`（3×3：原版 / 去趋势 / 两者之差）。

## 备注

- Task 1–4 文档与文献已在 `docs/`（来自导师 Google Drive `1hFPx8SuuVn399gsVioY4fCWsB1HON10F`）。
- 旧 MIROC6 版本示例输出在同一 Drive（task1-2_output 等文件夹）可作图式参照。
