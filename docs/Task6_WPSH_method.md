# JJA WPSH position: adaptation of Dong (2016), Figure 4(a)

This task plots the **5870 geopotential-metre contour of the 500-hPa JJA mean height**, separately for positive and negative PDO phases. The selected ACCESS-CM2 members are **r1i1p1f1, r2i1p1f1 and r5i1p1f1**. Both choices, and retaining all ENSO states in the PDO groups, follow the user's explicit instructions on 15 September 2026.

## Figures

- [Overview](../figures/task6_wpsh/Figure4_JJA_WPSH_overview.png): observation reference, historical, SSP245 and SSP585. Model panels are equal-weight means of the three member-specific phase composites.
- [Selected members](../figures/task6_wpsh/Figure4_JJA_WPSH_selected_members.png): historical, SSP245 and SSP585 in rows; r1, r2 and r5 in columns. Each panel contains its own positive and negative PDO composites.
- [Common-period comparison](../figures/task6_wpsh/Figure4_JJA_WPSH_common_period.png): observation and the three historical members restricted to the same 1948–2014 height-data window, while retaining the original PDO indices and classifications.

Each PNG has a matching editable vector SVG. Red solid lines denote PDO positive and blue dashed lines PDO negative. The line-style convention follows the paper. There are no June-only, July-only or August-only panels.

**An absent contour is a valid result at this fixed threshold.** If a phase composite never reaches 5870 gpm in the computed region, its panel explicitly reports the absence and the regional maximum. The threshold is not lowered and no contour is invented. This occurs for several historical composites; `summary.json` records contour availability and field ranges for every dataset and phase.

## Sources and departures from the reference paper

The reference is **Dong, X. (2016), “Influences of the Pacific Decadal Oscillation on the East Asian Summer Monsoon in non-ENSO years,” Atmospheric Science Letters, 17, 115–120, [doi:10.1002/asl.634](https://rmets.onlinelibrary.wiley.com/doi/10.1002/asl.634)**. Its Figure 4(a) uses JJA positive/negative PDO composites and the 5870-gpm contour. Its other panels show individual summer months. The paper uses ERA-40 circulation data and excludes ENSO decaying summers.

**This is an adaptation with the user's data and repository grouping, not a numerical reproduction of the paper.** The user explicitly requested all years grouped by PDO, without an ENSO filter. We therefore reuse the repository's EOF1 / nine-year smoothing method, rather than substituting the paper's event lists or PDO thresholds. A fixed absolute height contour in a warming future climate reflects background height changes as well as PDO-associated differences; these maps alone do not attribute changes to PDO or establish statistical significance.

`docs/PDO_EASM_reference.pdf`, which was already in this repository, is **Yu (2013)**, a different paper. It was not used as the Figure 4 template for this task.

### Model height data

The input is the 14 `zg_Amon_ACCESS-CM2_*.nc` files directly under `E:\Sandro Dachuang`. Every file is inventoried, with its experiment, variant, monthly coverage, byte size, tracking ID and pressure selection recorded in [input_manifest.json](../results/task6_wpsh/input_manifest.json).

- Historical: two files per member, covering 1850–1949 and 1950–2014. The analysis uses 1900–2014, matching repository Tasks 1–3.
- SSP245 and SSP585: 2015–2100, matching repository Task 4.
- The r1 SSP585 files for 2101–2200 and 2201–2300 are inventoried but excluded. The existing repository PDO products and the other selected members' scenario inputs end in 2100. Extending only r1 would change both the temporal comparison and the phase definition.

Thus **all 14 files are checked, and 12 supply the requested-period composites**. No member r3 or r4, no other model and no other SSP are used. The original files remain untouched.

The selected level is verified against **50000 Pa**, tolerating only the tiny floating-point offset in the file coordinates. CMIP6 `zg` has standard name `geopotential_height` and units `m`: it is already geopotential height. It must **not** be divided by gravity again.

### Observation reference

The 14 input files contain model geopotential height only. The root-level observational SST and precipitation files cannot supply 500-hPa height. To complete the observation panel, the script downloads a small, uninterpolated **NOAA PSL NCEP/NCAR Reanalysis 1** monthly geopotential-height subset for 1948–2014, 500 hPa, 90–200°E, 0–50°N. This is explicitly labelled a reanalysis observation reference, not ERA-40 or a direct station observation.

- [NOAA dataset and access services](https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis.derived/pressure/catalog.html?dataset=Datasets%2Fncep.reanalysis.derived%2Fpressure%2Fhgt.mon.mean.nc)
- [NOAA NetCDF subset service](https://psl.noaa.gov/thredds/ncss/grid/Datasets/ncep.reanalysis.derived/pressure/hgt.mon.mean.nc/dataset.html)

The exact download request and SHA-256 checksum are in the manifest. `hgt` is already in metres of geopotential height. The observational PDO phase comes from the repository's existing ERSSTv4-derived `results/task2_obs.nc`.

## Calculation

1. **Reuse each dataset's own PDO index and years.** Read `pdo_pc` from `results/task2_obs.nc`, `task2_r1.nc`, `task2_r2.nc`, `task2_r5.nc`, and the corresponding six `task4_ssp*_r*.nc` files. Rerun `pdo_phase.pdo_index_and_phases` and assert agreement with the saved smoothed index and positive/negative year lists. File hashes make these dependencies auditable. Do not assign observed PDO years to independent model simulations.
2. **Keep the original repository PDO recipe.** North Pacific SST EOF1, pointwise linear detrending, monthly climatology removal, square-root cosine latitude weighting and the existing EOF sign convention. PDO classification uses the sign of the JJA mean of a 108-month centred moving average, with 108 required monthly values. The pre-existing PDO products are reused unchanged; no new EOF is selected or fitted in this task.
3. **Calculate complete JJA height seasons.** Extract the 500-hPa field and require June, July and August exactly once for each year. Average the three monthly mean height fields with equal weights, consistent with the repository's seasonal averaging convention. No height anomaly subtraction, detrending, bias correction, grid interpolation or contour smoothing is applied.
4. **Match by year, not by monthly timestamp.** For each height dataset, intersect its years with its own saved PDO phase years. Average JJA height fields with equal year weights within the positive group and within the negative group. All ENSO states remain eligible; missing PDO smoothing endpoints and exact-zero phases are unclassified.
5. **Contour the composite field at 5870 gpm.** Calculate the contour after averaging absolute height, not an average of annual contour positions and not a contour of the positive-minus-negative anomaly. Heights and contours remain on each dataset's native grid.
6. **For the overview only, average member composites equally.** Compute each member's positive and negative composites first, then average the three positive fields and three negative fields separately. This avoids overweighting members with more years in one phase. The selected-member figure retains all nine individual model results.

The original phase routine takes the mean of whichever smoothed JJA index values are finite. Consequently its first eligible year can contain fewer than three finite **smoothed PDO** months. This existing endpoint behaviour is retained exactly. Every **height** season nevertheless requires all three JJA months. Nominal height-data periods appear in panel subtitles; `phase_years.csv` and NetCDF attributes contain the actual included years. The nine-year smoothing reduces the eligible endpoint coverage (historical PDO: 1904–2010; future PDO: 2019–2096).

The plotted domain is 100–180°E, 0–40°N. A curve meeting an edge continues outside the map; its intersection with the frame is not necessarily the true westernmost or northernmost WPSH point. Contour GeoJSON contains the larger computed domain, with segments split across the dateline. These maps do not add significance shading, because the reference Figure 4 is a positional-contour comparison.

## Reproduce

From the repository root, using Python with `numpy`, `pandas`, `xarray`, `netCDF4`, `requests`, `eofs`, `matplotlib`, `cartopy` and `contourpy`:

```powershell
python scripts/task6_wpsh.py --data-root 'E:\Sandro Dachuang'
python scripts/validate_task6_wpsh.py --data-root 'E:\Sandro Dachuang'
```

To redraw all three figures using only the submitted NetCDF results, without downloading the original model data:

```powershell
python scripts/task6_wpsh.py --redraw
```

On this workstation Python was invoked at `D:\Program Files\anaconda3\python.exe`. Versions used are saved in the manifest. Cartopy needs its Natural Earth 110m coastline and land datasets (automatically downloaded on first use if not cached). The script checks the 14-file inventory, complete monthly source timelines, model/member identities, pressure, units, complete JJA seasons, finite fields and saved PDO phases. It records contour availability, including legitimate absent contours.

Intermediate annual height caches and the downloaded observation subset are under the git-ignored `data/` directory. All final annual JJA height arrays, positive/negative composites, differences, phase indicators and year lists are saved under [results/task6_wpsh](../results/task6_wpsh). This allows users to redraw the submitted contours without the multi-gigabyte source files. `phase_years.csv`, `summary.json`, `input_manifest.json`, `contours_5870gpm.geojson` and `validation.json` document the numerical result and checks.

## 中文说明

本次按用户确认的 **r1、r2、r5** 绘制，沿用 GitHub 原有 PDO 分组，**不剔除 ENSO 年**。仅使用 JJA 和 500 hPa 的 **5870 gpm** 等值线。主图比较观测参照与历史、SSP245、SSP585 的三成员平均；逐成员图给出全部九组模式结果；另附 1948–2014 共同时间窗的观测—历史对照。

原始 14 个文件中有两份 r1 SSP585 的 2101–2300 延伸数据，超出仓库现有 PDO 指数和其他成员的覆盖，因此登记但不加入 2015–2100 合成。观测高度由 NOAA NCEP/NCAR 再分析补充，其时间范围为 1948–2014；与原论文的 ERA-40 数据不同。合成保留绝对高度，不能将 `zg` 再除以重力加速度，也不能先求高度距平后再画 5870 等值线。全部实际合成年份、输入清单和可重画的 NetCDF 结果均随图提供。
