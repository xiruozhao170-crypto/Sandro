# Task 7: WPSH from 850-hPa winds and streamfunction EOFs

## Delivered calculation

The primary result uses **ERA5 global monthly u and v at exactly 850 hPa,
1979–2020, on the original 0.25° latitude–longitude grid**. A separate
NCEP/NCAR Reanalysis 1 calculation at 2.5° provides a cross-dataset check.
Both have 504 monthly inputs and 42 complete JJA seasons.

| Dataset | EOF1 variance | EOF2 variance | PC1 vs regional SF | PC1 > 1.5 SD |
|---|---:|---:|---:|---|
| ERA5 | 77.8554% | 9.1243% | 0.95346 | 1995, 1998, 2010, 2020 |
| NCEP/NCAR R1 | 75.8475% | 9.0939% | 0.95695 | 1988, 1998, 2010, 2020 |

The reference study reports approximately 77% and 9.3% for its two ERA5
modes and identifies 1995, 1998, 2010 and 2020 as strong events. Its observed
analysis uses detrended JJA streamfunction over 10–40°N, 110–180°E. Its
reanalysis trend is linear; its long historical-plus-future model analysis
uses a quadratic trend. See [Yang et al. (2022), PNAS,
doi:10.1073/pnas.2120335119](https://www.pnas.org/doi/10.1073/pnas.2120335119).

The close agreement and matching ERA5 event years are a useful check, but
this is an independent implementation. The paper does not document every
numerical convention, including EOF latitude weighting and all detrending
ordering details. No result was tuned to reproduce a published percentage.
The NCEP event list differs from ERA5 and is retained as calculated.

## Files

- [ERA5 figure](../figures/task7_wpsh_eof/era5_1979_2020.png) and
  [vector SVG](../figures/task7_wpsh_eof/era5_1979_2020.svg).
- [ERA5 annual index CSV](../results/task7_wpsh_eof/era5_1979_2020_index.csv).
- [ERA5 NetCDF](../results/task7_wpsh_eof/era5_1979_2020.nc): absolute JJA
  streamfunction, reference climatology, anomalies, both detrending orders,
  PC1/PC2, physical regression maps, explained variance and event flags.
- [NCEP figure](../figures/task7_wpsh_eof/ncep_1979_2020.png) and
  [annual index](../results/task7_wpsh_eof/ncep_1979_2020_index.csv).
- [Result directory](../results/task7_wpsh_eof): JSON summaries, per-input
  SHA-256 hashes, package versions and independent numerical validation.
- [Calculation](../scripts/task7_wpsh_eof.py),
  [scientific routines](../scripts/wpsh_streamfunction.py),
  [validation](../scripts/validate_task7_wpsh.py), and
  [analytic tests](../tests/test_wpsh_streamfunction.py).

The original Task 1–6 PDO, precipitation and height-contour analyses are
separate. These annual WPSH PCs are not PDO PCs and do not use the PDO
nine-year smoother or PDO phase classifications.

## Data provenance

**ERA5:** [CDS monthly averaged pressure-level dataset](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels-monthly-means?tab=overview).
The request specifies monthly averaged reanalysis, u/v wind, 850 hPa,
all 12 months, 1979–2020, 00:00, global extent and a 0.25° output grid.
For this monthly product, 00:00 identifies the monthly mean product; this
is not an extraction of instantaneous daily 00:00 winds.
The downloaded NetCDF is about 1.83 GB and is kept in ignored `data/`.
The SHA-256 in the result manifest identifies the actual downloaded input.

**NCEP:** [NOAA PSL monthly pressure-level catalogue](https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis.derived/pressure/catalog.html).
The downloader selects `uwnd`/`vwnd`, 850 hPa, 1979–2020 and the complete
73 × 144 global grid through the NOAA NetCDF subset service. Exact request
URLs and checksums are saved with the download and copied to the submitted
provenance records. Neither dataset is clipped before streamfunction inversion.

**User-provided Drive:** the [Sandro folder](https://drive.google.com/drive/folders/1hFPx8SuuVn399gsVioY4fCWsB1HON10F)
was inspected on 17 October 2026. Its ACCESS-CM2 folders contain `tos`, `pr`
and `zg` files; no `ua`/`va` files were visible. Their model/member context
was used to investigate matching public wind data, rather than treating
geopotential height or SST as wind observations.

## Calculation and conventions

1. Select the exact 850-hPa level, converting Pa to hPa where required.
   Reject unknown units, absent levels, duplicate longitude endpoints,
   incomplete or duplicate monthly records, and mismatched u/v grids.
   Sort latitude north-to-south and longitude into an increasing 0–360° convention.
2. Derive the rotational streamfunction using a **global spherical-harmonic
   inversion** of the winds, with SPHEREPACK via `skyborn.windspharm` (or
   the compatible `windspharm` backend). The radius is 6,371,200 m, with full
   resolved spectral truncation and zero degree-zero streamfunction.
   The convention is

   ```text
   u_rot = -(1/a) ∂ψ/∂φ
   v_rot =  (1/(a cosφ)) ∂ψ/∂λ
   ∇²ψ = relative vorticity
   ```

   The divergent part of the wind is not forced into ψ. A positive local
   streamfunction anomaly corresponds to anticyclonic flow in this convention.
   [Windspharm documentation](https://ajdawson.github.io/windspharm/api/windspharm.standard.html)
   describes the global vector-wind operator. A regional Poisson inversion with
   arbitrary boundary conditions would define a different calculation.
3. Average June, July and August monthly means with equal weights. Subtract
   the 1979–2020 JJA climatology. This equals forming calendar-month anomalies
   against that period and then averaging JJA. The primary observational
   result removes a least-squares **linear trend from each annual JJA gridpoint**.
4. Extract the inclusive **10–40°N, 110–180°E** EOF domain. Remove its
   temporal mean and multiply each grid column by `sqrt(cos(latitude))`.
   Compute covariance EOFs, without standardizing individual gridpoint variances.
   The leading mode is always EOF1; there is no search for a preferred mode.
5. Normalize the PCs to sample standard deviation one (`ddof=1`). Orient PC1
   so that its covariance with the area-weighted streamfunction over
   15–25°N, 120–150°E is positive. This resolves the arbitrary EOF sign and
   defines positive PC1 as stronger WPSH circulation. The sign choice is
   applied simultaneously to the PC and its map.
6. Plot regression of the detrended physical streamfunction onto each
   full-period unit-SD PC, in m²/s. The surrounding map extends to
   0–50°N, 90–200°E, while the outlined box is the actual EOF domain.
   Classify strong events using the strict inequality `PC1 > 1.5`.

The code first forms temporal linear combinations of the **global** u/v
fields, then transforms them. Global inversion, fixed climatology subtraction,
averaging and polynomial least-squares coefficients are linear operators, so
this gives the same mathematical result as transforming every monthly field.
It saves computation without changing the spatial inversion domain. Analytic
tests compare the optimized and direct monthly paths, including a quadratic
trend and a reference period shorter than the analysis period. The raw ERA5
check transforms June, July and August 1979 separately before averaging and
agrees to approximately 8.3 × 10⁻⁷ relative RMS, within transform/storage precision.

## Processing sensitivities

The separate monthly-order calculation subtracts a calendar-month climatology,
fits a single polynomial across **all monthly anomalies**, and then averages
the detrended values into JJA. That procedure need not equal detrending JJA
itself because seasonal trend differences remain after an all-month fit.
Both fields and their PC1 series are retained.

| ERA5 setting | EOF1 variance | PC1 correlation with primary |
|---|---:|---:|
| Primary: seasonal linear detrend, area-weighted EOF | 77.8554% | 1 |
| Monthly linear detrend, then JJA | 77.4060% | 0.99846 |
| Seasonal linear detrend, unweighted EOF | 77.3389% | 0.99985 |

The primary CSV uses the user's stated detrended-JJA formulation. The
monthly-first series is available as `pc1_other_detrend_order` in NetCDF;
it should not be silently substituted for the main index. Event thresholds
are sensitive to methodological choices for years close to 1.5 SD.

`pc` is normalized over the full analysis period; `pc_reference` is separately
centered and normalized over the configured present/reference period.
They coincide for these 1979–2020 reanalysis runs. Spatial result fields are
stored as float32 to keep the repository compact; indices and statistics
are stored as float64.

## ACCESS-CM2 model calculation

The public NCI catalogue contains matching r1, r2 and r5 `ua`/`va` data for
historical, SSP245 and SSP585. The separate
[download script](../scripts/fetch_access_wpsh_winds.py) reads original
OPeNDAP arrays at **85000 Pa**, preserving below-ground missing values.
A one-year r1 historical wind subset was checked against its original pressure
coordinate. The initial NCSS trial exposed a service issue: it relabelled the
85000-Pa coordinate as the mean of pressure bounds, 83125 Pa. That route is
not used in the submitted ACCESS downloader.

The model has masked 850-hPa winds below high terrain. Approximately 10% of
native grid points are missing (a grid-point fraction, not an area fraction),
including the Antarctic polar cap. The original downloaded files retain these
missing values. The low-level streamfunction solver still rejects missing winds;
model runs explicitly select the numerical wind extension described below.

The model-analysis configuration joins historical and each scenario
for **1900–2099**, uses monthly quadratic detrending over the combined period,
interpolates the resulting streamfunction to a 1° grid, and fits a common
EOF across the combined historical and future record **within each member/scenario combination**. This is
six separate EOF fits, not an EOF pooled across members or scenarios. The present
reference is **1900–2014 (115 years)**, matching the historical period in the
earlier project tasks. The future comparison is **2015–2099 (85 years)**,
matching the available scenario winds. The reference paper uses a full-record
monthly baseline for its models; the configured present-period baseline is
therefore an explicit adaptation. It must not be described as an exact
multimodel reproduction. Future WPSH indices retain the present-period SD
for event comparisons, and the full-period-normalized PC is saved separately.
Summary JSONs also report strong-event counts using the full-period SD so the
two threshold conventions can be distinguished. Historical and future years
are disjoint, with the boundary between 2014 and 2015. Because the periods
have unequal lengths, strong-event frequency is reported as the percentage
of summers, `100 * event count / number of years`; raw counts are also retained.
The standard-deviation change is `100 * (SD_future / SD_historical - 1)`.
Climatological anomalies, reference-normalized indices, thresholds and period
statistics have all been recalculated for this baseline from the original winds.
This supersedes the earlier 1900–1999 versus 2000–2099 comparison. The EOF
still spans the combined 1900–2099 record, and the ERA5/NCEP analysis remains
1979–2020. Scenario-specific detrending and EOF fits can yield slightly different
historical indices for the same member despite shared original historical winds.

Model deliverables are the six `access_{ssp245,ssp585}_{r1,r2,r5}` result sets,
the corresponding `_nearest` sensitivity results,
[combined summary](../results/task7_wpsh_eof/access_summary.csv),
[index overview](../figures/task7_wpsh_eof/access_overview.png), and
[historical/future comparison](../figures/task7_wpsh_eof/access_period_comparison.png).
Gray shading in the model index figures marks the historical reference through
2014. The comparison figure uses percentages of summers for strong-event
frequency; black crosses in its variability panel are the nearest-extension
sensitivity results, not significance markers.
All are from one model with three members, not six independent climate models.
No multimodel significance or reproduction of the paper's 32-model aggregate
is claimed. Each input manifest records original OPeNDAP URL, tracking ID,
exact pressure, file size and SHA-256.
The [download integrity audit](../results/task7_wpsh_eof/access_download_validation.json)
checks all 102 subsets (1,069,982,063 bytes), their SHA-256 checksums, 145 × 192
global coordinates, pressure level and complete monthly coverage in all 18
variable/member/experiment groups. Original pressure coordinates differ from
85000 Pa by at most 1.01 × 10⁻⁸ Pa due to numeric encoding.

### Below-ground wind extension and its limitations

The primary **harmonic extension** in
[`wpsh_wind_fill.py`](../scripts/wpsh_wind_fill.py) is an explicit numerical
choice for this independent calculation; it is not a method documented by the
reference paper. Filling missing points while holding known boundary values is
also the general principle of NCAR's
[`poisson_grid_fill`](https://www.ncl.ucar.edu/Document/Functions/Built-in/poisson_grid_fill.shtml),
but this implementation is not an NCL port or an assertion of NCL equivalence.

1. At each month, convert the valid eastward/northward winds to three Earth-fixed
   Cartesian vector components. Both wind components must have identical missing
   masks; other mismatches or infinite values fail explicitly.
2. Construct a spherical finite-volume graph with cyclic longitude. All
   longitudes at each geographic pole share a single scalar node. For regular
   angular spacings, meridional edge conductance is
   `cos(mid-latitude) * dlon / dlat`, and zonal conductance is
   `dlat / (cos(latitude) * dlon)`.
3. For each component, solve `sum_j w_ij (c_i - c_j) = 0` only at missing nodes,
   with the observed nodes fixed as Dirichlet boundaries. Check the sparse linear
   system residual. Transform back to tangent eastward/northward components and
   write only the originally missing wind points. Every original finite u/v
   value is retained exactly, and polar vector continuity is enforced for filled
   poles. There is no zero-filled wind field passed to the global inversion.
4. Fill **each month before temporal averaging or trend estimation**. Masks vary
   with time, so filling an annual mean would be a different operation. The
   later linear temporal optimization operates on already completed winds.

The separate **nearest** branch copies the Cartesian wind vector at the nearest
valid spherical node, then projects it onto the local tangent plane. It uses
the same original masks and preserves the same finite observations. Both methods
are run for all six model cases, with PC correlations, index RMSE, period
variability changes and threshold-crossing disagreements in the summary.

The [ERA5 masking experiment](../results/task7_wpsh_eof/wind_fill_benchmark.json)
regrids complete 1979–2020 ERA5 JJA winds to the native ACCESS global grid,
artificially masks them using the successive 1900–1941 ACCESS r1 JJA masks, and
compares complete, harmonic and nearest branches. Years are paired by sequence,
not by matching weather; the experiment tests the geographical masking problem.
All three branches share the same inversion, 1° interpolation, JJA linear
detrending and EOF analysis. Its correlations and errors quantify reconstruction
sensitivity on ERA5. They do not validate the unobserved below-ground model
winds or the physical accuracy of a climate projection. Model results should
be cited as conditional on the documented extension.

| ERA5 masking branch | EOF1 variance | PC1 correlation with complete | PC1 RMSE (SD) | Regional SF anomaly relative RMSE | Changed strong-event years |
|---|---:|---:|---:|---:|---|
| Complete winds on the common grid | 77.6461% | 1 | 0 | 0 | — |
| Harmonic extension | 77.7529% | 0.999860 | 0.01650 | 1.8340% | None |
| Nearest extension | 77.7094% | 0.999864 | 0.01629 | 1.7792% | None |

All three branches identify 1995, 1998, 2010 and 2020. Their common-grid complete
branch differs slightly from the original 0.25° ERA5 calculation because the
benchmark first regrids winds to the ACCESS grid. The harmonic method was chosen
before inspecting the benchmark; neither method is claimed to outperform the
other based on these small differences. Regional errors are unweighted RMS over
all years and the 10–40°N, 110–180°E EOF domain.

## Reproduction

From the repository root, in an isolated Python environment:

```bash
python -m pip install -r requirements-wpsh.txt
python scripts/fetch_wpsh_winds.py era5
python scripts/fetch_wpsh_winds.py ncep
python scripts/task7_wpsh_eof.py --source era5
python scripts/task7_wpsh_eof.py --source ncep
python scripts/fetch_access_wpsh_winds.py --workers 4
python scripts/task7_wpsh_eof.py --source access --wind-fill harmonic
python scripts/task7_wpsh_eof.py --source access --wind-fill nearest
python scripts/validate_wpsh_fill.py
python scripts/summarize_access_wpsh.py
python -m pytest tests/test_wpsh_streamfunction.py tests/test_wpsh_wind_fill.py -q
python scripts/validate_task7_wpsh.py --raw-check
```

ERA5 download requires the user's configured CDS API access and accepted
dataset terms. Credentials are read by `cdsapi`; they are not included in
the repository. NOAA downloads need no account. Raw files and computational
caches are in git-ignored `data/wpsh_winds/`.
The ACCESS downloader retries truncated OPeNDAP responses and writes its manifest
as each atomic subset completes. Rerunning resumes completed subsets. To start
each case while other downloads continue, run `python scripts/run_access_wpsh.py`
alongside the downloader; it waits for complete year coverage, then runs both
extensions and the combined summary.

The figures and numerical audits can also be reproduced without raw downloads:

```bash
python scripts/task7_wpsh_eof.py --redraw
python scripts/validate_task7_wpsh.py
```

Scientific tests cover the analytic sign and magnitude of streamfunction,
removal of divergent winds, global coverage/missing-data rejection, pressure
selection in Pa/hPa, missing and duplicate months, equivalent linear/quadratic
temporal processing, and EOF recovery against an independent SVD. The result
audits additionally verify PC normalization/orthogonality, CSV consistency,
physical regression maps and independent monthly-first raw ERA5/ACCESS checks.
The additional extension tests check unchanged valid winds, finite polar vectors,
longitude-rotation equivariance, fixed-mask linearity and correct processing of
time-varying monthly masks. Period-comparison checks cover the 2014/2015 boundary,
strict event threshold and occurrence rates for unequal historical/future lengths.
The result audits verify all 115 historical and 85 future years, reference
climatology, manifest baseline, event counts, rates and standard-deviation changes.

## 中文说明

本轮已实际完成 **ERA5 1979–2020 年全球 0.25°、850 hPa u/v 风场**的流函数与
JJA EOF 计算，并使用 NCEP/NCAR 做独立资料对照。EOF 区域为 **10–40°N、
110–180°E**，主结果逐格点对 JJA 流函数距平线性去趋势，第一主成分标准化后为
WPSH 指数。ERA5 的 EOF1 方差解释率为 **77.86%**，强事件年份为
**1995、1998、2010、2020**，与论文列出的年份一致。

提供可直接使用的逐年 CSV、NetCDF、PNG/SVG、数据校验值、敏感性对照及验证代码。
ACCESS-CM2 使用公开 NCI 原始 85000 Pa 风场，补充 r1、r2、r5 在 SSP245、SSP585
下的六组 1900–2099 年结果。高地下方缺测采用明确记录的球面调和延拓，并提供最近邻
延拓对照和 ERA5 人工遮挡检验。补出的风值属于数值假设，不能解释为恢复了真实风场。
逐月二次去趋势后计算 JJA EOF；历史参考期统一为 1900–2014 年（115 年），未来比较期为
2015–2099 年（85 年）。气候平均态、参考期标准化和强事件阈值均按新的历史期重算，
同时保留全时段标准化结果。强事件图比较夏季发生率（次数÷时段年数），表中另保留实际次数。
六组结果仅代表一个模式的三个成员和两个情景，不代表多模式结论。
