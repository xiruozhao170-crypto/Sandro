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

## ACCESS-CM2 status and interface

The public NCI catalogue contains matching r1, r2 and r5 `ua`/`va` data for
historical, SSP245 and SSP585. The separate
[download script](../scripts/fetch_access_wpsh_winds.py) reads original
OPeNDAP arrays at **85000 Pa**, preserving below-ground missing values.
A one-year r1 historical wind subset was checked against its original pressure
coordinate. The initial NCSS trial exposed a service issue: it relabelled the
85000-Pa coordinate as the mean of pressure bounds, 83125 Pa. That route is
not used in the submitted ACCESS downloader.

The original model has masked 850-hPa winds below high terrain. Approximately
10% of points in the tested monthly global grid are missing (a grid-point
fraction, not an area fraction). A global streamfunction is not determined by
simply ignoring these points. The calculation therefore rejects missing winds.
**No ACCESS-CM2 WPSH numerical results are submitted in this delivery.**
They require a documented below-ground extrapolation/imputation decision or
a supplied complete global wind product. No missing winds were replaced by zero.

The prepared model-analysis configuration joins historical and each scenario
for 1900–2099, uses monthly quadratic detrending over the combined period,
interpolates the resulting streamfunction to a 1° grid, and fits a common
EOF across both centuries. The default present reference is 1900–1999, as
requested for a present-period baseline. The reference paper uses a full-record
monthly baseline for its models; the configured present-period baseline is
therefore an explicit adaptation. It must not be described as an exact
multimodel reproduction. Future WPSH indices retain the present-period SD
for event comparisons, and the full-period-normalized PC is saved separately.

## Reproduction

From the repository root, in an isolated Python environment:

```bash
python -m pip install -r requirements-wpsh.txt
python scripts/fetch_wpsh_winds.py era5
python scripts/fetch_wpsh_winds.py ncep
python scripts/task7_wpsh_eof.py --source era5
python scripts/task7_wpsh_eof.py --source ncep
python -m pytest tests/test_wpsh_streamfunction.py -q
python scripts/validate_task7_wpsh.py --raw-check
```

ERA5 download requires the user's configured CDS API access and accepted
dataset terms. Credentials are read by `cdsapi`; they are not included in
the repository. NOAA downloads need no account. Raw files and computational
caches are in git-ignored `data/wpsh_winds/`.

The figures and numerical audits can also be reproduced without raw downloads:

```bash
python scripts/task7_wpsh_eof.py --redraw
python scripts/validate_task7_wpsh.py
```

Eight scientific tests cover the analytic sign and magnitude of streamfunction,
removal of divergent winds, global coverage/missing-data rejection, pressure
selection in Pa/hPa, missing and duplicate months, equivalent linear/quadratic
temporal processing, and EOF recovery against an independent SVD. The result
audits additionally verify PC normalization/orthogonality, CSV consistency,
physical regression maps and the independent raw ERA5 calculation.

## 中文说明

本轮已实际完成 **ERA5 1979–2020 年全球 0.25°、850 hPa u/v 风场**的流函数与
JJA EOF 计算，并使用 NCEP/NCAR 做独立资料对照。EOF 区域为 **10–40°N、
110–180°E**，主结果逐格点对 JJA 流函数距平线性去趋势，第一主成分标准化后为
WPSH 指数。ERA5 的 EOF1 方差解释率为 **77.86%**，强事件年份为
**1995、1998、2010、2020**，与论文列出的年份一致。

提供可直接使用的逐年 CSV、NetCDF、PNG/SVG、数据校验值、敏感性对照及验证代码。
Drive 中现有 ACCESS 文件没有风场；公开模式风场在高地以下有缺测，补值方案会
影响全球流函数，因此本轮模型部分保留下载与计算接口，未发布经过未说明补值的模式结果。
