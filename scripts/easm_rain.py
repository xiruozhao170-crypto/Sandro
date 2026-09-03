"""Shared East-China summer rainfall machinery for Tasks 2-4.

Recipe (Task 2.docx, same as the old MIROC6 run):
  - region 15N-45N, 100E-130E
  - JJA seasonal value = sum of the three monthly totals per year
  - anomaly = seasonal value minus whole-period seasonal mean
  - composite = mean over PDO+ years minus mean over PDO- years
Model precipitation is interpolated onto the GPCC 0.5-degree grid (the
"ec" grid of the old outputs) so that all members and the observations
share one grid. Maps are expressed in mm/day (seasonal total / 92 days).
"""
from pathlib import Path

import numpy as np
import xarray as xr
from scipy import stats

REPO = Path(__file__).resolve().parents[1]

LAT_MIN, LAT_MAX = 15, 45
LON_MIN, LON_MAX = 100, 130
JJA_DAYS = 92.0
PAD = 3.0  # extra margin kept before regridding


def gpcc_grid(path=None):
    """Target 0.5-degree lat/lon (ascending) of the GPCC East-China box."""
    path = path or REPO / "data/precip.mon.total.v2018.nc"
    ds = xr.open_dataset(path)
    lat = ds["lat"].values
    lon = ds["lon"].values
    lat = np.sort(lat[(lat >= LAT_MIN) & (lat <= LAT_MAX)])
    lon = np.sort(lon[(lon >= LON_MIN) & (lon <= LON_MAX)])
    ds.close()
    return lat, lon


def load_east_china_monthly_mm(path, var, period, target_grid=None):
    """Monthly precipitation totals (mm/month) over East China.

    target_grid=(lat, lon): bilinearly interpolate onto that grid
    (used for the model members -> GPCC 0.5-degree grid).
    """
    ds = xr.open_dataset(path)
    da = ds[var]
    for cand, std in [("latitude", "lat"), ("longitude", "lon")]:
        if cand in da.dims:
            da = da.rename({cand: std})
    da = da.sel(time=period)
    lo0, lo1 = LON_MIN - PAD, LON_MAX + PAD
    la0, la1 = LAT_MIN - PAD, LAT_MAX + PAD
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.sel(lat=slice(la1, la0))
        da = da.isel(lat=slice(None, None, -1))  # ascending lat
    else:
        da = da.sel(lat=slice(la0, la1))
    da = da.sel(lon=slice(lo0, lo1)).load()

    units = da.attrs.get("units", "")
    if units in ("kg m-2 s-1", "kg/m2/s"):  # CMIP pr flux -> mm/month
        days = da.time.dt.days_in_month
        da = da * 86400.0 * days
        da.attrs["units"] = "mm/month"

    if target_grid is not None:
        tlat, tlon = target_grid
        da = da.interp(lat=tlat, lon=tlon, method="linear")
    else:
        da = da.sel(lat=slice(LAT_MIN, LAT_MAX), lon=slice(LON_MIN, LON_MAX))
    return da


def jja_anomaly_mmday(da_mm):
    """Yearly JJA total -> anomaly vs whole-period mean, in mm/day."""
    jja = da_mm.sel(time=da_mm.time.dt.month.isin([6, 7, 8]))
    tot = jja.groupby("time.year").sum("time", min_count=3) / JJA_DAYS
    anom = tot - tot.mean("year")
    anom.attrs["units"] = "mm/day"
    return anom


def composite(anom, pos_years, neg_years):
    """Composite maps + Welch t-test p-values for pos-vs-neg difference."""
    pos_years = [y for y in pos_years if y in anom.year.values]
    neg_years = [y for y in neg_years if y in anom.year.values]
    pos = anom.sel(year=pos_years)
    neg = anom.sel(year=neg_years)
    pos_mean = pos.mean("year")
    neg_mean = neg.mean("year")
    diff = pos_mean - neg_mean
    _, p = stats.ttest_ind(pos.values, neg.values, axis=0, equal_var=False,
                           nan_policy="omit")
    pval = xr.DataArray(np.asarray(p), coords=diff.coords, dims=diff.dims)
    return {"pos": pos_mean, "neg": neg_mean, "diff": diff, "pval": pval,
            "n_pos": len(pos_years), "n_neg": len(neg_years)}
