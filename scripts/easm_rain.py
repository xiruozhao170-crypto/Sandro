"""Shared East-China summer rainfall machinery for Tasks 2-4.

Recipe (Task 2.docx):
  - region 15N-45N, 100E-130E
  - JJA seasonal value = sum of the three monthly totals per year
  - anomaly = seasonal value minus whole-period seasonal mean
  - composite = mean over PDO+ years minus mean over PDO- years
Maps are expressed in mm/day (seasonal total / 92 days) for comparability.
"""
from pathlib import Path

import numpy as np
import xarray as xr
from scipy import stats

REPO = Path(__file__).resolve().parents[1]

LAT_MIN, LAT_MAX = 15, 45
LON_MIN, LON_MAX = 100, 130
JJA_DAYS = 92.0


def load_east_china_monthly_mm(path, var, period):
    """Monthly precipitation totals (mm/month) over East China."""
    ds = xr.open_dataset(path)
    da = ds[var]
    for cand, std in [("latitude", "lat"), ("longitude", "lon")]:
        if cand in da.dims:
            da = da.rename({cand: std})
    da = da.sel(time=period)
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.sel(lat=slice(LAT_MAX, LAT_MIN))
        da = da.isel(lat=slice(None, None, -1))  # ascending lat for plotting
    else:
        da = da.sel(lat=slice(LAT_MIN, LAT_MAX))
    da = da.sel(lon=slice(LON_MIN, LON_MAX)).load()

    units = da.attrs.get("units", "")
    if units in ("kg m-2 s-1", "kg/m2/s"):  # CMIP pr flux -> mm/month
        days = da.time.dt.days_in_month
        da = da * 86400.0 * days
        da.attrs["units"] = "mm/month"
    return da


def jja_anomaly_mmday(da_mm, detrend_deg=None):
    """Yearly JJA total -> anomaly vs whole-period mean, in mm/day.

    detrend_deg: if set, remove a per-gridpoint polynomial trend of that
    degree instead of the whole-period mean (Task 4 sensitivity check;
    mirrors the SST treatment: 1 for historical, 2 for scenarios).
    Returns DataArray (year, lat, lon).
    """
    jja = da_mm.sel(time=da_mm.time.dt.month.isin([6, 7, 8]))
    tot = jja.groupby("time.year").sum("time", min_count=3) / JJA_DAYS
    if detrend_deg:
        fit = tot.polyfit(dim="year", deg=detrend_deg)
        anom = tot - xr.polyval(tot.year, fit.polyfit_coefficients)
    else:
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
