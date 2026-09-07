"""Shared East-China summer rainfall machinery for Tasks 2-4.

Recipe (Task 2.docx, same as the old MIROC6 run):
  - region 15N-45N, 100E-130E
  - monthly total (mm/month) / days in month = monthly rate (mm/day);
    JJA seasonal value = mean of the three monthly rates per year
  - anomaly = seasonal value minus whole-period seasonal mean
  - composite = mean over PDO+ years minus mean over PDO- years
Model precipitation is interpolated onto the GPCC 0.5-degree grid (the
"ec" grid of the old outputs) so that all members and the observations
share one grid. Maps are expressed in mm/day.
"""
from pathlib import Path

import numpy as np
import xarray as xr
from scipy import stats

REPO = Path(__file__).resolve().parents[1]

LAT_MIN, LAT_MAX = 15, 45
LON_MIN, LON_MAX = 100, 130
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


def land_mask(target_grid=None, path=None):
    """Boolean land mask (True = land) on the GPCC 0.5-degree grid.

    GPCC is a gauge-based land-only dataset: ocean cells are NaN in every
    month.  Cells with at least one finite monthly value are land.
    """
    path = path or REPO / "data/precip.mon.total.v2018.nc"
    tlat, tlon = target_grid if target_grid is not None else gpcc_grid(path)
    ds = xr.open_dataset(path)
    da = ds["precip"].isel(time=slice(0, 12))
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.isel(lat=slice(None, None, -1))
    da = da.sel(lat=tlat, lon=tlon).load()
    ds.close()
    mask = np.isfinite(da.values).any(axis=0)
    return xr.DataArray(mask, coords={"lat": tlat, "lon": tlon},
                        dims=("lat", "lon"))


def load_east_china_monthly_mm(path, var, period, target_grid=None):
    """Monthly precipitation totals (mm/month) over East China.

    target_grid=(lat, lon): bilinearly interpolate onto that grid
    (used for the model members -> GPCC 0.5-degree grid) and mask out
    ocean cells so only land precipitation is analysed/plotted.
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
        da = da.where(land_mask(target_grid))
    else:
        da = da.sel(lat=slice(LAT_MIN, LAT_MAX), lon=slice(LON_MIN, LON_MAX))
    return da


def jja_anomaly_mmday(da_mm):
    """Monthly mm/month -> mm/day rates; yearly JJA mean -> anomaly.

    Same as the old MIROC6 notebooks: divide each monthly total by its
    number of days, average the three JJA rates per year, then subtract
    the whole-period mean of that seasonal value.
    """
    jja = da_mm.sel(time=da_mm.time.dt.month.isin([6, 7, 8]))
    rate = jja / jja.time.dt.days_in_month
    seas = rate.groupby("time.year").mean("time")
    anom = seas - seas.mean("year")
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
