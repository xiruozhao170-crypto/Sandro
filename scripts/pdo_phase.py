"""Shared PDO machinery for Tasks 2-4.

Follows the course documents exactly (same recipe as the old MIROC6 run):
  - PDO = EOF1 of monthly SST anomalies over the North Pacific
    (20N-60N, 110E-250E), sqrt(cos(lat)) weights, per-gridpoint LINEAR
    detrend, monthly climatology removed (Task 1.docx).
  - PDO index = normalized PC1; sign fixed so that the positive phase has
    cool SST in the central/western North Pacific core.
  - 9-yr (108-month) centred running mean of the index; each year is
    classified positive/negative from the JJA mean of the smoothed index
    (Task 2.docx).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from eofs.standard import Eof

REPO = Path(__file__).resolve().parents[1]

JJA = (6, 7, 8)
SMOOTH_MONTHS = 108  # 9 years

# Task 1 conventions
HIST_PERIOD = slice("1900-01-01", "2014-12-31")
LAT_MIN, LAT_MAX = 20, 60
LON_MIN, LON_MAX = 110, 250


def load_north_pacific(path, var, period=HIST_PERIOD):
    ds = xr.open_dataset(path)
    da = ds[var]
    for cand, std in [("latitude", "lat"), ("longitude", "lon")]:
        if cand in da.dims:
            da = da.rename({cand: std})
    da = da.sel(time=period)
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.sel(lat=slice(LAT_MAX, LAT_MIN))
    else:
        da = da.sel(lat=slice(LAT_MIN, LAT_MAX))
    return da.sel(lon=slice(LON_MIN, LON_MAX))


def compute_pdo_eof1(da):
    """Linear detrend -> anomalies -> EOF1 with sqrt(cos(lat)) weights.

    Returns (anom, eof1[lat,lon] as covariance map, pc1[t] normalized,
    varfrac[3]).
    """
    fit = da.polyfit(dim="time", deg=1)
    trend = xr.polyval(da.time, fit.polyfit_coefficients)
    da_dt = da - trend

    clim = da_dt.groupby("time.month").mean("time")
    anom = da_dt.groupby("time.month") - clim

    wgts = np.sqrt(np.cos(np.deg2rad(anom.lat.values)))[:, np.newaxis]
    solver = Eof(anom.values.astype("float64"), weights=wgts)
    eof1 = solver.eofsAsCovariance(neofs=1, pcscaling=1)[0]
    pc1 = solver.pcs(npcs=1, pcscaling=1)[:, 0]
    varfrac = solver.varianceFraction(neigs=3)

    # sign convention: positive PDO = cool central/western North Pacific
    la, lo = anom.lat.values, anom.lon.values
    core = np.nanmean(eof1[np.ix_((la >= 35) & (la <= 45),
                                  (lo >= 160) & (lo <= 200))])
    if core > 0:
        eof1, pc1 = -eof1, -pc1
    return anom, eof1, pc1, varfrac


def pattern_corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[m], b[m])[0, 1]


def pdo_index_and_phases(pc, time):
    """108-month centred running mean; classify each year by its JJA mean.

    Returns (monthly Series raw, monthly Series smoothed,
             yearly Series of smoothed JJA mean, pos_years, neg_years).
    """
    s = pd.Series(pc, index=pd.DatetimeIndex(time))
    sm = s.rolling(SMOOTH_MONTHS, center=True, min_periods=SMOOTH_MONTHS).mean()
    jja = sm[sm.index.month.isin(JJA)]
    yearly = jja.groupby(jja.index.year).mean().dropna()
    pos_years = yearly.index[yearly.values > 0].to_numpy()
    neg_years = yearly.index[yearly.values <= 0].to_numpy()
    return s, sm, yearly, pos_years, neg_years


def get_pdo(tos_path, var, period=HIST_PERIOD, ref_eof1=None):
    """Full pipeline for one dataset; the PDO is always EOF1.

    `ref_eof1` (optional, same grid) only adds a diagnostic pattern
    correlation to the output - it never changes the chosen mode.
    """
    da = load_north_pacific(tos_path, var, period).load()
    anom, eof1, pc1, varfrac = compute_pdo_eof1(da)
    corr = (pattern_corr(eof1, ref_eof1) if ref_eof1 is not None else np.nan)
    raw, sm, yearly, pos_years, neg_years = pdo_index_and_phases(
        pc1, anom.time.values)
    return {
        "lat": anom.lat.values, "lon": anom.lon.values, "time": anom.time.values,
        "eof": eof1, "pc": pc1, "varfrac": varfrac, "corr_ref": corr,
        "pc_raw": raw, "pc_smooth": sm, "yearly": yearly,
        "pos_years": pos_years, "neg_years": neg_years,
    }
