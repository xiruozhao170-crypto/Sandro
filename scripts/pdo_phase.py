"""Shared PDO machinery for Tasks 2-4.

Extends task1_pdo with:
  - PDO-mode selection: the EOF (of the first 3) whose pattern best matches a
    reference PDO pattern is taken as the PDO mode. Guards against mode
    swapping (e.g. historical r4, where the PDO appears as EOF2).
  - 9-yr (108-month) centred running mean of the PDO index and JJA-based
    classification of years into positive / negative phase (Task 2.docx).
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


def compute_pdo_modes(da, neofs=3, detrend_deg=1):
    """Detrend -> anomalies -> first `neofs` EOFs with sqrt(cos(lat)) weights.

    detrend_deg=2 for future scenarios, where the warming is non-linear.
    Returns (anom, eofs[k,lat,lon], pcs[t,k], varfrac[k]).
    """
    fit = da.polyfit(dim="time", deg=detrend_deg)
    trend = xr.polyval(da.time, fit.polyfit_coefficients)
    da_dt = da - trend

    clim = da_dt.groupby("time.month").mean("time")
    anom = da_dt.groupby("time.month") - clim

    wgts = np.sqrt(np.cos(np.deg2rad(anom.lat.values)))[:, np.newaxis]
    solver = Eof(anom.values.astype("float64"), weights=wgts)
    eofs = solver.eofsAsCovariance(neofs=neofs, pcscaling=1)
    pcs = solver.pcs(npcs=neofs, pcscaling=1)
    varfrac = solver.varianceFraction(neigs=neofs)
    return anom, eofs, pcs, varfrac


def pattern_corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[m], b[m])[0, 1]


def select_pdo_mode(eofs, pcs, ref_eof1):
    """Pick the EOF best matching `ref_eof1` (same grid); align its sign.

    Returns (mode_index, eof, pc, corr_with_ref).
    """
    corrs = [pattern_corr(eofs[k], ref_eof1) for k in range(eofs.shape[0])]
    k = int(np.argmax(np.abs(corrs)))
    sgn = 1.0 if corrs[k] >= 0 else -1.0
    return k, sgn * eofs[k], sgn * pcs[:, k], corrs[k] * sgn


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


def get_pdo(tos_path, var, ref_eof1=None, period=HIST_PERIOD, detrend_deg=1):
    """Full pipeline for one dataset. If `ref_eof1` is None, EOF1 is the PDO.

    Returns dict with anom coords, eof (pattern), pc (index), varfrac,
    mode (0-based), corr_ref, and phase info.
    """
    da = load_north_pacific(tos_path, var, period).load()
    anom, eofs, pcs, varfrac = compute_pdo_modes(da, detrend_deg=detrend_deg)
    if ref_eof1 is None:
        # sign convention: positive PDO = cool central/western North Pacific
        la, lo = anom.lat.values, anom.lon.values
        core = np.nanmean(eofs[0][np.ix_((la >= 35) & (la <= 45),
                                         (lo >= 160) & (lo <= 200))])
        sgn = -1.0 if core > 0 else 1.0
        k, eof, pc, corr = 0, sgn * eofs[0], sgn * pcs[:, 0], np.nan
    else:
        k, eof, pc, corr = select_pdo_mode(eofs, pcs, ref_eof1)
    raw, sm, yearly, pos_years, neg_years = pdo_index_and_phases(pc, anom.time.values)
    return {
        "lat": anom.lat.values, "lon": anom.lon.values, "time": anom.time.values,
        "eof": eof, "pc": pc, "varfrac": varfrac, "mode": k, "corr_ref": corr,
        "pc_raw": raw, "pc_smooth": sm, "yearly": yearly,
        "pos_years": pos_years, "neg_years": neg_years,
    }
