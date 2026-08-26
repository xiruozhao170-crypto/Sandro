"""Task 1: Compute the PDO as the leading EOF of North Pacific SST anomalies.

Recipe (Task 1.docx):
  1. Monthly SST, select Jan 1900 - Dec 2014.
  2. North Pacific box 20N-60N, 110E-250E (plot one month as a sanity check).
  3. sqrt(cos(lat)) area weighting.
  4. Per-gridpoint linear detrend (least squares).
  5. Anomalies: subtract the 1900-2014 monthly climatology.
  6. EOF1 = PDO spatial pattern, PC1 = PDO index.

Usage:  python task1_pdo.py <dataset>     dataset in {ersst, access}
Outputs: results/task1_pdo_<dataset>.nc, figures/task1_*_<dataset>.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from eofs.standard import Eof
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"
FIG = REPO / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

PERIOD = slice("1900-01-01", "2014-12-31")
LAT_MIN, LAT_MAX = 20, 60
LON_MIN, LON_MAX = 110, 250

DATASETS = {
    "ersst": {
        "path": REPO / "data/ersst.v4.sst.mnmean.nc",
        "var": "sst",
        "label": "ERSSTv4 (obs)",
    },
    "access": {
        "path": REPO / "data/processed/tos_2deg_historical_r1i1p1f1.nc",
        "var": "tos",
        "label": "ACCESS-CM2 historical r1",
    },
}
for _m in ["r1", "r2", "r3", "r4", "r5"]:
    DATASETS[f"access_{_m}"] = {
        "path": REPO / f"data/processed/tos_2deg_historical_{_m}i1p1f1.nc",
        "var": "tos",
        "label": f"ACCESS-CM2 historical {_m}",
    }


def load_north_pacific(path, var):
    ds = xr.open_dataset(path)
    da = ds[var]
    # normalise coordinate names
    for cand, std in [("latitude", "lat"), ("longitude", "lon")]:
        if cand in da.dims:
            da = da.rename({cand: std})
    da = da.sel(time=PERIOD)
    # handle ascending or descending latitude axes
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.sel(lat=slice(LAT_MAX, LAT_MIN))
    else:
        da = da.sel(lat=slice(LAT_MIN, LAT_MAX))
    da = da.sel(lon=slice(LON_MIN, LON_MAX))
    return da


def sanity_map(da, label, tag):
    """Plot one arbitrary month to confirm the regional selection."""
    fig = plt.figure(figsize=(9, 4))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    da.isel(time=7).plot(
        ax=ax, transform=ccrs.PlateCarree(), cmap="RdYlBu_r",
        cbar_kwargs={"label": "SST (degC)"},
    )
    ax.coastlines()
    ax.add_feature(cfeature.LAND, facecolor="0.85")
    ax.set_title(f"{label}: SST, month #8 of record (region check)")
    fig.savefig(FIG / f"task1_regioncheck_{tag}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def compute_pdo(da):
    """Detrend -> anomalies -> EOF1 with sqrt(cos(lat)) weights."""
    # per-gridpoint linear detrend
    fit = da.polyfit(dim="time", deg=1)
    trend = xr.polyval(da.time, fit.polyfit_coefficients)
    da_dt = da - trend

    # remove the monthly climatology of the full analysis period
    clim = da_dt.groupby("time.month").mean("time")
    anom = da_dt.groupby("time.month") - clim

    wgts = np.sqrt(np.cos(np.deg2rad(anom.lat.values)))[:, np.newaxis]
    solver = Eof(anom.values.astype("float64"), weights=wgts)
    varfrac = solver.varianceFraction(neigs=3)
    eof1 = solver.eofsAsCovariance(neofs=1, pcscaling=1)[0]  # degC per sd of PC1
    pc1 = solver.pcs(npcs=1, pcscaling=1)[:, 0]

    # sign convention: positive PDO = cool central/western North Pacific
    la, lo = anom.lat.values, anom.lon.values
    core = np.nanmean(eof1[np.ix_((la >= 35) & (la <= 45), (lo >= 160) & (lo <= 200))])
    if core > 0:
        eof1, pc1 = -eof1, -pc1

    return anom, eof1, pc1, varfrac


def plot_results(anom, eof1, pc1, varfrac, label, tag):
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(2, 1, 1, projection=ccrs.PlateCarree(central_longitude=180))
    lim = np.nanmax(np.abs(eof1))
    cf = ax.contourf(
        anom.lon, anom.lat, eof1, levels=np.linspace(-lim, lim, 21),
        cmap="RdBu_r", transform=ccrs.PlateCarree(), extend="both",
    )
    ax.coastlines()
    ax.add_feature(cfeature.LAND, facecolor="0.85")
    fig.colorbar(cf, ax=ax, label="SST (degC per s.d. of PC1)", shrink=0.9)
    ax.set_title(f"{label}: PDO pattern (EOF1, {varfrac[0]*100:.1f}% of variance)")

    ax2 = fig.add_subplot(2, 1, 2)
    t = pd.DatetimeIndex(anom.time.values)
    pc_s = pd.Series(pc1, index=t)
    ann = pc_s.resample("YE").mean()
    ax2.bar(t, pc1, width=31, color=np.where(pc1 > 0, "#d6604d", "#4393c3"), alpha=0.5)
    ax2.plot(ann.index, ann.values, color="k", lw=1.5, label="annual mean")
    ax2.axhline(0, color="k", lw=0.5)
    ax2.set_ylabel("PDO index (s.d.)")
    ax2.set_title(f"{label}: PDO index (PC1), 1900-2014")
    ax2.legend(frameon=False)
    fig.savefig(FIG / f"task1_pdo_{tag}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main(tag):
    cfg = DATASETS[tag]
    if cfg["path"] is None:
        sys.exit(f"dataset '{tag}' not configured yet (waiting for data)")
    da = load_north_pacific(cfg["path"], cfg["var"])
    print(f"[{tag}] input:", dict(da.sizes),
          str(da.time.values[0])[:10], "->", str(da.time.values[-1])[:10])
    da = da.load()
    sanity_map(da, cfg["label"], tag)

    anom, eof1, pc1, varfrac = compute_pdo(da)
    print(f"[{tag}] EOF1/2/3 variance fraction: "
          + ", ".join(f"{v*100:.1f}%" for v in varfrac))
    plot_results(anom, eof1, pc1, varfrac, cfg["label"], tag)

    out = xr.Dataset(
        {
            "eof1": (("lat", "lon"), eof1),
            "pc1": (("time",), pc1),
            "variance_fraction": (("mode",), np.asarray(varfrac)),
        },
        coords={"lat": anom.lat, "lon": anom.lon, "time": anom.time,
                "mode": [1, 2, 3]},
        attrs={"description": f"Task 1 PDO ({cfg['label']}), 1900-2014, "
                              "EOF of detrended monthly SST anomalies, "
                              "20N-60N 110E-250E, sqrt(cos(lat)) weights"},
    )
    out.to_netcdf(RES / f"task1_pdo_{tag}.nc")
    print(f"[{tag}] saved results/task1_pdo_{tag}.nc, "
          f"figures/task1_pdo_{tag}.png, figures/task1_regioncheck_{tag}.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "ersst")
