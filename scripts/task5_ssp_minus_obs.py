"""Task 5: future-scenario rainfall minus the historical observations.

For each SSP scenario (ssp245/ssp585, 2015-2100, all five ACCESS-CM2
members) the JJA-mean rainfall climatology (mm/day) is interpolated onto
the GPCC 0.5-degree grid (land-only, same mask as Tasks 2-4) and the GPCC
observed JJA climatology of the historical analysis period (1900-2014) is
subtracted.  Positive values = the scenario is wetter than the observed
historical climatology.

Figures:
  Task5_Figure1_SSP_minus_Obs_MemberMean.png   obs | ssp245-obs | ssp585-obs
  Task5_Figure2_SSP_minus_Obs_by_member.png    2x5 per-member differences
Outputs: results/task5_ssp_minus_obs.nc, results/task5_summary.json
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from pdo_phase import HIST_PERIOD
from easm_rain import gpcc_grid, load_east_china_monthly_mm
from easm_plots import _map_axis, _sym_levels

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
RES = REPO / "results"
FIG = REPO / "figures"

MEMBERS = ["r1", "r2", "r3", "r4", "r5"]
SCENARIOS = ["ssp245", "ssp585"]
FUT_PERIOD = slice("2015-01-01", "2100-12-31")


def jja_clim_mmday(da_mm):
    """Monthly totals (mm/month) -> long-term mean JJA rate (mm/day)."""
    jja = da_mm.sel(time=da_mm.time.dt.month.isin([6, 7, 8]))
    rate = jja / jja.time.dt.days_in_month
    return rate.mean("time")


def main():
    grid = gpcc_grid()

    obs_mm = load_east_china_monthly_mm(
        DATA / "precip.mon.total.v2018.nc", "precip", HIST_PERIOD)
    obs_clim = jja_clim_mmday(obs_mm)
    # GPCC native subset == target grid, but align explicitly to be safe
    obs_clim = obs_clim.sel(lat=grid[0], lon=grid[1])

    diffs = {}   # (scen, m) -> DataArray
    for scen in SCENARIOS:
        for m in MEMBERS:
            fut_mm = load_east_china_monthly_mm(
                DATA / f"ACCESS-CM2_data/pr_Amon_ACCESS-CM2_{scen}_{m}i1p1f1_gn_201501-210012.nc",
                "pr", FUT_PERIOD, target_grid=grid)
            diffs[(scen, m)] = jja_clim_mmday(fut_mm) - obs_clim
            print(f"[{scen} {m}] mean diff = "
                  f"{float(diffs[(scen, m)].mean()):+.2f} mm/day")

    ens = {scen: xr.concat([diffs[(scen, m)] for m in MEMBERS],
                           dim="member").mean("member")
           for scen in SCENARIOS}

    lat, lon = grid

    # Figure 1: obs climatology + member-mean differences
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6),
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    cf0 = axes[0].contourf(lon, lat, obs_clim, levels=np.arange(0, 13, 1),
                           cmap="YlGnBu", transform=ccrs.PlateCarree(),
                           extend="max")
    _map_axis(axes[0])
    axes[0].set_title("(a)  GPCC observed JJA rainfall\nclimatology 1900–2014")
    fig.colorbar(cf0, ax=axes[0], label="mm/day", shrink=0.8, pad=0.03)
    lev = _sym_levels(ens["ssp245"].values, ens["ssp585"].values)
    for ax, scen, tag in zip(axes[1:], SCENARIOS, ["(b)", "(c)"]):
        cf = ax.contourf(lon, lat, ens[scen], levels=lev, cmap="BrBG",
                         transform=ccrs.PlateCarree(), extend="both")
        _map_axis(ax)
        ax.set_title(f"{tag}  {scen} (2015–2100) $-$ obs\n5-member mean")
        fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.8, pad=0.03)
    fig.suptitle("Task 5: East China JJA Rainfall — SSP scenarios minus GPCC "
                 "observations (land only)\nACCESS-CM2 2015–2100 vs GPCC "
                 "1900–2014", fontsize=13)
    fig.savefig(FIG / "Task5_Figure1_SSP_minus_Obs_MemberMean.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Figure 2: per-member differences (rows = scenario, cols = member)
    lev = _sym_levels(*[diffs[c].values for c in diffs])
    fig, axes = plt.subplots(2, 5, figsize=(22, 8),
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for i, scen in enumerate(SCENARIOS):
        for j, m in enumerate(MEMBERS):
            ax = axes[i, j]
            cf = ax.contourf(lon, lat, diffs[(scen, m)], levels=lev,
                             cmap="BrBG", transform=ccrs.PlateCarree(),
                             extend="both")
            _map_axis(ax)
            ax.set_title(f"{scen} {m} $-$ obs", fontsize=11)
    fig.colorbar(cf, ax=axes.ravel().tolist(), label="mm/day",
                 orientation="horizontal", fraction=0.05, pad=0.03, aspect=50)
    fig.suptitle("Task 5: JJA Rainfall, SSP scenario (2015–2100) minus GPCC "
                 "observations (1900–2014), per member", fontsize=13)
    fig.savefig(FIG / "Task5_Figure2_SSP_minus_Obs_by_member.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    out = xr.Dataset(
        {"obs_jja_clim": (("lat", "lon"), obs_clim.values),
         **{f"diff_{scen}_{m}": (("lat", "lon"), diffs[(scen, m)].values)
            for scen in SCENARIOS for m in MEMBERS},
         **{f"diff_{scen}_ensmean": (("lat", "lon"), ens[scen].values)
            for scen in SCENARIOS}},
        coords={"lat": lat, "lon": lon},
        attrs={"description": "Task 5: ACCESS-CM2 SSP JJA rainfall "
                              "climatology (2015-2100) minus GPCC observed "
                              "JJA climatology (1900-2014), land only, "
                              "GPCC 0.5-deg grid, mm/day"})
    out.to_netcdf(RES / "task5_ssp_minus_obs.nc")

    summary = {scen: {
        "ens_mean_diff_mmday": round(float(ens[scen].mean()), 3),
        **{m: round(float(diffs[(scen, m)].mean()), 3) for m in MEMBERS}}
        for scen in SCENARIOS}
    (RES / "task5_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
