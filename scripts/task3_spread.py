"""Task 3: spread of the PDO pattern and the rainfall composites across the
ACCESS-CM2 historical members (r1-r5; the course doc assumes 10 MIROC6 runs,
we have the 5 members supplied with the new dataset).

Reads the per-member results of Task 2 (results/task2_<m>.nc) and computes,
grid-point by grid-point across members:
  - ensemble mean and variance of the PDO pattern
  - ensemble mean and variance of the rainfall composite difference
  - variance of the rainfall composite for each phase separately

Also ranks members by pattern correlation between their composite difference
and the observed one (GPCC interpolated to the model grid) -> input to Task 4.

Outputs: results/task3_spread.nc, results/task3_summary.json,
figures/task3_pdo_spread.png, task3_rain_spread.png, task3_members_diff.png
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from pdo_phase import pattern_corr
from easm_plots import _map_axis, _sym_levels

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"
FIG = REPO / "figures"

MEMBERS = ["r1", "r2", "r3", "r4", "r5"]


def stack_members():
    pdo_eofs, rain = [], {"pos": [], "neg": [], "diff": []}
    coords = {}
    for m in MEMBERS:
        ds = xr.open_dataset(RES / f"task2_{m}.nc")
        pdo_eofs.append(ds["pdo_eof"].values)
        for k in rain:
            rain[k].append(ds[f"rain_{k}"].values)
        coords = {"lat": ds.lat.values, "lon": ds.lon.values,
                  "rlat": ds.rlat.values, "rlon": ds.rlon.values}
    return (np.array(pdo_eofs),
            {k: np.array(v) for k, v in rain.items()}, coords)


def plot_pdo_spread(mean, var, coords, path):
    fig, axes = plt.subplots(
        1, 2, figsize=(14, 4.6), constrained_layout=True,
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)})
    lev = _sym_levels(mean)
    cf = axes[0].contourf(coords["lon"], coords["lat"], mean, levels=lev,
                          cmap="RdBu_r", transform=ccrs.PlateCarree(),
                          extend="both")
    cf2 = axes[1].contourf(coords["lon"], coords["lat"], var,
                           levels=15, cmap="viridis",
                           transform=ccrs.PlateCarree(), extend="max")
    for ax, ttl in zip(axes, ["(a)  Ensemble-mean PDO pattern",
                              "(b)  Across-member variance"]):
        ax.coastlines()
        ax.set_title(ttl)
    fig.colorbar(cf, ax=axes[0], label="SST (degC per s.d.)", shrink=0.85)
    fig.colorbar(cf2, ax=axes[1], label="variance (degC$^2$)", shrink=0.85)
    fig.suptitle("Task 3 — PDO pattern spread, ACCESS-CM2 historical r1–r5 "
                 "(1900–2014)")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_rain_spread(rain, coords, path):
    la, lo = coords["rlat"], coords["rlon"]
    mean = {k: np.nanmean(v, axis=0) for k, v in rain.items()}
    var = {k: np.nanvar(v, axis=0, ddof=1) for k, v in rain.items()}
    fig, axes = plt.subplots(2, 3, figsize=(19, 10.5), constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    lev_ab = _sym_levels(mean["pos"], mean["neg"])
    lev_c = _sym_levels(mean["diff"])
    tops = [("(a)  Ens-mean PDO+ anomaly", mean["pos"], lev_ab, "BrBG"),
            ("(b)  Ens-mean PDO$-$ anomaly", mean["neg"], lev_ab, "BrBG"),
            ("(c)  Ens-mean composite (Pos$-$Neg)", mean["diff"], lev_c, "BrBG")]
    bots = [("(d)  Variance of PDO+ anomaly", var["pos"]),
            ("(e)  Variance of PDO$-$ anomaly", var["neg"]),
            ("(f)  Variance of composite", var["diff"])]
    for ax, (ttl, fld, lev, cmap) in zip(axes[0], tops):
        cf = ax.contourf(lo, la, fld, levels=lev, cmap=cmap,
                         transform=ccrs.PlateCarree(), extend="both")
        _map_axis(ax)
        ax.set_title(ttl)
        fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.8, pad=0.03)
    for ax, (ttl, fld) in zip(axes[1], bots):
        cf = ax.contourf(lo, la, fld, levels=15, cmap="viridis",
                         transform=ccrs.PlateCarree(), extend="max")
        _map_axis(ax)
        ax.set_title(ttl)
        fig.colorbar(cf, ax=ax, label="(mm/day)$^2$", shrink=0.8, pad=0.03)
    fig.suptitle("Task 3 — East China JJA rainfall composites: ensemble mean "
                 "and across-member spread (ACCESS-CM2 r1–r5, 1900–2014)",
                 fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_members_diff(rain, coords, corrs, path):
    la, lo = coords["rlat"], coords["rlon"]
    lev = _sym_levels(*[rain["diff"][i] for i in range(len(MEMBERS))])
    fig, axes = plt.subplots(2, 3, figsize=(19, 10.5), constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for i, (m, ax) in enumerate(zip(MEMBERS, axes.flat)):
        cf = ax.contourf(lo, la, rain["diff"][i], levels=lev, cmap="BrBG",
                         transform=ccrs.PlateCarree(), extend="both")
        _map_axis(ax)
        ax.set_title(f"({chr(97+i)})  {m}   r(obs)={corrs[m]:.2f}")
        fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.8, pad=0.03)
    axes.flat[-1].axis("off")
    fig.suptitle("Task 3 — composite (Pos$-$Neg) per member, with pattern "
                 "correlation vs observed composite", fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    pdo_eofs, rain, coords = stack_members()

    # obs composite interpolated to the model rainfall grid for ranking
    obs = xr.open_dataset(RES / "task2_obs.nc")
    obs_diff = (obs["rain_diff"].rename({"rlat": "lat", "rlon": "lon"})
                .interp(lat=coords["rlat"], lon=coords["rlon"]))
    corrs = {m: float(pattern_corr(rain["diff"][i], obs_diff.values))
             for i, m in enumerate(MEMBERS)}
    best3 = sorted(corrs, key=lambda m: -corrs[m])[:3]

    out = xr.Dataset(
        {
            "pdo_eof_mean": (("lat", "lon"), pdo_eofs.mean(0)),
            "pdo_eof_var": (("lat", "lon"), pdo_eofs.var(0, ddof=1)),
            "rain_pos_mean": (("rlat", "rlon"), np.nanmean(rain["pos"], 0)),
            "rain_neg_mean": (("rlat", "rlon"), np.nanmean(rain["neg"], 0)),
            "rain_diff_mean": (("rlat", "rlon"), np.nanmean(rain["diff"], 0)),
            "rain_pos_var": (("rlat", "rlon"), np.nanvar(rain["pos"], 0, ddof=1)),
            "rain_neg_var": (("rlat", "rlon"), np.nanvar(rain["neg"], 0, ddof=1)),
            "rain_diff_var": (("rlat", "rlon"), np.nanvar(rain["diff"], 0, ddof=1)),
        },
        coords=coords,
        attrs={"members": ",".join(MEMBERS),
               "description": "Task 3 across-member mean/variance, ACCESS-CM2 "
                              "historical r1-r5, 1900-2014"},
    )
    out.to_netcdf(RES / "task3_spread.nc")

    plot_pdo_spread(pdo_eofs.mean(0), pdo_eofs.var(0, ddof=1), coords,
                    FIG / "task3_pdo_spread.png")
    plot_rain_spread(rain, coords, FIG / "task3_rain_spread.png")
    plot_members_diff(rain, coords, corrs, FIG / "task3_members_diff.png")

    summary = {"pattern_corr_vs_obs_composite": corrs, "best3": best3}
    (RES / "task3_summary.json").write_text(json.dumps(summary, indent=2))
    print("pattern corr vs obs composite:", corrs)
    print("best 3 members for Task 4:", best3)


if __name__ == "__main__":
    main()
