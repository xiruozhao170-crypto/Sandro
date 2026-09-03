"""Task 3: spread of the PDO pattern and the rainfall composites across the
ACCESS-CM2 historical members (r1-r5; the course doc assumes 10 MIROC6 runs,
the new dataset supplies 5).

Reads the per-member results of Task 2 (results/task2_<m>.nc) and computes,
grid-point by grid-point across members, the ensemble mean and variance of
the PDO EOF1 pattern and of the rainfall composites (diff and each phase).
Members are ranked by pattern correlation between their composite difference
and the observed one (same GPCC 0.5-deg grid) -> best 3 go to Task 4.

Figures (mirroring the old MIROC6 outputs):
  Task3_Figure1_PDO_EOF1_Spread.png        r1-r5 EOF1 + across-member variance
  Task3_Figure2_PC1.png                    r1-r5 PC1 + 9-yr running mean
  Task3_Figure3_Rainfall_Composite_Spread.png  r1-r5 composite + variance
  Task3_Figure4_Positive_Negative_Rainfall_Spread.png  variance per phase
Outputs: results/task3_spread.nc, results/task3_summary.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from pdo_phase import pattern_corr
from easm_plots import _map_axis, _sym_levels

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"
FIG = REPO / "figures"

MEMBERS = ["r1", "r2", "r3", "r4", "r5"]


def np_map_axis(ax):
    ax.coastlines(lw=0.6)
    ax.add_feature(cfeature.LAND, facecolor="0.8")


def stack_members():
    pdo_eofs, pcs, rain = [], [], {"pos": [], "neg": [], "diff": []}
    counts, coords = {}, {}
    for m in MEMBERS:
        ds = xr.open_dataset(RES / f"task2_{m}.nc")
        pdo_eofs.append(ds["pdo_eof"].values)
        pcs.append((ds["pdo_pc"].values, ds["pdo_pc_smooth"].values,
                    ds["time"].values))
        for k in rain:
            rain[k].append(ds[f"rain_{k}"].values)
        counts[m] = (int(ds.attrs["n_pos"]), int(ds.attrs["n_neg"]))
        coords = {"lat": ds.lat.values, "lon": ds.lon.values,
                  "rlat": ds.rlat.values, "rlon": ds.rlon.values}
    return (np.array(pdo_eofs), pcs,
            {k: np.array(v) for k, v in rain.items()}, counts, coords)


def figure1_eof_spread(eofs, coords, path):
    la, lo = coords["lat"], coords["lon"]
    lev = _sym_levels(*[e for e in eofs])
    var = eofs.var(0, ddof=1)
    fig, axes = plt.subplots(
        2, 3, figsize=(16.5, 7.4), constrained_layout=True,
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)})
    for i, (m, ax) in enumerate(zip(MEMBERS, axes.flat)):
        cf = ax.contourf(lo, la, eofs[i], levels=lev, cmap="RdBu_r",
                         transform=ccrs.PlateCarree(), extend="both")
        np_map_axis(ax)
        ax.set_title(f"{m} EOF1")
    ax = axes.flat[-1]
    cf2 = ax.contourf(lo, la, var, levels=15, cmap="YlOrRd",
                      transform=ccrs.PlateCarree(), extend="max")
    np_map_axis(ax)
    ax.set_title("Spread (Variance)")
    fig.colorbar(cf2, ax=axes[:, :].ravel().tolist(), orientation="horizontal",
                 fraction=0.05, pad=0.04, aspect=45,
                 label="Variance of PDO EOF1 pattern (°C$^2$)")
    fig.colorbar(cf, ax=axes[:, :].ravel().tolist(), orientation="horizontal",
                 fraction=0.05, pad=0.02, aspect=45,
                 label="Regression/covariance pattern (°C per normalized PC1)")
    fig.suptitle("Task 3 Figure 1: PDO EOF1 Covariance Patterns (r1–r5) + "
                 "Spread\nACCESS-CM2 Historical 1900–2014", fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def figure2_pc1(pcs, counts, path):
    fig, axes = plt.subplots(3, 2, figsize=(15, 10.5), sharex=True,
                             sharey=True)
    for i, (m, ax) in enumerate(zip(MEMBERS, axes.flat)):
        pc, sm, t = pcs[i]
        t = pd.DatetimeIndex(t)
        ax.fill_between(t, pc, 0, where=pc > 0, color="#e15759",
                        alpha=0.75, lw=0)
        ax.fill_between(t, pc, 0, where=pc <= 0, color="#4e79a7",
                        alpha=0.75, lw=0)
        ax.plot(t, sm, "k", lw=1.8,
                label="9-yr running mean" if i == 0 else None)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylim(-4, 4)
        ax.set_ylabel("PC1")
        npos, nneg = counts[m]
        ax.set_title(f"{m} PC1   (pos={npos}yr  neg={nneg}yr)")
        if i == 0:
            ax.legend(frameon=False, loc="upper right", fontsize=9)
        ax.grid(lw=0.3, alpha=0.5)
    axes.flat[-1].axis("off")
    for ax in axes[-1, :]:
        ax.set_xlabel("Year")
    fig.suptitle("Task 3 Figure 2: PDO PC1 and Smoothed PC1 (r1–r5)\n"
                 "ACCESS-CM2 Historical 1900–2014", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def figure3_composite_spread(rain, coords, corrs, path):
    la, lo = coords["rlat"], coords["rlon"]
    lev = _sym_levels(*[rain["diff"][i] for i in range(len(MEMBERS))])
    var = np.nanvar(rain["diff"], axis=0, ddof=1)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10),
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for i, (m, ax) in enumerate(zip(MEMBERS, axes.flat)):
        cf = ax.contourf(lo, la, rain["diff"][i], levels=lev, cmap="BrBG",
                         transform=ccrs.PlateCarree(), extend="both")
        _map_axis(ax)
        ax.set_title(f"{m} Composite (Pos$-$Neg)   r(obs)={corrs[m]:.2f}")
    ax = axes.flat[-1]
    cf2 = ax.contourf(lo, la, var, levels=15, cmap="YlOrRd",
                      transform=ccrs.PlateCarree(), extend="max")
    _map_axis(ax)
    ax.set_title("Spread (Variance)")
    fig.colorbar(cf, ax=axes[0, :].tolist(), label="mm/day", shrink=0.85,
                 pad=0.02)
    fig.colorbar(cf2, ax=axes[1, :].tolist(), label="(mm/day)$^2$",
                 shrink=0.85, pad=0.02)
    fig.suptitle("Task 3 Figure 3: East China JJA Rainfall Composite "
                 "(Pos$-$Neg, r1–r5) + Spread\nACCESS-CM2 Historical "
                 "1900–2014  |  r(obs) = pattern corr vs observed composite",
                 fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def figure4_phase_spread(rain, coords, path):
    la, lo = coords["rlat"], coords["rlon"]
    var_pos = np.nanvar(rain["pos"], axis=0, ddof=1)
    var_neg = np.nanvar(rain["neg"], axis=0, ddof=1)
    lev = np.linspace(0, max(np.nanmax(var_pos), np.nanmax(var_neg)), 16)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6),
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for ax, fld, ttl in zip(axes, [var_pos, var_neg],
                            ["(a)  Variance of PDO$+$ rainfall anomaly",
                             "(b)  Variance of PDO$-$ rainfall anomaly"]):
        cf = ax.contourf(lo, la, fld, levels=lev, cmap="YlOrRd",
                         transform=ccrs.PlateCarree(), extend="max")
        _map_axis(ax)
        ax.set_title(ttl)
        fig.colorbar(cf, ax=ax, label="(mm/day)$^2$", shrink=0.8, pad=0.03)
    fig.suptitle("Task 3 Figure 4: Across-member Spread of the Rainfall "
                 "Anomaly per PDO Phase (ACCESS-CM2 r1–r5, 1900–2014)",
                 fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    pdo_eofs, pcs, rain, counts, coords = stack_members()

    # obs composite (already on the same GPCC 0.5-deg grid) for ranking
    obs = xr.open_dataset(RES / "task2_obs.nc")
    obs_diff = obs["rain_diff"].values
    corrs = {m: float(pattern_corr(rain["diff"][i], obs_diff))
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
                              "historical r1-r5, 1900-2014, PDO = EOF1"},
    )
    out.to_netcdf(RES / "task3_spread.nc")

    figure1_eof_spread(pdo_eofs, coords, FIG / "Task3_Figure1_PDO_EOF1_Spread.png")
    figure2_pc1(pcs, counts, FIG / "Task3_Figure2_PC1.png")
    figure3_composite_spread(rain, coords, corrs,
                             FIG / "Task3_Figure3_Rainfall_Composite_Spread.png")
    figure4_phase_spread(rain, coords,
                         FIG / "Task3_Figure4_Positive_Negative_Rainfall_Spread.png")

    summary = {"pattern_corr_vs_obs_composite": corrs, "best3": best3}
    (RES / "task3_summary.json").write_text(json.dumps(summary, indent=2))
    print("pattern corr vs obs composite:", corrs)
    print("best 3 members for Task 4:", best3)


if __name__ == "__main__":
    main()
