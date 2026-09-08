"""Task 4: PDO-rainfall composites under future scenarios (ACCESS-CM2).

Uses the 3 prescribed members (r2, r3, r5; "best3" in
results/task3_summary.json) and repeats the Task 2/3 procedure for
ssp245 and ssp585 (2015-2100), exactly as in the historical run:
PDO = EOF1 of linearly detrended monthly SST anomalies; phases from the
JJA mean of the 9-yr smoothed PC1; rainfall on the GPCC 0.5-deg grid.

Figures (mirroring the old MIROC6 T4 outputs, adapted to 2 scenarios):
  Task4_T4-1_Future_PDO_EOF1_by_member.png       2x3 EOF1 maps
  Task4_T4-2_Future_PC1_by_member.png            2x3 PC1 + smoothed
  Task4_T4-3_Positive_Phase_Rainfall_by_member.png
  Task4_T4-4_Negative_Phase_Rainfall_by_member.png
  Task4_T4-5_Composite_Difference_3MemberMean.png  hist / ssp245 / ssp585
  Task4_T4-6_Composite_Difference_by_member.png
  Task4_T4-7_Scenario_Difference_SSP585_minus_SSP245.png
Outputs: results/task4_<scen>_<m>.nc, results/task4_summary.json
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

from pdo_phase import get_pdo, pattern_corr
from easm_rain import (gpcc_grid, load_east_china_monthly_mm,
                       jja_anomaly_mmday, composite)
from easm_plots import _map_axis, _sym_levels
from task3_spread import np_map_axis
from task2_composite import save_dataset

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
RES = REPO / "results"
FIG = REPO / "figures"

SCENARIOS = ["ssp245", "ssp585"]
FUT_PERIOD = slice("2015-01-01", "2100-12-31")


def grid_maps(fields, titles, path, suptitle, cmap, units, coords=None,
              east_china=True, levels=None):
    """2x3 grid of maps: rows = scenarios, cols = members."""
    nrow, ncol = 2, 3
    proj = (ccrs.PlateCarree() if east_china
            else ccrs.PlateCarree(central_longitude=180))
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(15, 8.6 if east_china else 6.6),
                             constrained_layout=True,
                             subplot_kw={"projection": proj})
    lev = levels if levels is not None else _sym_levels(
        *[np.asarray(f) for f, _, _ in fields])
    for ax, (fld, la, lo), ttl in zip(axes.flat, fields, titles):
        cf = ax.contourf(lo, la, fld, levels=lev, cmap=cmap,
                         transform=ccrs.PlateCarree(), extend="both")
        if east_china:
            _map_axis(ax)
        else:
            np_map_axis(ax)
        ax.set_title(ttl, fontsize=11)
    fig.colorbar(cf, ax=axes.ravel().tolist(), label=units,
                 orientation="horizontal", fraction=0.05, pad=0.03, aspect=45)
    fig.suptitle(suptitle, fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def pc_grid(pdos, titles, path, suptitle):
    fig, axes = plt.subplots(2, 3, figsize=(16, 7), sharey=True)
    for ax, pdo, ttl in zip(axes.flat, pdos, titles):
        raw, sm = pdo["pc_raw"], pdo["pc_smooth"]
        ax.fill_between(raw.index, raw.values, 0, where=raw.values > 0,
                        color="#e15759", alpha=0.75, lw=0)
        ax.fill_between(raw.index, raw.values, 0, where=raw.values <= 0,
                        color="#4e79a7", alpha=0.75, lw=0)
        ax.plot(sm.index, sm.values, "k", lw=1.6)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylim(-4, 4)
        ax.set_title(ttl, fontsize=11)
        ax.grid(lw=0.3, alpha=0.5)
    for ax in axes[:, 0]:
        ax.set_ylabel("PC1")
    for ax in axes[-1, :]:
        ax.set_xlabel("Year")
    fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    best3 = json.loads((RES / "task3_summary.json").read_text())["best3"]
    obs_eof = xr.open_dataset(RES / "task2_obs.nc")["pdo_eof"].values
    grid = gpcc_grid()
    print("members:", best3)

    # historical baseline: same 3 members from Task 2
    hist_diff = []
    for m in best3:
        ds = xr.open_dataset(RES / f"task2_{m}.nc")
        hist_diff.append(ds["rain_diff"].values)
        rl, rn = ds["rlat"].values, ds["rlon"].values
    hist_diff_mean = np.nanmean(hist_diff, axis=0)

    summary = {"members": best3, "scenarios": {}}
    results = {}  # (scen, m) -> dict(pdo=..., comp=...)

    for scen in SCENARIOS:
        scen_sum = {}
        for m in best3:
            pdo = get_pdo(
                DATA / f"processed/tos_2deg_{scen}_{m}i1p1f1.nc", "tos",
                period=FUT_PERIOD, ref_eof1=obs_eof)
            rain_mm = load_east_china_monthly_mm(
                DATA / f"ACCESS-CM2_data/pr_Amon_ACCESS-CM2_{scen}_{m}i1p1f1_gn_201501-210012.nc",
                "pr", FUT_PERIOD, target_grid=grid)
            comp = composite(jja_anomaly_mmday(rain_mm),
                             pdo["pos_years"], pdo["neg_years"])
            save_dataset(f"{scen}_{m}", pdo, comp, RES / f"task4_{scen}_{m}.nc")
            results[(scen, m)] = {"pdo": pdo, "comp": comp}
            scen_sum[m] = {
                "eof1_varfrac": round(float(pdo["varfrac"][0]), 4),
                "eof1_corr_with_obs_hist": round(float(pdo["corr_ref"]), 3),
                "n_pos": comp["n_pos"], "n_neg": comp["n_neg"]}
            print(f"[{scen} {m}] EOF1 var={pdo['varfrac'][0]*100:.1f}% "
                  f"corr_obs={pdo['corr_ref']:.2f} "
                  f"pos={comp['n_pos']}yr neg={comp['n_neg']}yr")
        summary["scenarios"][scen] = scen_sum

    cells = [(scen, m) for scen in SCENARIOS for m in best3]

    # T4-1 EOF1 maps
    grid_maps(
        [(results[c]["pdo"]["eof"], results[c]["pdo"]["lat"],
          results[c]["pdo"]["lon"]) for c in cells],
        [f"{scen} {m} EOF1 "
         f"({results[(scen, m)]['pdo']['varfrac'][0]*100:.0f}%)"
         for scen, m in cells],
        FIG / "Task4_T4-1_Future_PDO_EOF1_by_member.png",
        "Task 4-1: Future PDO EOF1 Covariance Patterns (best-3 members)\n"
        "ACCESS-CM2 2015–2100  (°C per normalized PC1)",
        "RdBu_r", "°C per normalized PC1", east_china=False)

    # T4-2 PC1
    pc_grid([results[c]["pdo"] for c in cells],
            [f"{scen} {m} PC1  (pos={results[(scen, m)]['comp']['n_pos']}yr  "
             f"neg={results[(scen, m)]['comp']['n_neg']}yr)"
             for scen, m in cells],
            FIG / "Task4_T4-2_Future_PC1_by_member.png",
            "Task 4-2: Future PDO PC1 and 9-yr Running Mean (best-3 members)\n"
            "ACCESS-CM2 2015–2100")

    # T4-3 / T4-4 phase rainfall maps
    for key, name, ttl in [("pos", "Task4_T4-3_Positive_Phase_Rainfall_by_member.png",
                            "PDO Positive Phase JJA Rainfall Anomaly"),
                           ("neg", "Task4_T4-4_Negative_Phase_Rainfall_by_member.png",
                            "PDO Negative Phase JJA Rainfall Anomaly")]:
        grid_maps(
            [(results[c]["comp"][key].values, rl, rn) for c in cells],
            [f"{scen} {m}" for scen, m in cells],
            FIG / name,
            f"Task 4: {ttl} (best-3 members)\nACCESS-CM2 2015–2100  |  "
            "9-yr running mean phases", "BrBG", "mm/day")

    # T4-6 composite difference by member
    grid_maps(
        [(results[c]["comp"]["diff"].values, rl, rn) for c in cells],
        [f"{scen} {m} (Pos$-$Neg)" for scen, m in cells],
        FIG / "Task4_T4-6_Composite_Difference_by_member.png",
        "Task 4-6: JJA Rainfall Composite Difference (Pos$-$Neg) by member\n"
        "ACCESS-CM2 2015–2100", "BrBG", "mm/day")

    # T4-5 member-mean composite difference: historical vs scenarios
    ens = {scen: np.nanmean([results[(scen, m)]["comp"]["diff"].values
                             for m in best3], axis=0) for scen in SCENARIOS}
    fields = [hist_diff_mean, ens["ssp245"], ens["ssp585"]]
    labels = ["historical 1900–2014", "ssp245 2015–2100", "ssp585 2015–2100"]
    lev = _sym_levels(*fields)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for ax, fld, ttl in zip(axes, fields, labels):
        cf = ax.contourf(rn, rl, fld, levels=lev, cmap="BrBG",
                         transform=ccrs.PlateCarree(), extend="both")
        _map_axis(ax)
        ax.set_title(ttl)
    fig.colorbar(cf, ax=axes.ravel().tolist(), label="mm/day", shrink=0.8,
                 pad=0.02)
    fig.suptitle("Task 4-5: Composite Difference (Pos$-$Neg), mean of best-3 "
                 "members\nACCESS-CM2: historical vs future scenarios",
                 fontsize=13)
    fig.savefig(FIG / "Task4_T4-5_Composite_Difference_3MemberMean.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    # T4-7 scenario difference of the member-mean composites
    dd = ens["ssp585"] - ens["ssp245"]
    lev = _sym_levels(dd)
    fig = plt.figure(figsize=(7.2, 5.8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    cf = ax.contourf(rn, rl, dd, levels=lev, cmap="BrBG",
                     transform=ccrs.PlateCarree(), extend="both")
    _map_axis(ax)
    fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.85)
    ax.set_title("Task 4-7: Composite Difference, SSP585 $-$ SSP245\n"
                 "(mean of best-3 members, 2015–2100)")
    fig.savefig(FIG / "Task4_T4-7_Scenario_Difference_SSP585_minus_SSP245.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    for scen in SCENARIOS:
        summary["scenarios"][scen]["ens_diff_corr_vs_hist"] = round(
            float(pattern_corr(ens[scen], hist_diff_mean)), 3)
    summary["ssp585_vs_ssp245_ens_diff_corr"] = round(
        float(pattern_corr(ens["ssp585"], ens["ssp245"])), 3)

    (RES / "task4_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
