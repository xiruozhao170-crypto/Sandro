"""Task 4: PDO-rainfall composites under future scenarios (ACCESS-CM2).

Uses the 3 members that best represent the historical rainfall composite
(results/task3_summary.json) and repeats the Task 2/3 procedure for
ssp245 and ssp585 (2015-2100). SST is detrended with a quadratic fit
(forced warming is non-linear there); the PDO mode is selected against the
observed historical PDO pattern, as in Task 2.

Outputs:
  results/task4_<scen>_<m>.nc, results/task4_summary.json
  figures/task4_<scen>_pos.png / _neg.png    (3 members + ensemble mean)
  figures/task4_<scen>_pdo_timeseries.png
  figures/task4_scenario_comparison.png      (hist vs ssp245 vs ssp585)
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from pdo_phase import get_pdo
from easm_rain import load_east_china_monthly_mm, jja_anomaly_mmday, composite
from easm_plots import _map_axis, _sym_levels, pdo_timeseries_panels
from task2_composite import save_dataset

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
RES = REPO / "results"
FIG = REPO / "figures"

SCENARIOS = ["ssp245", "ssp585"]
FUT_PERIOD = slice("2015-01-01", "2100-12-31")


def phase_maps_figure(members, comps, ens, phase, scen, path):
    """3 member panels + ensemble mean for one phase's rainfall anomaly."""
    sign = "Positive" if phase == "pos" else "Negative"
    fields = [comps[m][phase] for m in members] + [ens[phase]]
    titles = [f"({chr(97+i)})  {m}" for i, m in enumerate(members)] + \
             [f"({chr(97+len(members))})  ensemble mean"]
    lev = _sym_levels(*[np.asarray(f) for f in fields])
    la = ens[phase].lat if hasattr(ens[phase], "lat") else ens["lat"]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5), constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for ax, fld, ttl in zip(axes.flat, fields, titles):
        cf = ax.contourf(fld.lon, fld.lat, fld, levels=lev, cmap="BrBG",
                         transform=ccrs.PlateCarree(), extend="both")
        _map_axis(ax)
        ax.set_title(ttl)
        fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.8, pad=0.03)
    fig.suptitle(f"Task 4 — {scen}: JJA rainfall anomaly in PDO {sign} phase "
                 "(2015–2100, 9-yr running mean)", fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def comparison_figure(rows, path):
    """rows = [(label, {'pos','neg','diff'} ens-mean maps), ...] -> 3x3."""
    fig, axes = plt.subplots(len(rows), 3, figsize=(19, 5.2 * len(rows)),
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    lev_ab = _sym_levels(*[np.asarray(r[1][k]) for r in rows for k in ("pos", "neg")])
    lev_c = _sym_levels(*[np.asarray(r[1]["diff"]) for r in rows])
    for i, (label, ens) in enumerate(rows):
        for j, (key, ttl, lev) in enumerate([
                ("pos", "PDO+ anomaly", lev_ab),
                ("neg", "PDO$-$ anomaly", lev_ab),
                ("diff", "Composite (Pos$-$Neg)", lev_c)]):
            ax = axes[i, j]
            fld = ens[key]
            cf = ax.contourf(fld.lon, fld.lat, fld, levels=lev, cmap="BrBG",
                             transform=ccrs.PlateCarree(), extend="both")
            _map_axis(ax)
            ax.set_title(f"({chr(97+i*3+j)})  {label} — {ttl}", fontsize=11)
            fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.8, pad=0.03)
    fig.suptitle("Task 4 — PDO-phase rainfall composites: historical vs "
                 "future scenarios (ensemble mean of best-3 members)",
                 fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def ens_mean(comps, members):
    out = {}
    for key in ("pos", "neg", "diff"):
        stack = xr.concat([comps[m][key] for m in members], dim="member")
        out[key] = stack.mean("member")
    return out


def main():
    best3 = json.loads((RES / "task3_summary.json").read_text())["best3"]
    obs_eof = xr.open_dataset(RES / "task2_obs.nc")["pdo_eof"].values
    print("members:", best3)

    summary = {"members": best3, "scenarios": {}}
    rows = []

    # historical baseline: ensemble of the same 3 members (from Task 2)
    hist_comps = {}
    for m in best3:
        ds = xr.open_dataset(RES / f"task2_{m}.nc")
        hist_comps[m] = {k: ds[f"rain_{k}"].rename({"rlat": "lat", "rlon": "lon"})
                         for k in ("pos", "neg", "diff")}
    rows.append(("historical 1900–2014", ens_mean(hist_comps, best3)))

    for scen in SCENARIOS:
        comps, entries, scen_sum = {}, [], {}
        for i, m in enumerate(best3):
            pdo = get_pdo(
                DATA / f"processed/tos_2deg_{scen}_{m}i1p1f1.nc", "tos",
                ref_eof1=obs_eof, period=FUT_PERIOD, detrend_deg=2)
            rain_mm = load_east_china_monthly_mm(
                DATA / f"ACCESS-CM2_data/pr_Amon_ACCESS-CM2_{scen}_{m}i1p1f1_gn_201501-210012.nc",
                "pr", FUT_PERIOD)
            comp = composite(jja_anomaly_mmday(rain_mm),
                             pdo["pos_years"], pdo["neg_years"])
            save_dataset(f"{scen}_{m}", pdo, comp, RES / f"task4_{scen}_{m}.nc")
            comps[m] = {k: comp[k] for k in ("pos", "neg", "diff")}
            entries.append((f"({chr(97+i)})  {scen} {m} PDO "
                            f"(EOF{pdo['mode']+1})", pdo))
            scen_sum[m] = {"pdo_mode": int(pdo["mode"]) + 1,
                           "pdo_corr_with_obs": round(float(pdo["corr_ref"]), 3),
                           "n_pos": comp["n_pos"], "n_neg": comp["n_neg"]}
            print(f"[{scen} {m}] mode=EOF{pdo['mode']+1} "
                  f"corr_obs={pdo['corr_ref']:.2f} "
                  f"pos={comp['n_pos']}yr neg={comp['n_neg']}yr")
        ens = ens_mean(comps, best3)
        rows.append((scen + " 2015–2100", ens))
        phase_maps_figure(best3, comps, ens, "pos", scen,
                          FIG / f"task4_{scen}_pos.png")
        phase_maps_figure(best3, comps, ens, "neg", scen,
                          FIG / f"task4_{scen}_neg.png")
        pdo_timeseries_panels(entries, FIG / f"task4_{scen}_pdo_timeseries.png")
        summary["scenarios"][scen] = scen_sum

    comparison_figure(rows, FIG / "task4_scenario_comparison.png")

    # scenario-vs-historical pattern correlations of the ensemble composites
    from pdo_phase import pattern_corr
    hist_diff = rows[0][1]["diff"]
    for (label, ens) in rows[1:]:
        summary["scenarios"][label.split()[0]]["ens_diff_corr_vs_hist"] = round(
            float(pattern_corr(np.asarray(ens["diff"]),
                               np.asarray(hist_diff))), 3)

    (RES / "task4_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
