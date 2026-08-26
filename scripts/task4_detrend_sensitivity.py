"""Task 4 add-on: sensitivity of the rainfall composites to detrending.

The baseline Task 2/4 composites remove only the whole-period JJA mean, so
under the SSP scenarios the forced precipitation trend can leak into the
PDO+/PDO- composites (the smoothed phases cluster in time). Here the JJA
rainfall is additionally detrended per gridpoint before compositing —
degree 1 for historical, degree 2 for the scenarios, mirroring the SST
treatment — and the result is compared with the original composites.

PDO phase years are reused from the saved Task 2 / Task 4 datasets, so the
only change is the rainfall preprocessing.

Outputs:
  results/task4_detrend_<tag>.nc         member + ens composites (detrended)
  results/task4_detrend_summary.json
  figures/task4_detrend_sensitivity.png  rows: hist/ssp245/ssp585,
                                         cols: original, detrended, difference
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from pdo_phase import pattern_corr, HIST_PERIOD
from easm_rain import load_east_china_monthly_mm, jja_anomaly_mmday, composite
from easm_plots import _map_axis, _sym_levels
from task4_future import FUT_PERIOD

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
RES = REPO / "results"
FIG = REPO / "figures"

CASES = [
    # tag, label, prior-results prefix, pr file pattern, period, detrend deg
    ("historical", "historical 1900–2014", "task2",
     "pr_Amon_ACCESS-CM2_historical_{m}i1p1f1_gn_185001-201412.nc",
     HIST_PERIOD, 1),
    ("ssp245", "ssp245 2015–2100", "task4_ssp245",
     "pr_Amon_ACCESS-CM2_ssp245_{m}i1p1f1_gn_201501-210012.nc",
     FUT_PERIOD, 2),
    ("ssp585", "ssp585 2015–2100", "task4_ssp585",
     "pr_Amon_ACCESS-CM2_ssp585_{m}i1p1f1_gn_201501-210012.nc",
     FUT_PERIOD, 2),
]


def phase_years(ds):
    return ([int(y) for y in ds.attrs["pos_years"].split(",")],
            [int(y) for y in ds.attrs["neg_years"].split(",")])


def save_case(tag, members, comps, ens, deg, path):
    ds = xr.Dataset(
        {f"rain_{k}": (("member", "lat", "lon"),
                       np.stack([np.asarray(comps[m][k]) for m in members]))
         for k in ("pos", "neg", "diff", "pval")} |
        {f"ens_{k}": (("lat", "lon"), np.asarray(ens[k]))
         for k in ("pos", "neg", "diff")},
        coords={"member": members,
                "lat": comps[members[0]]["pos"].lat.values,
                "lon": comps[members[0]]["pos"].lon.values},
        attrs={"detrend_deg": deg,
               "description": f"Task 4 detrend-sensitivity composites ({tag}): "
                              "JJA rainfall detrended per gridpoint before "
                              "compositing; PDO phases as in task2/task4."},
    )
    ds.to_netcdf(path)


def sensitivity_figure(rows, path):
    """rows = [(label, orig_diff, dt_diff), ...] -> len(rows) x 3 maps."""
    fig, axes = plt.subplots(len(rows), 3, figsize=(19, 5.2 * len(rows)),
                             constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    lev_ab = _sym_levels(*[np.asarray(f) for r in rows for f in r[1:]])
    lev_d = _sym_levels(*[np.asarray(r[2] - r[1]) for r in rows])
    for i, (label, orig, dt) in enumerate(rows):
        r = pattern_corr(np.asarray(orig), np.asarray(dt))
        for j, (fld, ttl, lev, cmap) in enumerate([
                (orig, "original composite", lev_ab, "BrBG"),
                (dt, f"detrended composite  (r={r:.2f} vs original)",
                 lev_ab, "BrBG"),
                (dt - orig, "detrended $-$ original", lev_d, "RdBu")]):
            ax = axes[i, j]
            cf = ax.contourf(fld.lon, fld.lat, fld, levels=lev, cmap=cmap,
                             transform=ccrs.PlateCarree(), extend="both")
            _map_axis(ax)
            ax.set_title(f"({chr(97+i*3+j)})  {label} — {ttl}", fontsize=11)
            fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.8, pad=0.03)
    fig.suptitle("Task 4 sensitivity — ensemble-mean PDO composite "
                 "(Pos$-$Neg) with vs without rainfall detrending "
                 "(best-3 members)", fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    best3 = json.loads((RES / "task3_summary.json").read_text())["best3"]
    print("members:", best3)

    summary = {"members": best3, "cases": {}}
    rows, ens_dt_diffs = [], {}

    for tag, label, prefix, pr_pat, period, deg in CASES:
        comps, orig_diffs, case_sum = {}, [], {}
        for m in best3:
            prior = xr.open_dataset(RES / f"{prefix}_{m}.nc")
            pos_years, neg_years = phase_years(prior)
            orig_diffs.append(prior["rain_diff"]
                              .rename({"rlat": "lat", "rlon": "lon"}))
            rain_mm = load_east_china_monthly_mm(
                DATA / "ACCESS-CM2_data" / pr_pat.format(m=m), "pr", period)
            comp = composite(jja_anomaly_mmday(rain_mm, detrend_deg=deg),
                             pos_years, neg_years)
            comps[m] = comp
            r = pattern_corr(np.asarray(comp["diff"]),
                             np.asarray(orig_diffs[-1]))
            case_sum[m] = {"diff_corr_orig_vs_detrended": round(float(r), 3)}
            print(f"[{tag} {m}] corr(orig, detrended diff) = {r:.2f}")
        ens = {k: xr.concat([comps[m][k] for m in best3], dim="member")
               .mean("member") for k in ("pos", "neg", "diff")}
        ens_orig = xr.concat(orig_diffs, dim="member").mean("member")
        save_case(tag, best3, comps, ens, deg, RES / f"task4_detrend_{tag}.nc")
        ens_dt_diffs[tag] = ens["diff"]
        r_ens = pattern_corr(np.asarray(ens["diff"]), np.asarray(ens_orig))
        case_sum["ens_diff_corr_orig_vs_detrended"] = round(float(r_ens), 3)
        summary["cases"][tag] = case_sum
        rows.append((label, ens_orig, ens["diff"]))
        print(f"[{tag}] ens corr(orig, detrended) = {r_ens:.2f}")

    # does detrending change the scenario-vs-historical resemblance?
    for tag in ("ssp245", "ssp585"):
        summary["cases"][tag]["ens_diff_corr_vs_hist_detrended"] = round(
            float(pattern_corr(np.asarray(ens_dt_diffs[tag]),
                               np.asarray(ens_dt_diffs["historical"]))), 3)

    sensitivity_figure(rows, FIG / "task4_detrend_sensitivity.png")
    (RES / "task4_detrend_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
