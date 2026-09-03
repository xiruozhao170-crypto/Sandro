"""Task 2: composite East-China JJA rainfall for PDO+/- phases.

Datasets: observations (ERSSTv4 PDO + GPCC v2018 rainfall) and the five
ACCESS-CM2 historical members (tos 2-deg regrid + pr interpolated onto the
GPCC 0.5-degree grid). 1900-2014. The PDO is always EOF1 (Task 1.docx);
phases come from the JJA mean of the 9-yr smoothed PC1 (Task 2.docx).

Outputs per tag in {obs, r1..r5}:
  results/task2_<tag>.nc                 composites + PDO index/phases
  figures/task2_<tag>_three_panels.png
Plus figures/task2_pdo_timeseries.png (obs + all members).
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr

from pdo_phase import get_pdo, HIST_PERIOD
from easm_rain import (gpcc_grid, load_east_china_monthly_mm,
                       jja_anomaly_mmday, composite)
from easm_plots import rainfall_three_panels, pdo_timeseries_panels

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
RES = REPO / "results"
FIG = REPO / "figures"

MEMBERS = ["r1", "r2", "r3", "r4", "r5"]


def save_dataset(tag, pdo, comp, path):
    ds = xr.Dataset(
        {
            "pdo_eof": (("lat", "lon"), pdo["eof"]),
            "pdo_pc": (("time",), pdo["pc"]),
            "pdo_pc_smooth": (("time",), pdo["pc_smooth"].values),
            "pdo_jja_smooth": (("year_cls",), pdo["yearly"].values),
            "rain_pos": (("rlat", "rlon"), comp["pos"].values),
            "rain_neg": (("rlat", "rlon"), comp["neg"].values),
            "rain_diff": (("rlat", "rlon"), comp["diff"].values),
            "rain_pval": (("rlat", "rlon"), comp["pval"].values),
        },
        coords={
            "lat": pdo["lat"], "lon": pdo["lon"], "time": pdo["time"],
            "year_cls": pdo["yearly"].index.values,
            "rlat": comp["pos"].lat.values, "rlon": comp["pos"].lon.values,
        },
        attrs={
            "pdo_eof1_varfrac": float(pdo["varfrac"][0]),
            "pdo_corr_with_obs": float(pdo["corr_ref"])
            if np.isfinite(pdo["corr_ref"]) else -999.0,
            "pos_years": ",".join(map(str, pdo["pos_years"])),
            "neg_years": ",".join(map(str, pdo["neg_years"])),
            "n_pos": comp["n_pos"], "n_neg": comp["n_neg"],
            "description": f"Task 2 composites ({tag}), PDO = EOF1, "
                           "9-yr running-mean JJA phases, rainfall in mm/day "
                           "on the GPCC 0.5-deg grid",
        },
    )
    ds.to_netcdf(path)


def run(tag, tos_path, tos_var, rain_path, rain_var, label,
        ref_eof1=None, target_grid=None):
    pdo = get_pdo(tos_path, tos_var, period=HIST_PERIOD, ref_eof1=ref_eof1)
    rain_mm = load_east_china_monthly_mm(rain_path, rain_var, HIST_PERIOD,
                                         target_grid=target_grid)
    anom = jja_anomaly_mmday(rain_mm)
    comp = composite(anom, pdo["pos_years"], pdo["neg_years"])
    save_dataset(tag, pdo, comp, RES / f"task2_{tag}.nc")
    rainfall_three_panels(
        comp,
        f"East China JJA Rainfall — PDO Phase Analysis  |  {label}\n"
        "1900–2014  |  9-yr running mean",
        FIG / f"task2_{tag}_three_panels.png",
    )
    print(f"[{tag}] EOF1 var={pdo['varfrac'][0]*100:.1f}% "
          f"corr_obs={pdo['corr_ref']:.2f} "
          f"pos={comp['n_pos']}yr neg={comp['n_neg']}yr")
    return pdo, comp


def main():
    entries, summary = [], {}
    grid = gpcc_grid()

    obs_pdo, obs_comp = run(
        "obs", DATA / "ersst.v4.sst.mnmean.nc", "sst",
        DATA / "precip.mon.total.v2018.nc", "precip",
        "ERSSTv4 + GPCC v2018 (obs)")
    entries.append(("(a)  ERSSTv4 PDO PC1", obs_pdo))
    summary["obs"] = {"n_pos": obs_comp["n_pos"], "n_neg": obs_comp["n_neg"],
                      "eof1_varfrac": round(float(obs_pdo["varfrac"][0]), 4)}

    for i, m in enumerate(MEMBERS):
        pdo, comp = run(
            m,
            DATA / f"processed/tos_2deg_historical_{m}i1p1f1.nc", "tos",
            DATA / f"ACCESS-CM2_data/pr_Amon_ACCESS-CM2_historical_{m}i1p1f1_gn_185001-201412.nc",
            "pr",
            f"ACCESS-CM2 historical {m}",
            ref_eof1=obs_pdo["eof"], target_grid=grid)
        entries.append((f"({chr(98+i)})  ACCESS-CM2 {m} PDO PC1", pdo))
        summary[m] = {
            "eof1_varfrac": round(float(pdo["varfrac"][0]), 4),
            "eof1_corr_with_obs": round(float(pdo["corr_ref"]), 3),
            "n_pos": comp["n_pos"], "n_neg": comp["n_neg"],
        }

    pdo_timeseries_panels(entries, FIG / "task2_pdo_timeseries.png")
    (RES / "task2_summary.json").write_text(json.dumps(summary, indent=2))
    print("saved figures/task2_pdo_timeseries.png, results/task2_summary.json")


if __name__ == "__main__":
    main()
