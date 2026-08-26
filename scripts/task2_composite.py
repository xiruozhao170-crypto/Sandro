"""Task 2: composite East-China JJA rainfall for PDO+/- phases.

Datasets: observations (ERSSTv4 PDO + GPCC v2018 rainfall) and the five
ACCESS-CM2 historical members (tos 2-deg regrid + native pr). 1900-2014.

For the model members the PDO mode is selected as the EOF (of the first 3)
best matching the observed PDO pattern (handles the r4 EOF1/EOF2 swap).

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
from easm_rain import load_east_china_monthly_mm, jja_anomaly_mmday, composite
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
            "pdo_mode": int(pdo["mode"]) + 1,
            "pdo_corr_with_obs": float(pdo["corr_ref"])
            if np.isfinite(pdo["corr_ref"]) else -999.0,
            "pos_years": ",".join(map(str, pdo["pos_years"])),
            "neg_years": ",".join(map(str, pdo["neg_years"])),
            "n_pos": comp["n_pos"], "n_neg": comp["n_neg"],
            "description": f"Task 2 composites ({tag}), 1900-2014, "
                           "9-yr running-mean PDO phases, rainfall in mm/day",
        },
    )
    ds.to_netcdf(path)


def run(tag, tos_path, tos_var, rain_path, rain_var, ref_eof1, label):
    pdo = get_pdo(tos_path, tos_var, ref_eof1=ref_eof1, period=HIST_PERIOD)
    rain_mm = load_east_china_monthly_mm(rain_path, rain_var, HIST_PERIOD)
    anom = jja_anomaly_mmday(rain_mm)
    comp = composite(anom, pdo["pos_years"], pdo["neg_years"])
    save_dataset(tag, pdo, comp, RES / f"task2_{tag}.nc")
    rainfall_three_panels(
        comp,
        f"East China JJA Rainfall — PDO Phase Analysis  |  {label}\n"
        "1900–2014  |  9-yr running mean",
        FIG / f"task2_{tag}_three_panels.png",
    )
    print(f"[{tag}] PDO mode=EOF{pdo['mode']+1} corr_obs={pdo['corr_ref']:.2f} "
          f"pos={comp['n_pos']}yr neg={comp['n_neg']}yr")
    return pdo, comp


def main():
    entries, summary = [], {}

    obs_pdo, obs_comp = run(
        "obs", DATA / "ersst.v4.sst.mnmean.nc", "sst",
        DATA / "precip.mon.total.v2018.nc", "precip",
        None, "ERSSTv4 + GPCC v2018 (obs)")
    entries.append(("(a)  ERSSTv4 PDO", obs_pdo))
    summary["obs"] = {"n_pos": obs_comp["n_pos"], "n_neg": obs_comp["n_neg"]}

    for i, m in enumerate(MEMBERS):
        pdo, comp = run(
            m,
            DATA / f"processed/tos_2deg_historical_{m}i1p1f1.nc", "tos",
            DATA / f"ACCESS-CM2_data/pr_Amon_ACCESS-CM2_historical_{m}i1p1f1_gn_185001-201412.nc",
            "pr",
            obs_pdo["eof"], f"ACCESS-CM2 historical {m}")
        entries.append((f"({chr(98+i)})  ACCESS-CM2 {m} PDO (EOF{pdo['mode']+1})", pdo))
        summary[m] = {
            "pdo_mode": int(pdo["mode"]) + 1,
            "pdo_corr_with_obs": round(float(pdo["corr_ref"]), 3),
            "n_pos": comp["n_pos"], "n_neg": comp["n_neg"],
        }

    pdo_timeseries_panels(entries, FIG / "task2_pdo_timeseries.png")
    (RES / "task2_summary.json").write_text(json.dumps(summary, indent=2))
    print("saved figures/task2_pdo_timeseries.png, results/task2_summary.json")


if __name__ == "__main__":
    main()
