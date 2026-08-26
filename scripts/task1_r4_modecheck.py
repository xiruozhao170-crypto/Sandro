"""Task 1 addendum: member r4's leading EOF is not the canonical PDO.

For access_r4 the first three EOFs are computed and each is pattern-correlated
with the observed (ERSSTv4) PDO. Result (2026-08-25 run): EOF1 |r|=0.21,
EOF2 |r|=0.85, EOF3 |r|=0.33 -> in r4 the PDO appears as mode 2 (mode swap
with a Kuroshio-Oyashio-Extension-dominated mode; EOF1 17.3% vs EOF2 14.0%
are close, i.e. near-degenerate by North's rule).
Appends the finding to results/task1_summary.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from eofs.standard import Eof

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from task1_pdo import DATASETS, load_north_pacific  # noqa: E402
from task1_compare import pattern_corr  # noqa: E402


def main():
    cfg = DATASETS["access_r4"]
    da = load_north_pacific(cfg["path"], cfg["var"]).load()
    fit = da.polyfit(dim="time", deg=1)
    da_dt = da - xr.polyval(da.time, fit.polyfit_coefficients)
    clim = da_dt.groupby("time.month").mean("time")
    anom = da_dt.groupby("time.month") - clim
    w = np.sqrt(np.cos(np.deg2rad(anom.lat.values)))[:, None]
    solver = Eof(anom.values.astype("float64"), weights=w)
    eofs3 = solver.eofsAsCovariance(neofs=3, pcscaling=1)
    vf = solver.varianceFraction(neigs=3)

    obs = xr.open_dataset(REPO / "results/task1_pdo_ersst.nc")
    entry = {}
    for k in range(3):
        pat = xr.DataArray(eofs3[k], coords={"lat": anom.lat, "lon": anom.lon},
                           dims=("lat", "lon"))
        r = pattern_corr(pat, obs.eof1)
        entry[f"EOF{k+1}"] = {"variance_fraction_pct": round(float(vf[k]) * 100, 1),
                              "abs_pattern_corr_vs_obs": round(abs(r), 3)}
        print(f"r4 EOF{k+1}: {vf[k]*100:.1f}%  |r|={abs(r):.3f}")

    p = REPO / "results/task1_summary.json"
    summary = json.loads(p.read_text())
    summary["note_r4_mode_swap"] = {
        "finding": "In member r4 the canonical PDO appears as EOF2, not EOF1; "
                   "EOF1 is a KOE-dominated mode (EOF1/EOF2 near-degenerate).",
        "modes": entry,
    }
    p.write_text(json.dumps(summary, indent=2))
    print("updated", p.name)


if __name__ == "__main__":
    main()
