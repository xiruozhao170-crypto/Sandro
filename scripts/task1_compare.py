"""Task 1 wrap-up: compare the ACCESS-CM2 PDO against ERSSTv4 observations.

Reads results/task1_pdo_<tag>.nc produced by task1_pdo.py and reports, per
ensemble member: EOF1 variance fraction and the centered pattern correlation
with the observed PDO pattern (both are on the common 2-deg grid).
Writes results/task1_summary.json and a members overview figure.
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"
FIG = REPO / "figures"

MEMBERS = ["r1", "r2", "r3", "r4", "r5"]


def pattern_corr(a, b):
    """Centered pattern correlation with cos(lat) area weights."""
    a, b = xr.align(a.sortby("lat"), b.sortby("lat"))
    w = np.cos(np.deg2rad(a.lat)) * xr.ones_like(a.lon, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    aw, bw, ww = a.where(m), b.where(m), w.where(m)
    am = (aw * ww).sum() / ww.sum()
    bm = (bw * ww).sum() / ww.sum()
    ap, bp = aw - am, bw - bm
    num = (ww * ap * bp).sum()
    den = np.sqrt((ww * ap**2).sum() * (ww * bp**2).sum())
    return float(num / den)


def main():
    obs = xr.open_dataset(RES / "task1_pdo_ersst.nc")
    summary = {
        "ersst": {"variance_fraction_pct": round(float(obs.variance_fraction[0]) * 100, 1)}
    }

    fig = plt.figure(figsize=(16, 6.5))
    panels = [("ERSSTv4 (obs)", obs)] + [
        (f"ACCESS-CM2 {m}", xr.open_dataset(RES / f"task1_pdo_access_{m}.nc"))
        for m in MEMBERS
    ]
    for k, (label, ds) in enumerate(panels):
        vf = float(ds.variance_fraction[0]) * 100
        r = pattern_corr(ds.eof1, obs.eof1) if label != "ERSSTv4 (obs)" else 1.0
        if label != "ERSSTv4 (obs)":
            summary[label] = {
                "variance_fraction_pct": round(vf, 1),
                "pattern_corr_vs_obs": round(r, 3),
            }
        ax = fig.add_subplot(2, 3, k + 1,
                             projection=ccrs.PlateCarree(central_longitude=180))
        lim = float(np.nanmax(np.abs(ds.eof1.values)))
        cf = ax.contourf(ds.lon, ds.lat, ds.eof1,
                         levels=np.linspace(-lim, lim, 21), cmap="RdBu_r",
                         transform=ccrs.PlateCarree(), extend="both")
        ax.coastlines()
        ax.add_feature(cfeature.LAND, facecolor="0.85")
        title = f"{label}\nEOF1 {vf:.1f}%"
        if label != "ERSSTv4 (obs)":
            title += f", r(obs)={r:.2f}"
        ax.set_title(title, fontsize=10)
        fig.colorbar(cf, ax=ax, shrink=0.8)
    fig.suptitle("Task 1: PDO pattern, observations vs ACCESS-CM2 historical members "
                 "(1900-2014)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "task1_pdo_members.png", dpi=200, bbox_inches="tight")

    (RES / "task1_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
