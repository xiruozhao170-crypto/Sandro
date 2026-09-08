"""Task 2 overview: composite rainfall difference (Pos-Neg), obs + r1-r5.

Reads the saved Task 2 results (results/task2_<tag>.nc) and draws the six
composite maps in one 2x3 figure with a single shared colorbar; white dots
mark p<0.1 as in the three-panel figures.

Output: figures/task2_composite_obs_members.png
"""
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from easm_plots import _map_axis, _sym_levels

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"
FIG = REPO / "figures"

TAGS = [("obs", "ERSSTv4 + GPCC (obs)")] + [
    (m, f"ACCESS-CM2 {m}") for m in ["r1", "r2", "r3", "r4", "r5"]
]


def main():
    data = {tag: xr.open_dataset(RES / f"task2_{tag}.nc") for tag, _ in TAGS}
    lev = _sym_levels(*[data[tag]["rain_diff"].values for tag, _ in TAGS])

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for k, (ax, (tag, label)) in enumerate(zip(axes.flat, TAGS)):
        ds = data[tag]
        lo, la = ds["rlon"].values, ds["rlat"].values
        cf = ax.contourf(lo, la, ds["rain_diff"], levels=lev, cmap="BrBG",
                         transform=ccrs.PlateCarree(), extend="both")
        sig = ds["rain_pval"].values < 0.1
        LO, LA = np.meshgrid(lo, la)
        ax.plot(LO[sig], LA[sig], "o", color="white", ms=2.0, mew=0,
                transform=ccrs.PlateCarree())
        _map_axis(ax)
        ax.set_title(f"({chr(97 + k)})  {label}\n"
                     f"pos={ds.attrs['n_pos']}yr  neg={ds.attrs['n_neg']}yr",
                     fontsize=11)
    fig.colorbar(cf, ax=axes.ravel().tolist(), label="mm/day",
                 orientation="horizontal", fraction=0.05, pad=0.03, aspect=45)
    fig.suptitle("Task 2: East China JJA Rainfall Composite (Pos $-$ Neg), "
                 "obs + ACCESS-CM2 r1–r5\n1900–2014  |  9-yr running "
                 "mean phases  |  white dots: p<0.1", fontsize=13)
    fig.savefig(FIG / "task2_composite_obs_members.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("saved figures/task2_composite_obs_members.png")


if __name__ == "__main__":
    main()
