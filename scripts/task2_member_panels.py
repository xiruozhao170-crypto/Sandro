"""Redraw Task 2 R1-R5 as a 5-by-3 figure from the saved member results."""
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from easm_plots import _map_axis, _sym_levels

ROOT = Path(__file__).resolve().parents[1]
MEMBERS = [f"r{i}" for i in range(1, 6)]
OUTPUT = ROOT / "figures/task1-4/task2_members_three_panels.png"


def main():
    data = {}
    for member in MEMBERS:
        with xr.open_dataset(ROOT / f"results/task2_{member}.nc") as ds:
            data[member] = ds.load()
    variables = ["rain_pos", "rain_neg", "rain_diff"]
    levels = _sym_levels(*[data[m][v].values for m in MEMBERS for v in variables])
    fig, axes = plt.subplots(5, 3, figsize=(12.5, 20),
                             subplot_kw={"projection": ccrs.PlateCarree()})
    fig.subplots_adjust(left=.085, right=.97, bottom=.085, top=.925,
                        hspace=.24, wspace=.18)
    titles = ["PDO positive phase", "PDO negative phase", "Composite (Pos − Neg)"]
    for row, member in enumerate(MEMBERS):
        ds = data[member]
        lon, lat = ds.rlon.values, ds.rlat.values
        for col, variable in enumerate(variables):
            ax = axes[row, col]
            cf = ax.contourf(lon, lat, ds[variable], levels=levels, cmap="BrBG",
                             transform=ccrs.PlateCarree(), extend="both")
            if col == 2:
                sig = ds.rain_pval.values < .1
                xx, yy = np.meshgrid(lon, lat)
                ax.plot(xx[sig], yy[sig], "o", color="white", ms=1.7, mew=0,
                        transform=ccrs.PlateCarree())
            _map_axis(ax)
            ax.tick_params(labelsize=8, labelleft=col == 0)
            ax.set_title(f"({chr(97 + row*3 + col)})  {member.upper()}",
                         loc="left", fontsize=10, pad=5)
            if row == 0:
                ax.text(.5, 1.15, titles[col], transform=ax.transAxes,
                        ha="center", fontsize=12, fontweight="bold")
        axes[row, 0].text(-.23, .5,
                          f"{member.upper()}\npos={ds.attrs['n_pos']} yr\nneg={ds.attrs['n_neg']} yr",
                          transform=axes[row, 0].transAxes, ha="center", va="center",
                          fontsize=10, rotation=90)
    cax = fig.add_axes([.24, .048, .56, .012])
    fig.colorbar(cf, cax=cax, orientation="horizontal",
                 label="JJA rainfall anomaly / composite difference (mm/day)")
    fig.suptitle("Task 2: East China JJA rainfall | ACCESS-CM2 historical R1–R5\n"
                 "1900–2014 | 9-yr running-mean PDO phases | white dots: p < 0.1",
                 fontsize=14, y=.985)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    main()
