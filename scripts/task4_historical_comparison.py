"""Combine Task 4-7/4-8 from saved member results, retaining both original PNGs.

Rows: SSP245/SSP585; columns: R1-R5. Each field is the future PDO-phase
rainfall composite (positive minus negative) minus the same member's
historical composite. All ten maps use a shared symmetric color scale.
"""
from pathlib import Path

import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from easm_plots import _map_axis, _sym_levels

ROOT = Path(__file__).resolve().parents[1]
MEMBERS = [f"r{i}" for i in range(1, 6)]
SCENARIOS = ["ssp245", "ssp585"]
OUTPUT = ROOT / "figures/task1-4/Task4_T4-7-8_Composite_Difference_SSP_minus_Historical.png"


def load_differences():
    fields = {}
    for member in MEMBERS:
        with xr.open_dataset(ROOT / f"results/task2_{member}.nc") as ds:
            historical = ds.rain_diff.load()
        for scenario in SCENARIOS:
            with xr.open_dataset(ROOT / f"results/task4_{scenario}_{member}.nc") as ds:
                future, baseline = xr.align(ds.rain_diff, historical, join="exact")
                fields[scenario, member] = (future - baseline).load()
    return fields


def main():
    fields = load_differences()
    levels = _sym_levels(*[field.values for field in fields.values()])
    fig, axes = plt.subplots(2, 5, figsize=(22, 10.5),
                             subplot_kw={"projection": ccrs.PlateCarree()})
    fig.subplots_adjust(left=.035, right=.99, bottom=.14, top=.87,
                        wspace=.20, hspace=.28)
    for row, scenario in enumerate(SCENARIOS):
        for col, member in enumerate(MEMBERS):
            ax = axes[row, col]
            field = fields[scenario, member]
            cf = ax.contourf(field.rlon, field.rlat, field, levels=levels,
                             cmap="BrBG", transform=ccrs.PlateCarree(), extend="both")
            _map_axis(ax)
            ax.set_title(f"({chr(97 + row*5 + col)})  {member.upper()}  |  "
                         f"{scenario.upper()} − historical", fontsize=11)
    cax = fig.add_axes([.24, .065, .52, .025])
    fig.colorbar(cf, cax=cax, label="Change in composite difference (mm/day)",
                 orientation="horizontal")
    fig.suptitle("Task 4-7 / 4-8: Change in JJA rainfall composite (Pos − Neg) | ACCESS-CM2\n"
                 "Future (2015–2100) − same-member historical (1900–2014)", fontsize=15, y=.97)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    main()
