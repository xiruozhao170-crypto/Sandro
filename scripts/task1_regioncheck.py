"""Draw six Task 1 regional checks: ERSSTv4 observations and ACCESS R1-R5.

Run without arguments to redraw from the small committed NetCDF cache.
Use --data-root to rebuild ACCESS SST and --ersst-file for ERSSTv4 SST.
"""
import argparse
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

ROOT = Path(__file__).resolve().parents[1]
MEMBERS = [f"r{i}" for i in range(1, 6)]
CACHE = ROOT / "results/task1_regioncheck_access.nc"
OBS_CACHE = ROOT / "results/task1_regioncheck_ersst.nc"
OUTPUT = ROOT / "figures/task1_regioncheck_obs_members.png"


def save_fields(fields):
    """Save only the month used by the original region checks, in member order."""
    sst = xr.concat([fields[m] for m in MEMBERS if m in fields],
                    dim=xr.IndexVariable("member", [m for m in MEMBERS if m in fields]))
    sst.name = "sst"
    sst.attrs.update(units="degC", description="August 1900 SST; month #8 of 1900-2014")
    ds = sst.to_dataset()
    ds.attrs.update(description="Task 1 regional checks, ACCESS-CM2 historical r1-r5",
                    region="20-60N, 110-250E", month="1900-08",
                    regrid="IDW k=4, unit-sphere chord distances; regrid_tos.py")
    CACHE.parent.mkdir(exist_ok=True)
    ds.to_netcdf(CACHE)


def update_member(da, member):
    """Called by task1_pdo.py; update the member cache and combined figure."""
    fields = {}
    if CACHE.exists():
        with xr.open_dataset(CACHE) as ds:
            fields = {str(m): ds.sst.sel(member=m, drop=True).load()
                      for m in ds.member.values}
    fields[member] = da.isel(time=7, drop=True).load()
    save_fields(fields)
    plot_if_ready()


def save_observation(field):
    """Cache the same August 1900 observational field as the original check."""
    field = field.sortby("lat").sel(lat=slice(20, 60), lon=slice(110, 250))
    ds = field.rename("sst").to_dataset()
    ds.attrs.update(description="Task 1 regional check, ERSSTv4 observations",
                    region="20-60N, 110-250E", month="1900-08")
    OBS_CACHE.parent.mkdir(exist_ok=True)
    ds.to_netcdf(OBS_CACHE)


def plot_if_ready():
    """Allow Task 1 to run one dataset at a time before all six are cached."""
    if CACHE.exists() and OBS_CACHE.exists():
        with xr.open_dataset(CACHE) as ds:
            ready = all(m in ds.member.values for m in MEMBERS)
        if ready:
            plot_checks()


def rebuild_cache(data_root):
    """Use the existing regridding recipe, reading just the required month."""
    from regrid_tos import build_weights, TGT_LAT, TGT_LON

    fields = {}
    weights = None
    for member in MEMBERS:
        processed = data_root / f"processed/tos_2deg_historical_{member}i1p1f1.nc"
        if processed.exists():
            with xr.open_dataset(processed) as ds:
                field = ds.tos.sel(time="1900-08").isel(time=0, drop=True)
                fields[member] = field.sel(lat=slice(20, 60), lon=slice(110, 250)).load()
        else:
            raw = data_root / "ACCESS-CM2_data" / (
                f"tos_Omon_ACCESS-CM2_historical_{member}i1p1f1_gn_185001-201412.nc")
            if weights is None:
                weights = build_weights(raw)
            mask, idx, w = weights
            with xr.open_dataset(raw) as ds:
                month = ds.tos.sel(time="1900-08").isel(time=0).values
            values = np.einsum("pk,pk->p", month[mask][idx], w).astype(np.float32)
            field = xr.DataArray(values.reshape(len(TGT_LAT), len(TGT_LON)),
                                 dims=("lat", "lon"),
                                 coords={"lat": TGT_LAT, "lon": TGT_LON})
            fields[member] = field.sel(lat=slice(20, 60), lon=slice(110, 250))
        print(f"Loaded {member}: August 1900", flush=True)
    save_fields(fields)


def plot_checks():
    with xr.open_dataset(CACHE) as source:
        ds = source.sel(member=MEMBERS).load()
    with xr.open_dataset(OBS_CACHE) as source:
        obs = source.sst.load()
    panels = [("ERSSTv4 (obs)", obs)] + [
        (f"ACCESS-CM2 {m.upper()}", ds.sst.sel(member=m)) for m in MEMBERS]
    lo = np.floor(min(float(field.min()) for _, field in panels) / 5) * 5
    hi = np.ceil(max(float(field.max()) for _, field in panels) / 5) * 5
    fig = plt.figure(figsize=(14, 8.8))
    grid = fig.add_gridspec(3, 4, left=.055, right=.96, bottom=.16,
                           top=.88, hspace=.38, wspace=.38)
    slots = [grid[row, col:col+2] for row in range(3) for col in (0, 2)]
    for i, ((label, field), slot) in enumerate(zip(panels, slots)):
        ax = fig.add_subplot(slot, projection=ccrs.PlateCarree(central_longitude=180))
        cf = ax.pcolormesh(field.lon, field.lat, field, cmap="RdYlBu_r", vmin=lo, vmax=hi,
                          transform=ccrs.PlateCarree(), shading="auto")
        ax.coastlines(linewidth=.7)
        ax.add_feature(cfeature.LAND, facecolor="0.85")
        ax.set_extent([109, 251, 19, 61], crs=ccrs.PlateCarree())
        ax.set_xticks([120, 160, 200, 240], crs=ccrs.PlateCarree())
        ax.set_yticks([20, 40, 60], crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.tick_params(labelsize=9)
        ax.set_title(f"({chr(97+i)})  {label}", loc="left", fontsize=12)
    cax = fig.add_axes([.27, .075, .46, .025])
    fig.colorbar(cf, cax=cax, orientation="horizontal", label="SST (°C)",
                 ticks=np.arange(lo, hi + 1, 5))
    fig.suptitle("Task 1: North Pacific regional checks | ERSSTv4 + ACCESS-CM2 R1–R5\n"
                 "August 1900 (month #8) | 20–60°N, 110–250°E", fontsize=15, y=.97)
    OUTPUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUTPUT, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, help="Rebuild ACCESS SST from this data directory")
    parser.add_argument("--ersst-file", type=Path, help="Rebuild observation cache from ERSSTv4 monthly SST")
    args = parser.parse_args()
    if args.data_root:
        rebuild_cache(args.data_root)
    if args.ersst_file:
        with xr.open_dataset(args.ersst_file) as ds:
            save_observation(ds.sst.sel(time="1900-08").isel(time=0, drop=True).load())
    plot_checks()
