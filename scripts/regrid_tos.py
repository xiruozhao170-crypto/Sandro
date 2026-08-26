"""Regrid ACCESS-CM2 tos from the native tripolar ocean grid (j,i + 2-D
lat/lon) to the regular 2x2 deg ERSST-style grid, so the Task-1 recipe
(lat/lon box selection, cos(lat) weights) applies directly.

Method: inverse-distance weighting of the 4 nearest ocean grid points,
with distances measured as chords on the unit sphere (no longitude-wrap
issues). Weights are built once and reused for all 15 files. Target
points farther than ~2.5 deg from any ocean point are left as NaN (land).

Output: data/processed/tos_2deg_<experiment>_<member>.nc
"""
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "data/ACCESS-CM2_data"
OUT = REPO / "data/processed"
OUT.mkdir(exist_ok=True)

# ERSST-style global 2-deg grid
TGT_LAT = np.arange(-88, 89, 2.0)
TGT_LON = np.arange(0, 359, 2.0)
K = 4
MAX_DIST_DEG = 2.5  # reject targets with no ocean point within this radius


def to_xyz(lat, lon):
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    return np.column_stack(
        [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
    )


def build_weights(sample_file):
    ds = xr.open_dataset(sample_file)
    lat2d = ds["latitude"].values
    lon2d = ds["longitude"].values
    mask = np.isfinite(ds["tos"].isel(time=0).values)  # ocean mask (constant)
    src_xyz = to_xyz(lat2d[mask], lon2d[mask])
    tree = cKDTree(src_xyz)

    glon, glat = np.meshgrid(TGT_LON, TGT_LAT)
    tgt_xyz = to_xyz(glat.ravel(), glon.ravel())
    dist, idx = tree.query(tgt_xyz, k=K)

    w = 1.0 / np.maximum(dist, 1e-10)
    w /= w.sum(axis=1, keepdims=True)
    max_chord = 2 * np.sin(np.deg2rad(MAX_DIST_DEG) / 2)
    land = dist[:, 0] > max_chord
    w[land] = np.nan
    print(f"weights: {mask.sum()} ocean src pts -> {land.size} tgt pts "
          f"({land.sum()} land/NaN)")
    return mask, idx, w


def regrid_file(path, mask, idx, w):
    ds = xr.open_dataset(path)
    nt = ds.sizes["time"]
    out = np.empty((nt, TGT_LAT.size * TGT_LON.size), dtype=np.float32)
    for t0 in range(0, nt, 240):  # chunk over time to bound memory
        t1 = min(t0 + 240, nt)
        block = ds["tos"].isel(time=slice(t0, t1)).values  # (t, j, i)
        flat = block[:, mask]                              # (t, npts)
        out[t0:t1] = np.einsum("tpk,pk->tp", flat[:, idx], w)
    da = xr.DataArray(
        out.reshape(nt, TGT_LAT.size, TGT_LON.size),
        dims=("time", "lat", "lon"),
        coords={"time": ds.time, "lat": TGT_LAT, "lon": TGT_LON},
        name="tos",
        attrs={"units": ds["tos"].attrs.get("units", "degC"),
               "regrid": "IDW k=4 from ACCESS-CM2 native ocean grid"},
    )
    return da


def main():
    files = sorted(SRC.glob("tos_Omon_*.nc"))
    mask, idx, w = build_weights(files[0])
    for f in files:
        parts = f.name.split("_")          # tos Omon ACCESS-CM2 exp member gn dates
        exp, member = parts[3], parts[4]
        outfile = OUT / f"tos_2deg_{exp}_{member}.nc"
        if outfile.exists():
            print("skip (exists):", outfile.name)
            continue
        da = regrid_file(f, mask, idx, w)
        da.to_netcdf(outfile, encoding={"tos": {"zlib": True, "complevel": 4}})
        print("wrote", outfile.name, dict(da.sizes))


if __name__ == "__main__":
    main()
