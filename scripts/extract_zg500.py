"""Extract the 500-hPa level from a raw CMIP6 zg file (data prep for Task 6).

Usage: python extract_zg500.py <raw_zg.nc>

Writes data/processed/zg500_<experiment>_<member>_<t0>-<t1>.nc with only
plev=50000 Pa, lon 60-240E, lat 15S-70N (float32, compressed). The output
name is taken from the file's own global attributes, not the input filename.
"""
import sys
from pathlib import Path

import xarray as xr

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "processed"

LON = slice(60, 240)
LAT = slice(-15, 70)


def main(path):
    ds = xr.open_dataset(path)
    exp = ds.attrs["experiment_id"]
    member = ds.attrs["variant_label"]
    da = ds["zg"].sel(plev=50000.0, method="nearest").sel(lat=LAT, lon=LON)
    assert abs(float(da.plev) - 50000.0) < 1.0, float(da.plev)
    da = da.drop_vars("plev")
    t = ds["time"].dt.strftime("%Y%m").values
    out = OUT / f"zg500_{exp}_{member}_{t[0]}-{t[-1]}.nc"
    da.load().astype("float32").to_dataset(name="zg500").to_netcdf(
        out, encoding={"zg500": {"zlib": True, "complevel": 4}})
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1])
