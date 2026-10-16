"""Download global monthly 850-hPa winds; never subset to the EOF box first."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import requests
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]


def fetch_ncep(out: Path, first=1979, last=2020):
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for var in ("uwnd", "vwnd"):
        url = ("https://psl.noaa.gov/thredds/ncss/grid/Datasets/"
               f"ncep.reanalysis.derived/pressure/{var}.mon.mean.nc")
        params = dict(var=var, north=90, south=-90, west=0, east=357.5,
                      horizStride=1, vertCoord=850, timeStride=1,
                      time_start=f"{first}-01-01T00:00:00Z",
                      time_end=f"{last}-12-01T00:00:00Z", accept="netcdf4")
        path = out / f"ncep_{var}_850_{first}_{last}.nc"
        request_url = requests.Request("GET", url, params=params).prepare().url
        if not path.exists():
            print(f"Downloading {path.name}", flush=True)
            part = path.with_suffix(".part.nc")
            for attempt in range(4):
                try:
                    with requests.get(url, params=params, stream=True, timeout=(30, 300)) as r:
                        r.raise_for_status()
                        with part.open("wb") as f:
                            for chunk in r.iter_content(1024 * 1024):
                                f.write(chunk)
                    break
                except requests.RequestException:
                    if attempt == 3:
                        raise
                    time.sleep(2)
            with xr.open_dataset(part) as ds:
                if ds.sizes["time"] != (last - first + 1) * 12:
                    raise ValueError("Incomplete NOAA monthly download")
                if float(ds.level.item()) != 850 or ds.sizes["lon"] != 144:
                    raise ValueError("Wrong pressure or incomplete longitude coverage")
            part.replace(path)
        with path.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256").hexdigest()
        records.append(dict(file=path.name, url=request_url, sha256=digest,
                            source="NCEP/NCAR Reanalysis 1", pressure_hPa=850))
        print(f"Ready: {path.name}", flush=True)
    (out / "ncep_download_manifest.json").write_text(json.dumps(records, indent=2)+"\n")


def fetch_era5(out: Path, first=1979, last=2020, grid=0.25):
    import cdsapi
    out.mkdir(parents=True, exist_ok=True)
    request = dict(product_type=["monthly_averaged_reanalysis"],
                   variable=["u_component_of_wind", "v_component_of_wind"],
                   pressure_level=["850"], year=[str(y) for y in range(first, last+1)],
                   month=[f"{m:02d}" for m in range(1, 13)], time=["00:00"],
                   data_format="netcdf", download_format="unarchived", grid=[grid, grid])
    path = out / f"era5_uv850_{first}_{last}_{grid:g}deg.nc"
    if not path.exists():
        part = path.with_suffix(".part.nc")
        cdsapi.Client().retrieve("reanalysis-era5-pressure-levels-monthly-means",
                                 request, str(part))
        with xr.open_dataset(part) as ds:
            if not {"u", "v"}.issubset(ds):
                raise ValueError("ERA5 download does not contain both u and v")
        part.replace(path)
    (out / "era5_download_request.json").write_text(json.dumps(request, indent=2)+"\n")
    print(f"Ready: {path.name}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", choices=["ncep", "era5"])
    parser.add_argument("--out", type=Path, default=ROOT / "data/wpsh_winds")
    parser.add_argument("--first", type=int, default=1979)
    parser.add_argument("--last", type=int, default=2020)
    parser.add_argument("--grid", type=float, default=0.25, help="ERA5 output grid in degrees")
    args = parser.parse_args()
    if args.source == "ncep":
        fetch_ncep(args.out, args.first, args.last)
    else:
        fetch_era5(args.out, args.first, args.last, args.grid)
