"""Fetch global 850-hPa ACCESS-CM2 winds from NCI's original OPeNDAP arrays.

Raw multi-level files are not downloaded. Each independently restartable request
covers at most 20 years; source versions and subset URLs are recorded.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import time
import xml.etree.ElementTree as ET

import requests
import numpy as np
import xarray as xr

from wpsh_streamfunction import FUTURE_PERIOD, HISTORICAL_PERIOD

ROOT = Path(__file__).resolve().parents[1]
HOST = "https://thredds.nci.org.au/thredds"
NS = {"t": "http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0"}
XLINK = "{http://www.w3.org/1999/xlink}"


def get_xml(url):
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=(30, 90))
            r.raise_for_status()
            return ET.fromstring(r.content)
        except (requests.RequestException, ET.ParseError):
            if attempt == 3:
                raise
            time.sleep(2*(attempt+1))


def catalog(experiment, member, variable):
    activity = "CMIP" if experiment == "historical" else "ScenarioMIP"
    base = (f"{HOST}/catalog/fs38/publications/CMIP6/{activity}/CSIRO-ARCCSS/"
            f"ACCESS-CM2/{experiment}/{member}i1p1f1/Amon/{variable}/gn/")
    tree = get_xml(base + "catalog.xml")
    versions = [x.attrib[XLINK + "title"] for x in tree.findall(".//t:catalogRef", NS)
                if re.fullmatch(r"v\d{8}", x.attrib[XLINK + "title"])]
    if not versions:
        raise ValueError(f"No dated version in {base}")
    version = max(versions)
    tree = get_xml(base + version + "/catalog.xml")
    return [(version, x.attrib["urlPath"]) for x in tree.findall(".//t:dataset", NS)
            if x.attrib.get("urlPath", "").endswith(".nc")]


def download(job):
    path, variable, url_path, first, last, version = job
    url = f"{HOST}/dodsC/{url_path}"
    # Do not use this server's NCSS output: NetCDF-Java rewrites pressure as
    # the mean of its bounds (85000 Pa becomes 83125 Pa). OPeNDAP preserves
    # the original pressure and monthly timestamps.
    if not path.exists():
        part = path.with_suffix(".part.nc")
        for attempt in range(6):
            try:
                print(f"Downloading {path.name}", flush=True)
                with xr.open_dataset(url, engine="pydap") as remote:
                    iz = np.flatnonzero(np.isclose(remote.plev, 85000, rtol=0, atol=1))
                    if len(iz) != 1 or remote.plev.attrs.get("units") != "Pa":
                        raise ValueError("Original file does not have an exact 85000-Pa level")
                    da = remote[variable].isel(plev=[int(iz[0])]).sel(time=slice(str(first), str(last))).load()
                    # PyDAP may leave CMOR's 1e20 sentinel undecoded in memory.
                    da = da.where(abs(da) < 1e10)
                    subset = da.to_dataset()
                    subset.attrs.update(remote.attrs)
                    subset.attrs["source_url"] = url
                    subset.attrs["below_ground_winds"] = "Original missing values preserved; no imputation"
                    subset.to_netcdf(part, encoding={variable: dict(zlib=True, complevel=4)})
                with xr.open_dataset(part) as ds:
                    if ds.sizes["time"] != 12 * (last-first+1):
                        raise ValueError(f"Incomplete monthly subset {path.name}")
                    if abs(float(ds.plev.item()) - 85000) > 1:
                        raise ValueError("Subset selected wrong pressure")
                    if ds.sizes["lon"] != 192 or ds.sizes["lat"] != 145:
                        raise ValueError("Incomplete ACCESS-CM2 global grid")
                part.replace(path)
                break
            # PyDAP raises IndexError/ValueError on truncated binary responses.
            # Retry the whole atomic subset; never accept an incomplete array.
            except (requests.RequestException, OSError, IndexError, ValueError) as exc:
                if attempt == 5:
                    raise
                print(f"Retry {path.name}: {type(exc).__name__}: {exc}", flush=True)
                time.sleep(2*(attempt+1))
    with path.open("rb") as f:
        digest = hashlib.file_digest(f, "sha256").hexdigest()
    print(f"Ready {path.name}", flush=True)
    return dict(file=path.name, version=version, source_file=url_path,
                opendap_url=url, first_year=first, last_year=last,
                sha256=digest, pressure_Pa=85000,
                missing_values="Preserved; explicitly select --wind-fill for a documented extension")


def subset_year_ranges(first, last):
    """Keep 2100 separate so existing subsets ending in 2099 remain reusable."""
    for a, b in ((first, min(last, 2099)), (max(first, 2100), last)):
        for start in range(a, b+1, 20):
            yield start, min(b, start+19)


def execute_downloads(jobs, out, workers, existing=()):
    # NetCDF/HDF5 library access is not thread safe on many Windows builds.
    from concurrent.futures import ProcessPoolExecutor, as_completed
    records = {r["file"]: r for r in existing}
    failures = []
    print(f"Prepared {len(jobs)} subsets", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download, job): job for job in jobs}
        for future in as_completed(futures):
            try:
                record = future.result()
                records[record["file"]] = record
                manifest_part = out / "download_manifest.part.json"
                manifest_part.write_text(json.dumps(sorted(records.values(), key=lambda r: r["file"]), indent=2)+"\n")
                manifest_part.replace(out / "download_manifest.json")
            except Exception as exc:
                failures.append((str(futures[future][0]), repr(exc)))
                print(f"FAILED {failures[-1]}", flush=True)
    if failures:
        raise RuntimeError(f"Subsets failed after retries; rerun to resume: {failures}")


def extend_to_2100(out, workers):
    """Append 2100 from the exact same source versions as the existing winds."""
    records = json.loads((out / "download_manifest.json").read_text())
    tails = [r for r in records if r["last_year"] == 2099]
    expected = {(v, s, m) for v in ("ua", "va") for s in ("ssp245", "ssp585") for m in ("r1", "r2", "r5")}
    jobs, found = [], set()
    for record in tails:
        stem = Path(record["file"]).stem.split("_")
        variable, scenario, member = stem[0][:2], stem[1], stem[2]
        found.add((variable, scenario, member))
        match = re.search(r"_(\d{4})\d{2}-(\d{4})\d{2}\.nc$", record["source_file"])
        if not match or not int(match[1]) <= 2100 <= int(match[2]):
            raise ValueError(f"Pinned source does not cover 2100: {record['source_file']}")
        path = out / f"{variable}850_{scenario}_{member}_2100_2100.nc"
        jobs.append((path, variable, record["source_file"], 2100, 2100, record["version"]))
    if found != expected or len(jobs) != len(expected):
        raise ValueError("Need all twelve existing scenario u/v subsets ending in 2099")
    execute_downloads(jobs, out, workers, existing=records)


def main(out, workers):
    out.mkdir(parents=True, exist_ok=True)
    jobs = []
    # Catalog discovery is pure HTTP/XML: threads avoid serial round trips.
    from concurrent.futures import ThreadPoolExecutor
    keys = [(e, m, v) for e in ("historical", "ssp245", "ssp585")
            for m in ("r1", "r2", "r5") for v in ("ua", "va")]
    with ThreadPoolExecutor(max_workers=6) as pool:
        catalogs = dict(zip(keys, pool.map(lambda key: catalog(*key), keys)))
    for experiment in ("historical", "ssp245", "ssp585"):
        first, last = HISTORICAL_PERIOD if experiment == "historical" else FUTURE_PERIOD
        for member in ("r1", "r2", "r5"):
            for variable in ("ua", "va"):
                for version, source in catalogs[(experiment, member, variable)]:
                    match = re.search(r"_(\d{4})\d{2}-(\d{4})\d{2}\.nc$", source)
                    a, b = max(first, int(match[1])), min(last, int(match[2]))
                    for start, end in subset_year_ranges(a, b):
                        path = out / f"{variable}850_{experiment}_{member}_{start}_{end}.nc"
                        jobs.append((path, variable, source, start, end, version))
    execute_downloads(jobs, out, workers)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=ROOT / "data/wpsh_winds/access")
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--extend-to-2100", action="store_true",
                   help="Append 2100 using source versions pinned in the existing download manifest")
    args = p.parse_args()
    if args.extend_to_2100:
        extend_to_2100(args.out, args.workers)
    else:
        main(args.out, args.workers)
