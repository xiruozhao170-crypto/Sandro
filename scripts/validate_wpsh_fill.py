"""Artificially mask ERA5 with ACCESS topographic masks and measure WPSH errors.

ERA5 1979-2020 JJA is regridded to ACCESS's native global grid. Successive
ACCESS historical r1 JJA masks from 1900-1941 are applied month by month.
All three branches use identical global inversion, interpolation and EOF steps.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
import xarray as xr

from task7_wpsh_eof import open_winds
from wpsh_streamfunction import canonical_wind, eof_analysis, sha256, streamfunction
from wpsh_wind_fill import FILL_VERSION, WindFiller

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/task7_wpsh_eof"


def main():
    raw = ROOT / "data/wpsh_winds/era5_uv850_1979_2020_0.25deg.nc"
    files = sorted((ROOT / "data/wpsh_winds/access").glob("ua850_historical_r1_*.nc"))
    files = [p for p in files if ".part." not in p.name and int(p.stem.split("_")[-2]) <= 1941]
    masks, opened = open_winds(files, "ua", 1900, 1941)
    lat, lon = masks.lat.values, masks.lon.values
    filler = WindFiller(lat, lon)
    methods = ["complete", "harmonic", "nearest"]
    fields = {m: [] for m in methods}
    mask_fractions = []
    iy = np.flatnonzero((lat >= 0) & (lat <= 50))
    ix = np.flatnonzero((lon >= 90) & (lon <= lon[lon >= 200].min()))
    with xr.open_dataset(raw) as ds:
        u, v = canonical_wind(ds, "u"), canonical_wind(ds, "v")
        for year in range(1979, 2021):
            take = (u.time.dt.year == year) & u.time.dt.month.isin([6, 7, 8])
            uw = u.where(take, drop=True).interp(lat=lat, lon=lon).values.astype(float)
            vw = v.where(take, drop=True).interp(lat=lat, lon=lon).values.astype(float)
            mask_year = year-1979+1900
            missing = np.isnan(masks.where((masks.time.dt.year == mask_year) & masks.time.dt.month.isin([6, 7, 8]), drop=True).values)
            if missing.shape != uw.shape:
                raise ValueError("Incomplete benchmark masks")
            mask_fractions.extend(missing.mean(axis=(1, 2)).tolist())
            batch_u, batch_v = [uw], [vw]
            for method in methods[1:]:
                filled = [filler.fill(np.where(missing[k], np.nan, uw[k]),
                                      np.where(missing[k], np.nan, vw[k]), method) for k in range(3)]
                batch_u.append(np.stack([x[0] for x in filled]))
                batch_v.append(np.stack([x[1] for x in filled]))
            sf = streamfunction(np.concatenate(batch_u), np.concatenate(batch_v), lat, lon)
            for k, method in enumerate(methods):
                fields[method].append(sf[k*3:k*3+3].mean(axis=0)[iy][:, ix])
            print(f"ERA5 masking benchmark: {year}", flush=True)
    for ds in opened:
        ds.close()
    years = np.arange(1979, 2021)
    design = np.polynomial.polynomial.polyvander(np.linspace(-1, 1, len(years)), 1)
    products, eofs = {}, {}
    for method in methods:
        values = np.stack(fields[method]).astype(float)
        flat = values.reshape(len(years), -1)
        values = (flat-design @ np.linalg.lstsq(design, flat, rcond=None)[0]).reshape(values.shape)
        da = xr.DataArray(values, dims=("year", "lat", "lon"), coords=dict(year=years, lat=lat[iy], lon=lon[ix]))
        products[method] = da.sortby("lat").interp(lat=np.arange(0, 51), lon=np.arange(90, 201))
        eofs[method] = eof_analysis(products[method], [1979, 2020])
    complete = eofs["complete"]
    report = dict(description=__doc__, fill_version=FILL_VERSION,
                  missing_grid_fraction_mean=float(np.mean(mask_fractions)),
                  inputs=[dict(file=p.name, sha256=sha256(p)) for p in [raw]+files], results={})
    table = dict(year=years)
    original_pc = complete.pc.sel(mode=1).values
    original_years = set(years[original_pc > 1.5].tolist())
    for method in methods:
        eof = eofs[method]
        pc = eof.pc.sel(mode=1).values
        table[f"pc1_{method}"] = pc
        field = products[method].sel(lat=slice(10, 40), lon=slice(110, 180))
        truth = products["complete"].sel(lat=field.lat, lon=field.lon)
        error = np.sqrt(float(((field-truth)**2).mean()))
        report["results"][method] = dict(eof1_variance_percent=float(eof.variance_fraction.sel(mode=1)*100),
            pc1_correlation_with_complete=float(np.corrcoef(pc, original_pc)[0, 1]),
            pc1_rmse_sd=float(np.sqrt(np.mean((pc-original_pc)**2))),
            sf_anomaly_rmse_m2_s=error, sf_anomaly_relative_rmse=error/np.sqrt(float((truth**2).mean())),
            strong_event_years=years[pc > 1.5].tolist(),
            changed_event_years=sorted(set(years[pc > 1.5].tolist()) ^ original_years))
    (OUT / "wind_fill_benchmark.json").write_text(json.dumps(report, indent=2)+"\n")
    pd.DataFrame(table).to_csv(OUT / "wind_fill_benchmark_index.csv", index=False)
    print(json.dumps(report["results"], indent=2), flush=True)


if __name__ == "__main__":
    main()
