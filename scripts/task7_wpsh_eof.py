"""WPSH streamfunction/EOF calculation following Yang et al. (2022), PNAS.

Run --source era5, ncep, or access after the corresponding download script.
Use --redraw to regenerate figures using only the committed result NetCDFs.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import xarray as xr

from wpsh_streamfunction import (FUTURE_PERIOD, HISTORICAL_PERIOD, canonical_wind,
                                compare_periods, eof_analysis, regional_products, sha256)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/task7_wpsh_eof"
FIGURES = ROOT / "figures/task7_wpsh_eof"
REFERENCE = "https://doi.org/10.1073/pnas.2120335119"


def open_winds(paths, variable, first, last):
    arrays, opened = [], []
    for path in paths:
        ds = xr.open_dataset(path)
        opened.append(ds)
        arrays.append(canonical_wind(ds, variable))
    if not arrays:
        raise FileNotFoundError(f"No input files for {variable}")
    wind = xr.concat(arrays, dim="time", join="exact", coords="minimal", compat="equals").sortby("time")
    wind = wind.sel(time=slice(str(first), str(last)))
    return wind, opened


def run(tag, label, ufiles, vfiles, uvar, vvar, first, last, baseline,
        degree, detrend_order="seasonal", regrid=False, wind_fill="reject"):
    RESULTS.mkdir(parents=True, exist_ok=True)
    cache = ROOT / f"data/wpsh_winds/cache/{tag}.nc"
    cache.parent.mkdir(parents=True, exist_ok=True)
    inputs = [dict(file=str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else p.name,
                   sha256=sha256(p), bytes=p.stat().st_size) for p in sorted(set(ufiles+vfiles))]
    if tag.startswith("access_"):
        for item, path in zip(inputs, sorted(set(ufiles+vfiles))):
            with xr.open_dataset(path) as source:
                item["source_url"] = source.attrs["source_url"]
                item["tracking_id"] = source.attrs.get("tracking_id", "")
                item["pressure_Pa"] = float(source.plev.item())
    config = dict(first=first, last=last, baseline=baseline, degree=degree, inputs=inputs)
    if wind_fill != "reject":
        from wpsh_wind_fill import FILL_VERSION
        config.update(wind_fill=wind_fill, fill_version=FILL_VERSION)
    fingerprint = json.dumps(config, sort_keys=True)
    if cache.exists():
        with xr.open_dataset(cache) as cached:
            valid = cached.attrs.get("input_configuration") == fingerprint
            products = cached.load() if valid else None
    else:
        products = None
    if products is None:
        u, uds = open_winds(ufiles, uvar, first, last)
        v, vds = open_winds(vfiles, vvar, first, last)
        try:
            products = regional_products(u, v, first, last, baseline, degree, wind_fill=wind_fill)
        finally:
            for ds in uds+vds:
                ds.close()
        products.attrs["input_configuration"] = fingerprint
        products.to_netcdf(cache, encoding={x: dict(zlib=True, complevel=4) for x in products.data_vars})
    if regrid:
        # Paper's common 1-degree grid for models, after the global inversion.
        products = products.interp(lat=np.arange(0, 51), lon=np.arange(90, 201))
        products.attrs["output_regrid"] = "Bilinear regional streamfunction interpolation to 1 degree"
    field = products[f"sf_detrended_{detrend_order}"]
    eof = eof_analysis(field, baseline)
    result = xr.merge([products, eof], combine_attrs="drop_conflicts")
    result.attrs.update(products.attrs)
    result.attrs.update(eof.attrs)
    result.attrs.update(label=label, tag=tag, reference=REFERENCE, detrend_order=detrend_order,
                        main_index="PC1; positive is stronger anticyclonic WPSH",
                        note="First EOF is retained without mode searching; weighting/order sensitivities are reported separately")
    # The paper does not specify EOF area weights or every processing detail.
    # Preserve these sensitivity calculations rather than tuning to 77%.
    other_order = "monthly" if detrend_order == "seasonal" else "seasonal"
    alt = eof_analysis(products[f"sf_detrended_{other_order}"], baseline)
    unweighted = eof_analysis(field, baseline, weighting="none")
    pc = result.pc.sel(mode=1).values
    summary = dict(tag=tag, label=label, first=first, last=last, baseline=baseline,
                   wind_fill=wind_fill,
                   detrend_degree=degree, detrend_order=detrend_order,
                   eof1_variance_percent=float(result.variance_fraction.sel(mode=1)*100),
                   eof2_variance_percent=float(result.variance_fraction.sel(mode=2)*100),
                   pc1_regional_sf_correlation=float(np.corrcoef(pc, result.regional_sf)[0, 1]),
                   strong_years_full_sd=result.year.values[pc > 1.5].astype(int).tolist(),
                   strong_years_reference_sd=result.year.values[result.pc_reference.sel(mode=1).values > 1.5].astype(int).tolist(),
                   sensitivity=dict(other_detrend_order=other_order,
                       other_order_eof1_variance_percent=float(alt.variance_fraction.sel(mode=1)*100),
                       other_order_pc1_correlation=float(np.corrcoef(pc, alt.pc.sel(mode=1))[0, 1]),
                       unweighted_eof1_variance_percent=float(unweighted.variance_fraction.sel(mode=1)*100),
                       unweighted_pc1_correlation=float(np.corrcoef(pc, unweighted.pc.sel(mode=1))[0, 1])))
    if first == 1900 and last == 2099:
        summary["period_comparison"] = compare_periods(result.pc_reference.sel(mode=1))
        summary["period_comparison"]["threshold"] = "PC1 > 1.5 historical-period sample SD, after historical-period centering"
        summary["period_comparison_full_period_threshold"] = compare_periods(result.pc.sel(mode=1))
        summary["period_comparison_full_period_threshold"]["threshold"] = "PC1 > 1.5 full-1900-2099-period sample SD"
    if wind_fill != "reject":
        summary["missing_winds"] = {k: products.attrs[k] for k in
            ("missing_grid_fraction_mean", "missing_grid_fraction_min", "missing_grid_fraction_max")}
        summary["missing_winds"]["caveat"] = "Below-ground winds are a numerical extension, not observed or recovered model winds"
    result["pc1_other_detrend_order"] = alt.pc.sel(mode=1, drop=True)
    result["pc1_unweighted"] = unweighted.pc.sel(mode=1, drop=True)
    result["strong_event"] = (result.pc_reference.sel(mode=1, drop=True) > 1.5).astype("int8")
    result.strong_event.attrs["definition"] = "PC1 standardized against the reference period > 1.5"
    encoding = {x: dict(zlib=True, complevel=4) for x in result.data_vars}
    for name in ("sf_jja", "sf_climatology", "sf_anomaly", "sf_detrended_monthly", "sf_detrended_seasonal", "regression"):
        encoding[name]["dtype"] = "float32"
    result.attrs["storage_precision"] = "Spatial fields float32; PC indices and statistics float64"
    result.to_netcdf(RESULTS / f"{tag}.nc", encoding=encoding)
    pd.DataFrame(dict(year=result.year.values, pc1=pc, pc2=result.pc.sel(mode=2).values,
                      pc1_reference=result.pc_reference.sel(mode=1).values,
                      regional_sf_m2_s=result.regional_sf.values,
                      regional_sf_reference=result.regional_sf_reference.values,
                      strong_event=result.strong_event.values)).to_csv(RESULTS / f"{tag}_index.csv", index=False)
    (RESULTS / f"{tag}_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    versions = {}
    for name in ("numpy", "xarray", "scipy", "netCDF4", "skyborn", "windspharm", "matplotlib"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    manifest = dict(configuration=config, source=label, reference=REFERENCE,
                    python=platform.python_version(), packages=versions)
    if tag.startswith("era5"):
        request_path = ufiles[0].parent / "era5_download_request.json"
        manifest["source_catalogue"] = "https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels-monthly-means"
    elif tag.startswith("ncep"):
        request_path = ufiles[0].parent / "ncep_download_manifest.json"
        manifest["source_catalogue"] = "https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis.derived/pressure/catalog.html"
    else:
        request_path = ufiles[0].parent / "download_manifest.json"
    if request_path.exists():
        manifest["download_request"] = json.loads(request_path.read_text())
        if tag.startswith("access_"):
            used_names = {p.name for p in ufiles+vfiles}
            manifest["download_request"] = [r for r in manifest["download_request"] if r["file"] in used_names]
    (RESULTS / f"{tag}_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(summary, indent=2), flush=True)
    if wind_fill != "nearest":
        plot_result(result)


def map_panel(fig, spec, data, title, limit):
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
    ax = fig.add_subplot(spec, projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_extent([90, 200, 0, 50], ccrs.PlateCarree())
    mesh = ax.contourf(data.lon, data.lat, data/1e6, levels=np.linspace(-limit, limit, 17),
                       cmap="RdBu_r", extend="both", transform=ccrs.PlateCarree())
    ax.coastlines(resolution="110m", linewidth=.6, color="#444444")
    ax.add_patch(Rectangle((110, 10), 70, 30, fill=False, linewidth=1.2,
                          transform=ccrs.PlateCarree(), edgecolor="#202020"))
    ax.set_xticks([100, 120, 140, 160, 180, 200], crs=ccrs.PlateCarree())
    ax.set_yticks([0, 10, 20, 30, 40, 50], crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(labelsize=8)
    ax.set_title(title, fontsize=11, loc="left", pad=9)
    return ax, mesh


def plot_result(ds):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.8, 8.2), layout="constrained")
    grid = fig.add_gridspec(3, 2, height_ratios=[1, .10, .86])
    limit = float(np.ceil(abs(ds.regression.values/1e6).max()*2)/2)
    for i in range(2):
        _, mesh = map_panel(fig, grid[0, i], ds.regression.isel(mode=i),
            f"({'ab'[i]}) EOF{i+1}: {float(ds.variance_fraction.isel(mode=i))*100:.1f}% variance", limit)
    cax = fig.add_subplot(grid[1, :])
    fig.colorbar(mesh, cax=cax, orientation="horizontal", label=r"Streamfunction regression per 1 SD PC ($10^6$ m$^2$ s$^{-1}$)")
    ax = fig.add_subplot(grid[2, :])
    pc = ds.pc_reference.sel(mode=1)
    if ds.attrs["tag"].startswith("access_"):
        ax.axvspan(HISTORICAL_PERIOD[0], HISTORICAL_PERIOD[1]+.5, color="#eeeeee", zorder=0)
        ax.axvline(HISTORICAL_PERIOD[1]+.5, color="#777777", lw=.7)
    ax.plot(ds.year, pc, color="#202020", lw=1.5, label="WPSH PC1")
    ax.plot(ds.year, ds.regional_sf_reference, color="#b44831", lw=1.05, alpha=.9,
            label="Regional streamfunction (15-25 N, 120-150 E)")
    ax.axhline(0, color="#999999", lw=.6)
    ax.axhline(1.5, color="#6b5598", ls="--", lw=1.1, label="Strong event: PC1 > 1.5 reference SD")
    strong = pc > 1.5
    ax.scatter(ds.year.where(strong), pc.where(strong), s=17, color="#6b5598", zorder=4)
    if len(ds.year) < 60:
        for y, value in zip(ds.year.values[strong], pc.values[strong]):
            ax.annotate(str(int(y)), (y, value), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=8)
    ax.set_ylim(min(-3, float(pc.min())-.3, float(ds.regional_sf_reference.min())-.3),
                max(3.5, float(pc.max())+.6, float(ds.regional_sf_reference.max())+.3))
    ax.set_xlim(int(ds.year.min())-1, int(ds.year.max())+1)
    ax.set_ylabel("Reference-period SD")
    ax.set_title("(c) JJA WPSH index", loc="left", fontsize=11)
    ax.legend(loc="lower left", fontsize=8, ncol=3)
    ax.grid(axis="y", alpha=.15)
    degree = "Linear" if ds.attrs["detrend_degree"] == 1 else "Quadratic"
    fig.suptitle(f"{ds.attrs['label']} | 850-hPa JJA WPSH\n"
                 f"{int(ds.year.min())}-{int(ds.year.max())}; {degree.lower()} {ds.attrs['detrend_order']} detrending; "
                 f"reference {ds.attrs['baseline_start']}-{ds.attrs['baseline_end']}", fontsize=14)
    for ext in ("png", "svg"):
        target = FIGURES / f"{ds.attrs['tag']}.{ext}"
        fig.savefig(target, dpi=200, bbox_inches="tight")
        if ext == "svg":
            target.write_text("\n".join(line.rstrip() for line in target.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", choices=["era5", "ncep", "access"], default="era5")
    p.add_argument("--data", type=Path, default=ROOT / "data/wpsh_winds")
    p.add_argument("--members", nargs="+", default=["r1", "r2", "r5"], choices=["r1", "r2", "r5"])
    p.add_argument("--scenarios", nargs="+", default=["ssp245", "ssp585"], choices=["ssp245", "ssp585"])
    p.add_argument("--redraw", action="store_true")
    p.add_argument("--wind-fill", choices=["reject", "harmonic", "nearest"], default="reject",
                   help="Explicit below-ground extension, only used for ACCESS")
    args = p.parse_args()
    if args.redraw:
        for file in sorted(RESULTS.glob("*.nc")):
            with xr.open_dataset(file) as ds:
                plot_result(ds.load())
        return
    if args.source == "access":
        data = args.data / "access"
        for scenario in args.scenarios:
            for member in args.members:
                ufiles = sorted(data.glob(f"ua850_historical_{member}_*.nc")) + sorted(data.glob(f"ua850_{scenario}_{member}_*.nc"))
                vfiles = sorted(data.glob(f"va850_historical_{member}_*.nc")) + sorted(data.glob(f"va850_{scenario}_{member}_*.nc"))
                ufiles = [p for p in ufiles if ".part." not in p.name]
                vfiles = [p for p in vfiles if ".part." not in p.name]
                suffix = "" if args.wind_fill == "harmonic" else f"_{args.wind_fill}"
                run(f"access_{scenario}_{member}{suffix}",
                    f"ACCESS-CM2 {member} | historical + {scenario.upper()} | {args.wind_fill} wind extension",
                    ufiles, vfiles, "ua", "va", HISTORICAL_PERIOD[0], FUTURE_PERIOD[1],
                    list(HISTORICAL_PERIOD), 2, "monthly", regrid=True,
                    wind_fill=args.wind_fill)
    elif args.source == "ncep":
        run("ncep_1979_2020", "NCEP/NCAR R1 (2.5-degree grid)",
            [args.data / "ncep_uwnd_850_1979_2020.nc"], [args.data / "ncep_vwnd_850_1979_2020.nc"],
            "uwnd", "vwnd", 1979, 2020, [1979, 2020], 1)
    else:
        files = [args.data / "era5_uv850_1979_2020_0.25deg.nc"]
        run("era5_1979_2020", "ERA5 (0.25-degree grid)", files, files, "u", "v", 1979, 2020, [1979, 2020], 1)


if __name__ == "__main__":
    main()
