"""JJA 500-hPa / 5870-gpm WPSH contours, following Dong (2016), Fig. 4(a).

PDO phases are reused and independently checked against pdo_phase.py.
All ENSO states are retained, as explicitly requested for this adaptation.
Raw CMIP6 files stay outside git. Run with --data-root pointing to the 14 files.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import netCDF4
import numpy as np
import pandas as pd
import requests
import xarray as xr

from pdo_phase import pdo_index_and_phases

REPO = Path(__file__).resolve().parents[1]
MEMBERS = ("r1", "r2", "r5")
EXPERIMENTS = ("historical", "ssp245", "ssp585")
LEVEL = 5870.0
EXTENT = (90, 200, 0, 50)
NCSS = ("https://psl.noaa.gov/thredds/ncss/grid/Datasets/"
        "ncep.reanalysis.derived/pressure/hgt.mon.mean.nc")
POS, NEG = "#b52230", "#2459a6"


def sha256(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def download_observations(path):
    """Public NOAA NCEP/NCAR Reanalysis 1, subset without interpolation."""
    params = dict(var="hgt", north=50, west=90, east=200, south=0,
                  horizStride=1, time_start="1948-01-01T00:00:00Z",
                  time_end="2014-12-01T00:00:00Z", timeStride=1,
                  vertCoord=500, accept="netcdf4")
    url = requests.Request("GET", NCSS, params=params).prepare().url
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        response = requests.get(url, timeout=180)
        response.raise_for_status()
        temporary = path.with_suffix(".part.nc")
        temporary.write_bytes(response.content)
        with xr.open_dataset(temporary) as ds:
            assert "hgt" in ds and float(ds.level.item()) == 500
            assert ds.sizes["time"] == 804
        temporary.replace(path)
    return dict(source="NOAA PSL NCEP/NCAR Reanalysis 1", url=url,
                file=path.name, sha256=sha256(path))


def inventory(data_root):
    records = []
    paths = sorted(data_root.glob("zg_Amon_ACCESS-CM2_*.nc"))
    if len(paths) != 14:
        raise ValueError(f"Expected the 14 specified zg files; found {len(paths)}")
    for path in paths:
        with netCDF4.Dataset(path) as ds:
            var, t = ds["zg"], ds["time"]
            dates = netCDF4.num2date(t[:], t.units, t.calendar)
            ym = np.array([d.year * 12 + d.month for d in dates])
            assert np.all(np.diff(ym) == 1), path.name
            pressure = ds["plev"][:]
            assert ds["plev"].units == "Pa"
            assert np.count_nonzero(np.isclose(pressure, 50000, atol=.001, rtol=0)) == 1
            assert var.units == "m" and var.standard_name == "geopotential_height"
            assert ds.source_id == "ACCESS-CM2"
            assert ds.experiment_id in EXPERIMENTS
            assert ds.variant_label in [m + "i1p1f1" for m in MEMBERS]
            first, last = dates[0].year, dates[-1].year
            used = first <= 2100
            records.append(dict(file=path.name, bytes=path.stat().st_size,
                                experiment=ds.experiment_id, member=ds.variant_label,
                                first_month=f"{first:04d}-{dates[0].month:02d}",
                                last_month=f"{last:04d}-{dates[-1].month:02d}",
                                n_months=len(dates), units=var.units,
                                selected_pressure_Pa=float(pressure[np.argmin(abs(pressure-50000))]),
                                tracking_id=getattr(ds, "tracking_id", ""),
                                used_for_composite=used,
                                note="" if used else "Beyond 2100; no matching repository PDO phases or other-member scenario coverage."))
    return records


def phase_series(path):
    """Reuse repository PC1, asserting exact agreement with its phase function."""
    with xr.open_dataset(path) as ds:
        _, smooth, yearly, pos, neg = pdo_index_and_phases(ds.pdo_pc.values, ds.time.values)
        np.testing.assert_allclose(smooth.values, ds.pdo_pc_smooth.values, equal_nan=True)
        np.testing.assert_array_equal(yearly.index.values, ds.year_cls.values)
        np.testing.assert_allclose(yearly.values, ds.pdo_jja_smooth.values)
        assert list(pos) == list(map(int, ds.attrs["pos_years"].split(",")))
        assert list(neg) == list(map(int, ds.attrs["neg_years"].split(",")))
    return yearly, sha256(path)


def model_jja(data_root, experiment, member, cache):
    """Read only JJA at 500 hPa, using each monthly chunk just once."""
    first, last = (1900, 2014) if experiment == "historical" else (2015, 2100)
    paths = sorted(data_root.glob(f"zg_Amon_ACCESS-CM2_{experiment}_{member}i1p1f1_*.nc"))
    used_paths = [p for p in paths if int(p.stem.split("_")[-1][:4]) <= last]
    fingerprint = json.dumps([(p.name, p.stat().st_size, p.stat().st_mtime_ns)
                              for p in used_paths])
    fingerprint = hashlib.sha256((fingerprint + "jja500-v1").encode()).hexdigest()
    if cache.exists():
        with xr.open_dataset(cache) as ds:
            if ds.attrs.get("source_fingerprint") == fingerprint:
                print(f"Cache: {experiment} {member}", flush=True)
                return ds.z500_jja.load()
    annual, years = [], []
    for path in used_paths:
        print(f"Read JJA / 500 hPa: {path.name}", flush=True)
        with netCDF4.Dataset(path) as ds:
            t = ds["time"]
            dates = netCDF4.num2date(t[:], t.units, t.calendar)
            lat, lon = ds["lat"][:], ds["lon"][:]
            ilat = np.flatnonzero((lat >= -2.5) & (lat <= 52.5))
            ilon = np.flatnonzero((lon >= 88) & (lon <= 202))
            iy = slice(ilat[0], ilat[-1] + 1)
            ix = slice(ilon[0], ilon[-1] + 1)
            iz = int(np.argmin(abs(ds["plev"][:] - 50000)))
            for year in sorted({d.year for d in dates if first <= d.year <= last}):
                indices = [i for i, d in enumerate(dates) if d.year == year and d.month in (6, 7, 8)]
                assert [dates[i].month for i in indices] == [6, 7, 8]
                months = np.stack([np.ma.filled(ds["zg"][i, iz, iy, ix], np.nan)
                                   for i in indices]).astype("float64")
                assert np.isfinite(months).all(), (path.name, year)
                assert 4000 < months.min() < months.max() < 7000
                annual.append(months.mean(axis=0))
                years.append(year)
    np.testing.assert_array_equal(years, np.arange(first, last + 1))
    result = xr.DataArray(np.stack(annual), dims=("year", "lat", "lon"),
                          coords=dict(year=years, lat=lat[iy], lon=lon[ix]), name="z500_jja",
                          attrs=dict(units="m", standard_name="geopotential_height",
                                     averaging="Arithmetic mean of June, July, August monthly means"))
    ds = result.to_dataset()
    ds.attrs["source_fingerprint"] = fingerprint
    cache.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(cache, encoding={"z500_jja": dict(zlib=True, complevel=4)})
    return result


def observation_jja(path):
    with xr.open_dataset(path) as ds:
        assert ds.hgt.attrs["units"] == "m" and float(ds.level.item()) == 500
        z = ds.hgt.sel(level=500).sortby("lat")
        z = z.where(z.time.dt.month.isin([6, 7, 8]), drop=True).load().astype("float64")
        np.testing.assert_array_equal(z.groupby("time.year").count("time"),
                                      np.full((67, z.sizes["lat"], z.sizes["lon"]), 3))
        out = z.groupby("time.year").mean("time")
    assert np.isfinite(out).all()
    return out.rename("z500_jja")


def composite(z, phases, tag, phase_path, phase_hash):
    aligned = phases.reindex(z.year.values)
    pos = aligned.index[aligned > 0].to_numpy(dtype=int)
    neg = aligned.index[aligned < 0].to_numpy(dtype=int)
    assert len(pos) and len(neg)
    ds = xr.Dataset(dict(z500_jja=z, z500_pos=z.sel(year=pos).mean("year"),
                         z500_neg=z.sel(year=neg).mean("year"),
                         pdo_jja_smooth=("year", aligned.values),
                         phase=("year", np.where(aligned > 0, 1, np.where(aligned < 0, -1, 0)))))
    ds["z500_diff"] = ds.z500_pos - ds.z500_neg
    for name in ("z500_jja", "z500_pos", "z500_neg", "z500_diff"):
        ds[name].attrs.update(units="m", long_name=name.replace("_", " "))
    ds.phase.attrs["definition"] = "+1=PDO positive, -1=PDO negative, 0=unclassified (smoothing edge or exact zero)"
    ds.pdo_jja_smooth.attrs["units"] = "1"
    ds.attrs.update(dataset=tag, pressure_hPa=500., wpsh_contour_gpm=LEVEL,
                    season="JJA", n_pos=len(pos), n_neg=len(neg),
                    pos_years=",".join(map(str, pos)), neg_years=",".join(map(str, neg)),
                    phase_file=str(phase_path.relative_to(REPO)), phase_file_sha256=phase_hash,
                    phase_method="Repository EOF1 PC; 108-month centered mean; sign of JJA mean",
                    enso_filter="None (explicit user instruction)",
                    height_processing="Absolute geopotential height; no detrending, anomalies or bias correction",
                    temporal_averaging="Equal monthly weights within JJA; equal year weights within each phase",
                    requested_period=f"{int(z.year.min())}-{int(z.year.max())}",
                    used_period=f"{min(pos.min(),neg.min())}-{max(pos.max(),neg.max())}")
    assert np.isfinite(ds.z500_pos).all() and np.isfinite(ds.z500_neg).all()
    return ds


def ensemble(datasets, tag):
    """Equal mean of the three member-specific phase composites, not pooled years."""
    ds = xr.Dataset({v: xr.concat([d[v] for d in datasets], dim="member").mean("member")
                     for v in ("z500_pos", "z500_neg", "z500_diff")})
    ds.attrs.update(dataset=tag, pressure_hPa=500., wpsh_contour_gpm=LEVEL,
                    members=",".join(MEMBERS), averaging="Equal member weights AFTER each member's PDO phase compositing",
                    requested_period=datasets[0].attrs["requested_period"],
                    used_period=datasets[0].attrs["used_period"],
                    enso_filter="None (explicit user instruction)")
    return ds


def map_axis(ax):
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="#f0f0ec", zorder=0)
    ax.coastlines(resolution="110m", linewidth=.65, color="#555555")
    ax.set_xticks([100, 120, 140, 160, 180, 200], crs=ccrs.PlateCarree())
    ax.set_yticks([0, 10, 20, 30, 40, 50], crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(labelsize=9, length=3)
    ax.gridlines(xlocs=[100, 120, 140, 160, 180, 200], ylocs=[0, 10, 20, 30, 40, 50],
                 linewidth=.35, color="gray", alpha=.4, linestyle=":")


def draw_panel(ax, ds, title, phases=("pos", "neg")):
    map_axis(ax)
    unavailable = []
    for phase in phases:
        name = f"z500_{phase}"
        color = POS if phase == "pos" else NEG
        # Separate-phase panels use continuous lines; overlay retains the paper's styles.
        style = (0, (6, 4)) if phase == "neg" and len(phases) == 2 else "solid"
        field = ds[name]
        if not float(field.min()) < LEVEL < float(field.max()):
            label = "PDO+" if name.endswith("pos") else "PDO-"
            unavailable.append((f"{label}: no 5870-gpm contour\n"
                                f"regional maximum = {float(field.max()):.1f} gpm", color))
            continue
        ax.contour(field.lon, field.lat, field, levels=[LEVEL], colors=[color],
                   linestyles=[style], linewidths=1.35 if phase == "neg" else 1.65,
                   transform=ccrs.PlateCarree(), zorder=3)
    for i, (message, color) in enumerate(unavailable):
        ax.text(.97, .10 + i * .17, message, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=9, color=color, zorder=5,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.88, pad=3))
    period = ds.attrs["requested_period"]
    counts = (f"PDO+ n={ds.attrs['n_pos']}  |  PDO- n={ds.attrs['n_neg']}"
              if "n_pos" in ds.attrs else "r1/r2/r5 mean")
    if len(phases) == 1 and "n_pos" in ds.attrs:
        counts = f"n={ds.attrs['n_' + phases[0]]} years"
    ax.set_title(f"{title}\n{period}  |  {counts}", loc="left", fontsize=10.5, pad=9)


def save_figure(fig, path, title, subtitle, separated=False):
    fig.suptitle(title, fontsize=17, fontweight="semibold", y=.98)
    fig.text(.5, .935, subtitle, ha="center", fontsize=10, color="#444444")
    handles = [Line2D([], [], color=POS, lw=1.65, label="PDO positive: 5870 gpm", linestyle="-"),
               Line2D([], [], color=NEG, lw=1.35, label="PDO negative: 5870 gpm",
                      linestyle="-" if separated else (0, (6, 4)))]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .035),
               ncol=2, frameon=False, fontsize=11)
    fig.text(.5, .016, "5870-gpm height contours | All ENSO states | Lines meeting the map frame continue outside the view",
             ha="center", fontsize=9, color="#555555")
    fig.subplots_adjust(left=.055, right=.975, bottom=.11, top=.87, wspace=.23, hspace=.49)
    path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".svg"):
        destination = path.with_suffix(suffix)
        metadata = {"Date": None} if suffix == ".svg" else None
        fig.savefig(destination, dpi=300, facecolor="white", metadata=metadata)
        if suffix == ".svg":
            destination.write_text("\n".join(line.rstrip() for line in destination.read_text(encoding="utf-8").splitlines()) + "\n",
                                   encoding="utf-8")
    plt.close(fig)


def contour_display_audit(collections):
    """Audit physical absence, missing values and clipping without changing fields."""
    import contourpy
    from shapely.geometry import LineString, box

    rows = []
    boxes = {"old": box(100, 0, 180, 40), "new": box(90, 0, 200, 50)}
    for key, ds in collections.items():
        for phase in ("pos", "neg"):
            field = ds[f"z500_{phase}"]
            assert np.isfinite(field).all(), (key, phase)
            segments = contourpy.contour_generator(x=field.lon.values, y=field.lat.values,
                                                   z=field.values).lines(LEVEL)
            lines = [LineString(segment) for segment in segments if len(segment) > 1]
            bounds = [field.lon.min().item(), field.lon.max().item(),
                      field.lat.min().item(), field.lat.max().item()]
            status = ("below_threshold" if float(field.max()) < LEVEL else
                      "above_threshold" if float(field.min()) > LEVEL else
                      "contour_present" if lines else "no_resolved_contour")
            row = dict(dataset=key, phase=phase, min_gpm=float(field.min()), max_gpm=float(field.max()),
                       missing_values=0, source_contour_segments=len(lines), status=status,
                       field_west=bounds[0], field_east=bounds[1], field_south=bounds[2], field_north=bounds[3])
            for view, region in boxes.items():
                row[f"{view}_clipped_length_degrees"] = sum(line.difference(region).length for line in lines)
                row[f"{view}_visible_length_degrees"] = sum(line.intersection(region).length for line in lines)
                row[f"{view}_frame_crossings"] = sum(len(line.intersection(region.boundary).geoms)
                                                     if hasattr(line.intersection(region.boundary), "geoms")
                                                     else int(not line.intersection(region.boundary).is_empty)
                                                     for line in lines)
            assert row["new_visible_length_degrees"] + 1e-8 >= row["old_visible_length_degrees"]
            rows.append(row)
    destination = REPO / "results/task6_wpsh/contour_display_audit.csv"
    pd.DataFrame(rows).to_csv(destination, index=False)


def phase_separated_figure(collections, entries, path, title):
    """Identical maps for each phase avoid obscuring nearly coincident contours."""
    fig, axes = plt.subplots(len(entries), 2, figsize=(14, 3.25 * len(entries) + 1.8),
                             subplot_kw={"projection": ccrs.PlateCarree(central_longitude=140)})
    for i, (key, label) in enumerate(entries):
        for j, (phase, phase_label) in enumerate((("pos", "PDO positive"), ("neg", "PDO negative"))):
            draw_panel(axes[i, j], collections[key],
                       f"({chr(97 + i * 2 + j)}) {label} | {phase_label}", phases=(phase,))
    save_figure(fig, path, title,
                "Left: PDO positive | Right: PDO negative | JJA, 500 hPa, 5870 gpm",
                separated=True)


def contour_features(ds):
    """Save native-grid contour vertices; GeoJSON uses longitudes [-180,180]."""
    import contourpy
    features = []
    for phase in ("pos", "neg"):
        field = ds[f"z500_{phase}"]
        cg = contourpy.contour_generator(x=field.lon.values, y=field.lat.values, z=field.values)
        segments = cg.lines(LEVEL)
        # Below-threshold composites legitimately have no WPSH contour.
        for k, segment in enumerate(segments):
            # Split at the dateline rather than creating a line across the world.
            lon = (segment[:, 0] + 180) % 360 - 180
            split = np.flatnonzero(abs(np.diff(lon)) > 180) + 1
            for part in np.split(np.column_stack([lon, segment[:, 1]]), split):
                if len(part) >= 2:
                    features.append(dict(type="Feature", properties=dict(dataset=ds.attrs["dataset"],
                                         phase=phase, height_gpm=LEVEL, pressure_hPa=500,
                                         segment=k, period=ds.attrs["requested_period"]),
                                         geometry=dict(type="LineString", coordinates=part.tolist())))
    return features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--redraw", action="store_true", help="Redraw from the submitted composite NetCDF files only")
    args = parser.parse_args()
    out, figures = REPO / "results/task6_wpsh", REPO / "figures/task6_wpsh"
    out.mkdir(parents=True, exist_ok=True)
    if args.redraw:
        collections = {}
        for path in out.glob("*.nc"):
            with xr.open_dataset(path) as ds:
                collections[path.stem] = ds.load()
        render(collections, figures)
        return
    if args.data_root is None:
        parser.error("--data-root is required unless --redraw is used")
    records = inventory(args.data_root)
    obs_path = REPO / "data/observations/ncep_hgt500_monthly_1948_2014.nc"
    obs_provenance = download_observations(obs_path)
    collections, phases_table = {}, []
    phase_path = REPO / "results/task2_obs.nc"
    phases, digest = phase_series(phase_path)
    collections["obs"] = composite(observation_jja(obs_path), phases, "obs", phase_path, digest)
    for experiment in EXPERIMENTS:
        for member in MEMBERS:
            key = f"{experiment}_{member}"
            phase_path = REPO / "results" / (f"task2_{member}.nc" if experiment == "historical"
                                               else f"task4_{experiment}_{member}.nc")
            phases, digest = phase_series(phase_path)
            z = model_jja(args.data_root, experiment, member,
                          REPO / f"data/task6_cache/{key}_jja500.nc")
            collections[key] = composite(z, phases, key, phase_path, digest)
            if experiment == "historical":
                common = f"historical_common_{member}"
                collections[common] = composite(z.sel(year=slice(1948, 2014)), phases, common, phase_path, digest)
    for experiment in EXPERIMENTS:
        key = f"{experiment}_member_mean"
        collections[key] = ensemble([collections[f"{experiment}_{m}"] for m in MEMBERS], key)
    features, summary = [], {}
    for key, ds in collections.items():
        encoding = {name: dict(zlib=True, complevel=4) for name in ds.data_vars}
        ds.to_netcdf(out / f"{key}.nc", encoding=encoding)
        summary[key] = dict(ds.attrs)
        summary[key]["composite_ranges_gpm"] = {
            phase: [float(ds[f"z500_{phase}"].min()), float(ds[f"z500_{phase}"].max())]
            for phase in ("pos", "neg")}
        summary[key]["contour_exists_in_computed_domain"] = {
            phase: bool(float(ds[f"z500_{phase}"].min()) < LEVEL < float(ds[f"z500_{phase}"].max()))
            for phase in ("pos", "neg")}
        features.extend(contour_features(ds))
        if "year" in ds.dims:
            for year, phase, value in zip(ds.year.values, ds.phase.values, ds.pdo_jja_smooth.values):
                phases_table.append(dict(dataset=key, year=int(year), pdo_phase=int(phase),
                                         pdo_jja_smooth=float(value), included=bool(phase)))
        print(f"Composite {key}: {ds.attrs.get('n_pos','member mean')} / {ds.attrs.get('n_neg','member mean')}", flush=True)
    pd.DataFrame(phases_table).to_csv(out / "phase_years.csv", index=False)
    (out / "contours_5870gpm.geojson").write_text(json.dumps(dict(type="FeatureCollection", features=features)), encoding="utf-8")
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    manifest = dict(model_data_root=str(args.data_root.resolve()), model_files=records,
                    observations=obs_provenance, repository_base_commit="7fb16afb1522f7f96ed32987a0ba3d46648b52ec",
                    input_file_count=14, files_used_for_composites=12,
                    excluded_files=2, selected_members=list(MEMBERS),
                    phase_method="Unchanged repository PDO grouping; all ENSO states retained",
                    model_height_conversion="zg in m already is geopotential height; no division by gravity",
                    observation_height_conversion="hgt in m already is geopotential height",
                    versions={"python": platform.python_version(), "numpy": np.__version__,
                              "xarray": xr.__version__, "matplotlib": matplotlib.__version__,
                              "netCDF4": netCDF4.__version__})
    (out / "input_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    render(collections, figures)
    print("Completed figures and reproducible composite data.", flush=True)


def render(collections, figures):
    contour_display_audit(collections)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                         "svg.hashsalt": "Sandro-JJA-WPSH-5870"})
    projection = ccrs.PlateCarree(central_longitude=140)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), subplot_kw={"projection": projection})
    entries = [("obs", "(a) Observation reference: NCEP/NCAR"),
               ("historical_member_mean", "(b) ACCESS-CM2 historical"),
               ("ssp245_member_mean", "(c) ACCESS-CM2 SSP245"),
               ("ssp585_member_mean", "(d) ACCESS-CM2 SSP585")]
    for ax, (key, title) in zip(axes.flat, entries):
        draw_panel(ax, collections[key], title)
    save_figure(fig, figures / "Figure4_JJA_WPSH_overview",
                "JJA western Pacific subtropical high", "500 hPa | 5870 gpm | PDO phases from the repository's 9-year smoothed PC1")
    fig, axes = plt.subplots(3, 3, figsize=(18, 11.6), subplot_kw={"projection": projection})
    for i, experiment in enumerate(EXPERIMENTS):
        for j, member in enumerate(MEMBERS):
            draw_panel(axes[i, j], collections[f"{experiment}_{member}"],
                       f"({chr(97+i*3+j)}) {experiment.upper()}  {member}i1p1f1")
    save_figure(fig, figures / "Figure4_JJA_WPSH_selected_members",
                "JJA WPSH position by selected ACCESS-CM2 member", "Rows: historical, SSP245, SSP585 | Columns: r1, r2, r5 | 500 hPa / 5870 gpm")
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), subplot_kw={"projection": projection})
    draw_panel(axes[0, 0], collections["obs"], "(a) Observation reference: NCEP/NCAR")
    for i, (ax, member) in enumerate(zip(axes.flat[1:], MEMBERS)):
        draw_panel(ax, collections[f"historical_common_{member}"],
                   f"({chr(98+i)}) ACCESS-CM2 historical {member}i1p1f1")
    save_figure(fig, figures / "Figure4_JJA_WPSH_common_period",
                "JJA WPSH: observation and historical members", "Common height-data window: 1948-2014 | Original repository PDO indices and phase definitions retained")
    phase_separated_figure(collections,
                           [("obs", "Observation: NCEP/NCAR")] +
                           [(f"historical_{m}", f"Historical {m}i1p1f1") for m in MEMBERS],
                           figures / "Figure4_JJA_WPSH_observation_historical_phase_separated",
                           "JJA WPSH: observation and historical phases")
    for experiment in ("ssp245", "ssp585"):
        phase_separated_figure(collections,
                               [(f"{experiment}_{m}", f"{experiment.upper()} {m}i1p1f1") for m in MEMBERS],
                               figures / f"Figure4_JJA_WPSH_{experiment}_phase_separated",
                               f"JJA WPSH: {experiment.upper()} phases by member")
    print("Rendered six PNG/SVG figures and the contour display audit.", flush=True)


if __name__ == "__main__":
    main()
