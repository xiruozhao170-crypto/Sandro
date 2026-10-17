"""Independent numerical audits of submitted WPSH products and raw ERA5 winds."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import detrend
import xarray as xr

from wpsh_streamfunction import canonical_wind, sha256, streamfunction, timeline

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/task7_wpsh_eof"


def validate_access_downloads():
    directory = ROOT / "data/wpsh_winds/access"
    records = json.loads((directory / "download_manifest.json").read_text())
    if len(records) != 102 or len({r["file"] for r in records}) != 102:
        raise AssertionError("Expected 102 distinct complete ACCESS subsets")
    coverage = {}
    total_bytes = 0
    pressure_errors = []
    for record in records:
        path = directory / record["file"]
        assert sha256(path) == record["sha256"], f"Checksum mismatch: {path.name}"
        with xr.open_dataset(path) as ds:
            pressure = float(ds.plev.item())
            # Original CMOR coordinates can encode 85000 as 85000.00000001.
            assert ds.plev.attrs["units"] == "Pa" and np.isclose(pressure, 85000, rtol=0, atol=1e-4)
            pressure_errors.append(abs(pressure-85000))
            np.testing.assert_allclose(ds.lat, np.linspace(-90, 90, 145), atol=1e-8)
            np.testing.assert_allclose(ds.lon, np.arange(192)*1.875, atol=1e-8)
            years, months = timeline(ds.time, record["first_year"], record["last_year"])
            assert ds.attrs["source_url"] == record["opendap_url"]
            key = "_".join(path.stem.split("_")[:3])
            coverage.setdefault(key, []).extend((years*12+months).tolist())
        total_bytes += path.stat().st_size
    groups = []
    for key, dates in sorted(coverage.items()):
        historical = "_historical_" in key
        first, last = (1900, 2014) if historical else (2015, 2099)
        np.testing.assert_array_equal(sorted(dates), np.arange(first*12+1, (last+1)*12+1))
        groups.append(dict(group=key, first_year=first, last_year=last, months=len(dates)))
    assert len(groups) == 18
    report = dict(subsets=102, total_bytes=total_bytes, sha256_verified=True,
                  nominal_pressure_Pa=85000, pressure_max_abs_deviation_Pa=max(pressure_errors),
                  global_grid=[145, 192],
                  all_months_exactly_once=True, groups=groups)
    (RESULTS / "access_download_validation.json").write_text(json.dumps(report, indent=2)+"\n")
    return report


def validate(path, raw_check=False):
    with xr.open_dataset(path) as opened:
        ds = opened.load()
    for name in ds.data_vars:
        if not np.isfinite(ds[name].values).all():
            raise AssertionError(f"Nonfinite result: {path.name}/{name}")
    np.testing.assert_allclose(ds.pc.mean("year"), 0, atol=1e-12)
    np.testing.assert_allclose(ds.pc.std("year", ddof=1), 1, atol=1e-12)
    np.testing.assert_allclose(np.cov(ds.pc.values.T), np.eye(2), atol=1e-12)
    ref = ds.pc_reference.sel(year=slice(ds.attrs["baseline_start"], ds.attrs["baseline_end"]))
    np.testing.assert_allclose(ref.mean("year"), 0, atol=1e-12)
    np.testing.assert_allclose(ref.std("year", ddof=1), 1, atol=1e-12)
    np.testing.assert_array_equal(ds.strong_event, ds.pc_reference.sel(mode=1) > 1.5)
    csv = pd.read_csv(path.with_name(path.stem+"_index.csv"))
    np.testing.assert_array_equal(csv.year, ds.year)
    np.testing.assert_allclose(csv.pc1, ds.pc.sel(mode=1), atol=1e-14)
    np.testing.assert_allclose(csv.pc2, ds.pc.sel(mode=2), atol=1e-14)
    np.testing.assert_allclose(csv.pc1_reference, ds.pc_reference.sel(mode=1), atol=1e-14)
    np.testing.assert_array_equal(csv.strong_event, ds.strong_event)
    summary = json.loads(path.with_name(path.stem+"_summary.json").read_text())
    np.testing.assert_allclose(summary["eof1_variance_percent"], float(ds.variance_fraction.sel(mode=1))*100)
    assert summary["strong_years_reference_sd"] == ds.year.values[ds.strong_event.values.astype(bool)].tolist()
    if "century_comparison" in summary:
        pc = ds.pc_reference.sel(mode=1)
        for label, first, last in (("present", 1900, 1999), ("future", 2000, 2099)):
            subset = pc.sel(year=slice(first, last))
            assert summary["century_comparison"][f"strong_events_{label}"] == int((subset > 1.5).sum())
            np.testing.assert_allclose(summary["century_comparison"][f"sd_{label}"], float(subset.std(ddof=1)))
    if ds.attrs["detrend_degree"] == 1:
        direct = detrend(ds.sf_anomaly.values.astype(float), axis=0)
        np.testing.assert_allclose(ds.sf_detrended_seasonal, direct, atol=3, rtol=2e-6)
    else:
        # Independently refit the stored seasonal anomalies, including degree 2.
        design = np.polynomial.polynomial.polyvander(np.linspace(-1, 1, ds.sizes["year"]), ds.attrs["detrend_degree"])
        flat = ds.sf_anomaly.values.astype(float).reshape(ds.sizes["year"], -1)
        direct = (flat-design @ np.linalg.lstsq(design, flat, rcond=None)[0]).reshape(ds.sf_anomaly.shape)
        np.testing.assert_allclose(ds.sf_detrended_seasonal, direct, atol=3, rtol=2e-6)
    # Independent SVD on the weighted data matrix, rather than the implemented
    # time-covariance eigensolver.
    field = ds[f"sf_detrended_{ds.attrs['detrend_order']}"]
    box = field.sel(lat=slice(10, 40), lon=slice(110, 180))
    centered = box.values.astype(float)-box.values.mean(axis=0, dtype=float)
    weighted = centered*np.sqrt(np.cos(np.deg2rad(box.lat.values)))[None, :, None]
    left, singular, _ = np.linalg.svd(weighted.reshape(len(box.year), -1), full_matrices=False)
    fractions = (singular**2/(singular**2).sum())[:2]
    np.testing.assert_allclose(fractions, ds.variance_fraction, atol=1e-7, rtol=1e-6)
    correlations = [abs(float(np.corrcoef(left[:, i], ds.pc.isel(mode=i))[0, 1])) for i in range(2)]
    if min(correlations) < 1-1e-10:
        raise AssertionError("Independent SVD PC disagreement")
    # Check the complete physical regression maps, including outside the EOF box.
    x = field.values.astype(float)-field.values.mean(axis=0, dtype=float)
    for i in range(2):
        pc = ds.pc.isel(mode=i).values
        regression = np.einsum("y,yij->ij", pc, x)/np.dot(pc, pc)
        np.testing.assert_allclose(regression, ds.regression.isel(mode=i), atol=3, rtol=2e-6)
    audit = dict(file=path.name, n_years=ds.sizes["year"], finite_fields=True,
                 pc_zero_mean_unit_sample_sd=True, reference_pc_normalization=True,
                 csv_matches_netcdf=True, pc_orthogonality=True,
                 independent_svd_pc_correlations=correlations,
                 independent_svd_variance_fractions=fractions.tolist(),
                 physical_regression_maps_verified=True)
    if raw_check and path.stem == "era5_1979_2020":
        raw = ROOT / "data/wpsh_winds/era5_uv850_1979_2020_0.25deg.nc"
        with xr.open_dataset(raw) as winds:
            u = canonical_wind(winds, "u").sel(time=slice("1979-06", "1979-08"))
            v = canonical_wind(winds, "v").sel(time=slice("1979-06", "1979-08"))
            # Transform each raw month FIRST, then average; independent of the
            # optimized annual aggregation used in production.
            sf = streamfunction(u.values.astype(float), v.values.astype(float), u.lat, u.lon)
            actual = xr.DataArray(sf.mean(axis=0), dims=("lat", "lon"), coords=dict(lat=u.lat, lon=u.lon))
            actual = actual.sel(lat=ds.lat, lon=ds.lon)
        expected = ds.sf_jja.sel(year=1979)
        np.testing.assert_allclose(actual, expected, atol=20, rtol=3e-6)
        audit["raw_1979_JJA_monthly_first_max_abs_error_m2_s"] = float(abs(actual-expected).max())
        audit["raw_1979_JJA_relative_rms_error"] = float(np.sqrt(((actual-expected)**2).mean())/np.sqrt((actual**2).mean()))
    if raw_check and path.stem.startswith("access_"):
        from wpsh_wind_fill import WindFiller
        member = path.stem.split("_")[2]
        data = ROOT / "data/wpsh_winds/access"
        with xr.open_dataset(data / f"ua850_historical_{member}_1900_1919.nc") as ud, xr.open_dataset(data / f"va850_historical_{member}_1900_1919.nc") as vd:
            u = canonical_wind(ud, "ua").sel(time=slice("1900-06", "1900-08"))
            v = canonical_wind(vd, "va").sel(time=slice("1900-06", "1900-08"))
            filler = WindFiller(u.lat.values, u.lon.values)
            uv = [filler.fill(u.values[i], v.values[i], ds.attrs["wind_fill"]) for i in range(3)]
            sf = streamfunction(np.stack([a for a, _ in uv]), np.stack([b for _, b in uv]), u.lat, u.lon)
            actual = xr.DataArray(sf.mean(axis=0), dims=("lat", "lon"), coords=dict(lat=u.lat, lon=u.lon))
            actual = actual.sortby("lat").interp(lat=ds.lat, lon=ds.lon)
        expected = ds.sf_jja.sel(year=1900)
        np.testing.assert_allclose(actual, expected, atol=20, rtol=3e-6)
        audit["raw_1900_JJA_monthly_first_max_abs_error_m2_s"] = float(abs(actual-expected).max())
        audit["raw_1900_JJA_relative_rms_error"] = float(np.sqrt(((actual-expected)**2).mean())/np.sqrt((actual**2).mean()))
    return audit


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-check", action="store_true", help="Also check the original 0.25-degree ERA5 winds")
    args = p.parse_args()
    if args.raw_check and list(RESULTS.glob("access_*.nc")):
        validate_access_downloads()
    audits = [validate(path, args.raw_check) for path in sorted(RESULTS.glob("*.nc"))]
    if not audits:
        raise SystemExit("No WPSH result files")
    target = RESULTS / ("validation_raw.json" if args.raw_check else "validation.json")
    target.write_text(json.dumps(audits, indent=2)+"\n")
    if args.raw_check:
        core_audits = [{k: v for k, v in audit.items() if not k.startswith("raw_")} for audit in audits]
        (RESULTS / "validation.json").write_text(json.dumps(core_audits, indent=2)+"\n")
    print(json.dumps(audits, indent=2))
