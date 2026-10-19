"""Scientific checks against analytic fields and independent direct calculations."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from wpsh_streamfunction import (RADIUS, canonical_wind, eof_analysis,
                                regional_products, streamfunction,
                                temporal_weights, validate_global_grid, compare_periods)


def analytic_grid(nlat=37, nlon=72):
    lat, lon = np.linspace(90, -90, nlat), np.arange(nlon)*360/nlon
    phi, lam = np.deg2rad(lat)[:, None], np.deg2rad(lon)[None, :]
    amplitude = 2e7
    psi = amplitude*(np.sin(phi)+.6*np.cos(phi)*np.cos(lam))
    u = amplitude/RADIUS*(-np.cos(phi)+.6*np.sin(phi)*np.cos(lam))
    v = -amplitude/RADIUS*.6*np.broadcast_to(np.sin(lam), psi.shape)
    return lat, lon, phi, lam, psi, u, v


def test_global_inversion_sign_magnitude_and_divergent_winds():
    lat, lon, phi, lam, psi, u, v = analytic_grid()
    # Add a known velocity-potential gradient: it must not affect psi.
    divergent_u = 1.2*np.broadcast_to(np.cos(lam), u.shape)
    divergent_v = -1.2*np.sin(phi)*np.sin(lam)
    actual = streamfunction((u+divergent_u)[None], (v+divergent_v)[None], lat, lon)[0]
    np.testing.assert_allclose(actual, psi, atol=30, rtol=2e-6)


def test_global_domain_and_missing_winds_are_required():
    lat, lon, _, _, _, u, v = analytic_grid()
    with pytest.raises(ValueError):
        validate_global_grid(lat, lon[:20])
    u[10, 10] = np.nan
    with pytest.raises(ValueError, match="missing/nonfinite"):
        streamfunction(u[None], v[None], lat, lon)


@pytest.mark.parametrize("unit,level", [("Pa", 85000), ("hPa", 850)])
def test_pressure_selection_is_exact_and_handles_coordinate_order(unit, level):
    da = xr.DataArray(np.ones((2, 2, 3, 4)), dims=("time", "plev", "latitude", "longitude"),
                      coords=dict(time=pd.date_range("2000-01-01", periods=2, freq="MS"),
                                  plev=[level, level*.6], latitude=[-90, 0, 90], longitude=[-180, -90, 0, 90]),
                      attrs=dict(units="m s-1"))
    da.plev.attrs["units"] = unit
    out = canonical_wind(da.to_dataset(name="ua"), "ua")
    np.testing.assert_equal(out.lat.values, [90, 0, -90])
    np.testing.assert_equal(out.lon.values, [0, 90, 180, 270])
    with pytest.raises(ValueError, match="exact 850"):
        canonical_wind(da.isel(plev=[1]).to_dataset(name="ua"), "ua")


def test_missing_or_duplicate_months_fail():
    dates = xr.DataArray(pd.date_range("2000-01-01", "2005-12-01", freq="MS"), dims="time")
    with pytest.raises(ValueError, match="every month"):
        temporal_weights(dates.isel(time=slice(1, None)), 2000, 2005, [2000, 2005], 1)
    with pytest.raises(ValueError, match="every month"):
        temporal_weights(xr.concat([dates, dates[-1:]], dim="time"), 2000, 2005, [2000, 2005], 1)


@pytest.mark.parametrize("degree", [1, 2])
def test_optimized_temporal_order_matches_direct_monthly_solver(degree):
    lat, lon, _, _, psi, u, v = analytic_grid(19, 36)
    dates = pd.date_range("2000-01-01", "2007-12-01", freq="MS")
    t = np.linspace(-1, 1, len(dates))
    amp = np.sin(np.arange(len(dates))*.17)+t+2*t*t+np.cos(dates.month.to_numpy()*np.pi/6)
    uw = xr.DataArray(amp[:, None, None]*u, dims=("time", "lat", "lon"),
                      coords=dict(time=dates, lat=lat, lon=lon))
    vw = xr.DataArray(amp[:, None, None]*v, dims=uw.dims, coords=uw.coords)
    result = regional_products(uw, vw, 2000, 2007, [2000, 2003], degree)
    direct = xr.DataArray(streamfunction(uw.values, vw.values, lat, lon), dims=uw.dims, coords=uw.coords)
    clim = direct.sel(time=slice("2000", "2003")).groupby("time.month").mean()
    anom = direct.groupby("time.month")-clim
    design = np.polynomial.polynomial.polyvander(t, degree)
    values = anom.values.reshape(len(t), -1)
    residual = values-design @ np.linalg.lstsq(design, values, rcond=None)[0]
    direct = xr.DataArray(residual.reshape(anom.shape), dims=anom.dims, coords=anom.coords)
    jja = direct.where(direct.time.dt.month.isin([6, 7, 8]), drop=True).groupby("time.year").mean()
    jja = jja.sortby("lat").sel(lat=result.lat, lon=result.lon)
    np.testing.assert_allclose(result.sf_detrended_monthly, jja, rtol=2e-5, atol=30)


def test_covariance_eof_matches_independent_svd():
    rng = np.random.default_rng(20261017)
    years, lat, lon = np.arange(1979, 2021), np.arange(10, 41, 2.5), np.arange(110, 181, 2.5)
    signal = rng.normal(size=len(years))
    values = signal[:, None, None]*1e6+rng.normal(size=(len(years), len(lat), len(lon)))*2e5
    field = xr.DataArray(values, dims=("year", "lat", "lon"), coords=dict(year=years, lat=lat, lon=lon))
    ds = eof_analysis(field, [1979, 1999])
    matrix = ((values-values.mean(axis=0))*np.sqrt(np.cos(np.deg2rad(lat)))[None, :, None]).reshape(len(years), -1)
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    np.testing.assert_allclose(ds.variance_fraction, (s*s/(s*s).sum())[:2], rtol=1e-12)
    assert abs(np.corrcoef(ds.pc.sel(mode=1), u[:, 0])[0, 1]) > .999999999
    assert np.corrcoef(ds.pc.sel(mode=1), signal)[0, 1] > .99
    np.testing.assert_allclose(ds.pc.std("year", ddof=1), 1, atol=1e-12)
    ref = ds.pc_reference.sel(year=slice(1979, 1999))
    np.testing.assert_allclose(ref.mean("year"), 0, atol=1e-12)
    np.testing.assert_allclose(ref.std("year", ddof=1), 1, atol=1e-12)


def test_historical_future_boundary_and_unequal_period_event_rates():
    years = np.arange(1900, 2100)
    values = np.zeros(len(years))
    values[[0, 114, 115, 199]] = 2  # Include 2014 in historical and 2015 in future.
    values[1] = 1.5  # The threshold is strict, not >=.
    pc = xr.DataArray(values, dims="year", coords={"year": years})
    result = compare_periods(pc)
    assert result["n_years_historical"] == 115
    assert result["n_years_future"] == 85
    assert result["strong_events_historical"] == result["strong_events_future"] == 2
    np.testing.assert_allclose(result["strong_event_rate_historical_percent"], 200/115)
    np.testing.assert_allclose(result["strong_event_rate_future_percent"], 200/85)
    expected_ratio = np.std(values[115:], ddof=1)/np.std(values[:115], ddof=1)
    np.testing.assert_allclose(result["sd_change_percent"], (expected_ratio-1)*100)
    with pytest.raises(ValueError, match="Incomplete annual coverage"):
        compare_periods(pc.drop_sel(year=2014))
