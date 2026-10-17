"""Physical-coordinate, boundary-preservation and processing-order checks."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from wpsh_wind_fill import WindFiller
from wpsh_streamfunction import regional_products, streamfunction


def example():
    lat, lon = np.linspace(90, -90, 37), np.arange(72)*5.
    phi, lam = np.deg2rad(lat)[:, None], np.deg2rad(lon)[None, :]
    u = 5*np.cos(phi)+np.sin(phi)*np.cos(lam)
    v = -np.broadcast_to(np.sin(lam), u.shape).copy()
    mask = np.zeros(u.shape, bool)
    mask[8:12, :3] = True
    mask[8:12, -3:] = True  # crosses the cyclic seam
    mask[-5:] = True        # includes an entirely missing pole
    return lat, lon, u, v, mask


@pytest.mark.parametrize("method", ["harmonic", "nearest"])
def test_valid_winds_untouched_and_finite_polar_vectors(method):
    lat, lon, u, v, mask = example()
    fill = WindFiller(lat, lon)
    a, b = fill.fill(np.where(mask, np.nan, u), np.where(mask, np.nan, v), method)
    np.testing.assert_array_equal(a[~mask], u[~mask])
    np.testing.assert_array_equal(b[~mask], v[~mask])
    assert np.isfinite(a).all() and np.isfinite(b).all()
    cartesian = a[-1, :, None]*fill.east[-1] + b[-1, :, None]*fill.north[-1]
    np.testing.assert_allclose(cartesian, np.broadcast_to(cartesian[0], cartesian.shape), atol=1e-12)
    with pytest.raises(ValueError, match="entirely missing"):
        fill.fill(np.full_like(u, np.nan), np.full_like(v, np.nan), method)


def test_harmonic_rotation_equivariance_and_fixed_mask_linearity():
    lat, lon, u, v, mask = example()
    fill = WindFiller(lat, lon)
    uw, vw = np.where(mask, np.nan, u), np.where(mask, np.nan, v)
    a, b = fill.fill(uw, vw)
    c, d = fill.fill(np.roll(uw, 13, axis=1), np.roll(vw, 13, axis=1))
    np.testing.assert_allclose(c, np.roll(a, 13, axis=1), atol=1e-12)
    np.testing.assert_allclose(d, np.roll(b, 13, axis=1), atol=1e-12)
    c, d = fill.fill(-2.7*uw, -2.7*vw)
    np.testing.assert_allclose(c, -2.7*a, atol=1e-12)
    np.testing.assert_allclose(d, -2.7*b, atol=1e-12)


def test_partly_missing_pole_and_mismatched_masks():
    lat, lon, u, v, _ = example()
    fill = WindFiller(lat, lon)
    uw, vw = u.copy(), v.copy()
    uw[0, :12] = np.nan
    with pytest.raises(ValueError, match="masks differ"):
        fill.fill(uw, vw)
    vw[0, :12] = np.nan
    for method in ("harmonic", "nearest"):
        a, b = fill.fill(uw, vw, method)
        np.testing.assert_allclose(a, u, atol=1e-12)
        np.testing.assert_allclose(b, v, atol=1e-12)


def test_monthly_changing_masks_filled_before_aggregation():
    lat, lon, u, v, mask = example()
    dates = pd.date_range("2000-01-01", "2004-12-01", freq="MS")
    amplitude = np.sin(np.arange(len(dates))*.21)
    uw, vw = amplitude[:, None, None]*u, amplitude[:, None, None]*v
    for k in range(len(dates)):
        varying = np.roll(mask, k % 4, axis=1)
        uw[k, varying] = np.nan
        vw[k, varying] = np.nan
    coords = dict(time=dates, lat=lat, lon=lon)
    ua = xr.DataArray(uw, dims=("time", "lat", "lon"), coords=coords)
    va = xr.DataArray(vw, dims=ua.dims, coords=coords)
    product = regional_products(ua, va, 2000, 2004, [2000, 2004], 1, wind_fill="harmonic")
    filler = WindFiller(lat, lon)
    filled = [filler.fill(uw[k], vw[k]) for k in range(len(dates))]
    sf = streamfunction(np.stack([x[0] for x in filled]), np.stack([x[1] for x in filled]), lat, lon)
    sf = xr.DataArray(sf, dims=ua.dims, coords=coords)
    direct = sf.where(sf.time.dt.month.isin([6, 7, 8]), drop=True).groupby("time.year").mean()
    direct = direct.sel(lat=product.lat, lon=product.lon)
    np.testing.assert_allclose(product.sf_jja, direct, rtol=2e-5, atol=15)
