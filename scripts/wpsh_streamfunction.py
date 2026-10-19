"""Global wind -> streamfunction -> detrended JJA covariance EOFs.

The streamfunction operator is linear. Temporal averaging and polynomial
coefficients can therefore be formed in wind space before the expensive global
spherical-harmonic transforms, with the same result as transforming every month.
The tests compare both orders explicitly.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import xarray as xr

RADIUS = 6_371_200.0
DOMAIN = (110, 180, 10, 40)
CORE = (120, 150, 15, 25)
HISTORICAL_PERIOD = (1900, 2014)
FUTURE_PERIOD = (2015, 2100)


def compare_periods(pc, historical=HISTORICAL_PERIOD, future=FUTURE_PERIOD):
    """Compare annual indices using a fixed 1.5-SD threshold and annual rates.

    The input is already standardized with the intended reference period.
    Rates are percentages of summers, so unequal period lengths are comparable.
    """
    if historical[1] >= future[0]:
        raise ValueError("Historical and future periods must not overlap")
    result = {}
    for label, (first, last) in (("historical", historical), ("future", future)):
        subset = pc.sel(year=slice(first, last))
        if not np.array_equal(subset.year.values, np.arange(first, last+1)):
            raise ValueError(f"Incomplete annual coverage for {label}")
        if not np.isfinite(subset.values).all():
            raise ValueError(f"Nonfinite index in {label}")
        n = subset.sizes["year"]
        count = int((subset > 1.5).sum())
        result.update({label: f"{first}-{last}", f"n_years_{label}": n,
                       f"sd_{label}": float(subset.std(ddof=1)),
                       f"strong_events_{label}": count,
                       f"strong_event_rate_{label}_percent": 100*count/n})
    if result["sd_historical"] <= 0:
        raise ValueError("Historical index must have nonzero variability")
    result["future_to_historical_sd_ratio"] = result["sd_future"]/result["sd_historical"]
    result["sd_change_percent"] = (result["future_to_historical_sd_ratio"]-1)*100
    return result


def sha256(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def canonical_wind(ds, variable):
    """Require an explicitly identified 850-hPa eastward/northward wind."""
    da = ds[variable]
    expected = "eastward_wind" if variable in ("u", "ua", "uwnd") else "northward_wind"
    if da.attrs.get("standard_name", expected) != expected:
        raise ValueError(f"{variable} is not the expected {expected}")
    rename = {old: new for old, new in [("valid_time", "time"), ("latitude", "lat"),
                                       ("longitude", "lon")] if old in da.dims}
    da = da.rename(rename)
    level_names = [x for x in ("plev", "level", "pressure_level", "isobaricInhPa")
                   if x in da.coords]
    if len(level_names) != 1:
        raise ValueError("Wind must have exactly one explicit pressure coordinate")
    level = level_names[0]
    units = str(da[level].attrs.get("units", "")).lower()
    if units in ("pa", "pascal", "pascals"):
        pressure = np.atleast_1d(da[level].values).astype(float) / 100
    elif units in ("hpa", "millibar", "millibars", "mb", "mbar"):
        pressure = np.atleast_1d(da[level].values).astype(float)
    else:
        raise ValueError(f"Unknown pressure units: {units!r}")
    match = np.flatnonzero(np.isclose(pressure, 850, rtol=0, atol=0.01))
    if len(match) != 1:
        raise ValueError("An exact 850-hPa pressure level is required")
    if level in da.dims:
        da = da.isel({level: int(match[0])}, drop=True)
    else:
        da = da.drop_vars(level)
    if da.attrs.get("units", "").replace(" ", "") not in ("m/s", "ms-1", "ms**-1", "ms^-1"):
        raise ValueError(f"Expected wind in m/s, got {da.attrs.get('units')}")
    for dim in list(da.dims):
        if dim not in ("time", "lat", "lon"):
            if da.sizes[dim] != 1:
                raise ValueError(f"Select {dim} explicitly before running")
            da = da.isel({dim: 0}, drop=True)
    if set(da.dims) != {"time", "lat", "lon"}:
        raise ValueError("Only rectilinear time/latitude/longitude winds are supported")
    da = da.assign_coords(lon=da.lon % 360).sortby("lon").sortby("lat", ascending=False)
    if len(np.unique(da.lon)) != da.sizes["lon"]:
        raise ValueError("Duplicate cyclic longitude; remove the duplicate endpoint")
    return da.transpose("time", "lat", "lon")


def validate_global_grid(lat, lon):
    lat, lon = np.asarray(lat), np.asarray(lon)
    if len(lat) < 5 or len(lon) < 8:
        raise ValueError("Insufficient global grid")
    spacing = 360 / len(lon)
    if not np.allclose(np.diff(lon), spacing, atol=1e-4, rtol=0):
        raise ValueError("Streamfunction requires complete, equally spaced global longitude")
    try:
        from skyborn.windspharm._common import inspect_gridtype
    except ImportError:
        from windspharm._common import inspect_gridtype
    return inspect_gridtype(lat)


def streamfunction(u, v, lat, lon):
    """u=-1/a dpsi/dphi, v=1/(a cos(phi)) dpsi/dlambda; global zero gauge."""
    if u.shape != v.shape or u.ndim != 3:
        raise ValueError("u and v must have matching (field, lat, lon) shapes")
    if not np.isfinite(u).all() or not np.isfinite(v).all():
        raise ValueError("Global winds contain missing/nonfinite values; no zero filling is allowed")
    gridtype = validate_global_grid(lat, lon)
    try:
        from skyborn.windspharm.standard import VectorWind
    except ImportError:
        from windspharm.standard import VectorWind
    w = VectorWind(np.moveaxis(u, 0, -1), np.moveaxis(v, 0, -1),
                   gridtype=gridtype, rsphere=RADIUS, legfunc="computed")
    return np.moveaxis(w.streamfunction(), -1, 0)


def timeline(time, first, last):
    year, month = time.dt.year.values, time.dt.month.values
    expected = np.arange(first * 12 + 1, (last + 1) * 12 + 1)
    if not np.array_equal(year * 12 + month, expected):
        raise ValueError("Require every month exactly once, in order, for the requested full years")
    return year, month


def temporal_weights(time, first, last, baseline, degree):
    year, month = timeline(time, first, last)
    ref = (year >= baseline[0]) & (year <= baseline[1])
    if baseline[0] < first or baseline[1] > last or baseline[0] > baseline[1]:
        raise ValueError("Climatological baseline must lie inside the analysis period")
    if last-first+1 <= degree+2:
        raise ValueError("Insufficient years for detrending and EOFs")
    years = np.arange(first, last+1)
    design = np.polynomial.polynomial.polyvander(np.linspace(-1, 1, len(time)), degree)
    # Least-squares coefficients of monthly anomalies, fitted over ALL months.
    pinv = np.linalg.pinv(design)
    anomaly_pinv = pinv.copy()
    for m in range(1, 13):
        selected = ref & (month == m)
        anomaly_pinv[:, selected] -= pinv[:, month == m].sum(axis=1)[:, None] / selected.sum()
    weights = np.zeros((len(years)+1+degree+1, len(time)))
    jja_design = []
    for i, y in enumerate(years):
        take = (year == y) & np.isin(month, [6, 7, 8])
        weights[i, take] = 1/3
        jja_design.append(design[take].mean(axis=0))
    ref_jja = ref & np.isin(month, [6, 7, 8])
    weights[len(years), ref_jja] = 1/ref_jja.sum()
    weights[len(years)+1:] = anomaly_pinv
    return weights, np.array(jja_design)


def aggregate_winds(u, v, first, last, baseline, degree, batch=12, wind_fill="reject"):
    u, v = xr.align(u, v, join="exact")
    if not np.array_equal(u.time.dt.year, v.time.dt.year) or not np.array_equal(u.time.dt.month, v.time.dt.month):
        raise ValueError("u/v monthly dates differ")
    weights, design = temporal_weights(u.time, first, last, baseline, degree)
    validate_global_grid(u.lat.values, u.lon.values)
    shape = (weights.shape[0], u.sizes["lat"], u.sizes["lon"])
    ua, va = np.zeros(shape), np.zeros(shape)
    filler = None
    if wind_fill != "reject":
        from wpsh_wind_fill import WindFiller
        filler = WindFiller(u.lat.values, u.lon.values)
    missing_counts = []
    for start in range(0, u.sizes["time"], batch):
        end = min(start+batch, u.sizes["time"])
        w = weights[:, start:end]
        rows = np.flatnonzero(np.any(w != 0, axis=1))
        ub = u.isel(time=slice(start, end)).values.astype("float64")
        vb = v.isel(time=slice(start, end)).values.astype("float64")
        missing_counts.extend(np.isnan(ub).mean(axis=(1, 2)).tolist())
        if filler is not None:
            for k in range(len(ub)):
                ub[k], vb[k] = filler.fill(ub[k], vb[k], wind_fill)
        for values, dest in ((ub, ua), (vb, va)):
            if not np.isfinite(values).all():
                raise ValueError(f"Missing global winds in months {start}:{end}")
            dest[rows] += np.einsum("ft,tij->fij", w[rows], values, optimize=True)
        if start % 120 == 0:
            print(f"  Read monthly winds {start+1}-{end}/{u.sizes['time']}", flush=True)
    return ua, va, design, missing_counts


def regional_products(u, v, first, last, baseline, degree, batch=4, wind_fill="reject"):
    ua, va, design, missing_counts = aggregate_winds(u, v, first, last, baseline, degree, wind_fill=wind_fill)
    lat, lon = u.lat.values, u.lon.values
    # Keep a wider map than the EOF box for regression diagnostics.
    # Retain bracketing nodes so a subsequent 1-degree interpolation has no
    # missing eastern edge on ACCESS's 1.875-degree longitude grid.
    south, north = lat[lat <= 0].max(), lat[lat >= 50].min()
    west, east = lon[lon <= 90].max(), lon[lon >= 200].min()
    iy = np.flatnonzero((lat >= south) & (lat <= north))
    ix = np.flatnonzero((lon >= west) & (lon <= east))
    fields = []
    for start in range(0, len(ua), batch):
        stop = min(start+batch, len(ua))
        sf = streamfunction(ua[start:stop], va[start:stop], lat, lon)
        fields.append(sf[:, iy][:, :, ix].astype("float64"))
        print(f"  Global streamfunction {start+1}-{stop}/{len(ua)}", flush=True)
    fields = np.concatenate(fields)
    years = np.arange(first, last+1)
    n = len(years)
    raw, climatology, coefficients = fields[:n], fields[n], fields[n+1:]
    anomaly = raw-climatology
    monthly_residual = anomaly-np.einsum("yk,kij->yij", design, coefficients)
    seasonal_design = np.polynomial.polynomial.polyvander(np.linspace(-1, 1, n), degree)
    flat = anomaly.reshape(n, -1)
    seasonal_residual = (flat-seasonal_design @ np.linalg.lstsq(seasonal_design, flat, rcond=None)[0]).reshape(raw.shape)
    dims = ("year", "lat", "lon")
    ds = xr.Dataset(dict(sf_jja=(dims, raw), sf_climatology=(("lat", "lon"), climatology),
                         sf_anomaly=(dims, anomaly),
                         sf_detrended_monthly=(dims, monthly_residual),
                         sf_detrended_seasonal=(dims, seasonal_residual)),
                    coords=dict(year=years, lat=lat[iy], lon=lon[ix]))
    for name in ds.data_vars:
        ds[name].attrs["units"] = "m2 s-1"
    ds.attrs.update(pressure_hPa=850, source_nlat=len(lat), source_nlon=len(lon),
                    wind_fill=wind_fill,
                    missing_grid_fraction_mean=float(np.mean(missing_counts)),
                    missing_grid_fraction_min=float(np.min(missing_counts)),
                    missing_grid_fraction_max=float(np.max(missing_counts)),
                    baseline_start=baseline[0], baseline_end=baseline[1],
                    detrend_degree=degree, earth_radius_m=RADIUS,
                    solver="SPHEREPACK global spherical harmonics, full resolved truncation, l=0 streamfunction zero",
                    averaging="Equal-weight June, July, August monthly means",
                    streamfunction_convention="u_rot=-dpsi/dlatitude/a; v_rot=dpsi/dlongitude/(a*cos(latitude))")
    return ds.sortby("lat")


def area_average(da, box=CORE):
    west, east, south, north = box
    region = da.sel(lat=slice(south, north), lon=slice(west, east))
    return region.weighted(np.cos(np.deg2rad(region.lat))).mean(("lat", "lon"))


def eof_analysis(field, baseline, weighting="coslat", modes=2):
    west, east, south, north = DOMAIN
    region = field.sel(lat=slice(south, north), lon=slice(west, east))
    if region.sizes["lat"] < 2 or region.sizes["lon"] < 2:
        raise ValueError("EOF box not covered")
    x = region.values.astype(float)
    if not np.isfinite(x).all():
        raise ValueError("EOF field contains missing values")
    x = x-x.mean(axis=0)
    weights = np.ones(x.shape[1:])
    if weighting == "coslat":
        weights *= np.sqrt(np.cos(np.deg2rad(region.lat.values)))[:, None]
    elif weighting != "none":
        raise ValueError("Unknown EOF weighting")
    matrix = (x*weights).reshape(len(x), -1)
    # Time covariance avoids a huge spatial covariance matrix on ERA5's 0.25-degree grid.
    eigenvalues, vectors = np.linalg.eigh(matrix @ matrix.T / (len(x)-1))
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, vectors = eigenvalues[order], vectors[:, order]
    eigenvalues = np.maximum(eigenvalues, 0)
    if eigenvalues[0] <= 0 or eigenvalues[modes-1] < eigenvalues[0]*1e-12:
        raise ValueError("Insufficient nonzero EOF modes")
    pcs = vectors[:, :modes]*np.sqrt(len(x)-1)
    whole = field.values-field.values.mean(axis=0)
    pattern = np.einsum("ym,yij->mij", pcs, whole)/(len(x)-1)
    core = area_average(field).values
    # EOF signs are arbitrary: orient mode 1 toward enhanced positive NW Pacific SF.
    for k in range(modes):
        sign = np.sign(np.dot(pcs[:, k], core-core.mean()))
        if sign == 0:
            raise ValueError("EOF sign cannot be anchored to regional streamfunction")
        pcs[:, k] *= sign
        pattern[k] *= sign
    ref = (field.year.values >= baseline[0]) & (field.year.values <= baseline[1])
    if ref.sum() < 3:
        raise ValueError("At least three reference years required for PC normalization")
    pcs_ref = (pcs-pcs[ref].mean(axis=0))/pcs[ref].std(axis=0, ddof=1)
    regional = (core-core[ref].mean())/core[ref].std(ddof=1)
    ds = xr.Dataset(dict(pc=(("year", "mode"), pcs), pc_reference=(("year", "mode"), pcs_ref),
                         regression=(("mode", "lat", "lon"), pattern),
                         variance_fraction=("mode", eigenvalues[:modes]/eigenvalues.sum()),
                         regional_sf=("year", core), regional_sf_reference=("year", regional)),
                    coords=dict(year=field.year, mode=np.arange(1, modes+1), lat=field.lat, lon=field.lon))
    ds.pc.attrs.update(units="1", normalization="Sample SD=1 (ddof=1) over full analysis period")
    ds.pc_reference.attrs.update(units="1", normalization=f"Mean=0, sample SD=1 over {baseline[0]}-{baseline[1]}")
    ds.regression.attrs.update(units="m2 s-1", definition="Regression of detrended JJA SF onto full-period unit-SD PC")
    ds.regional_sf.attrs.update(units="m2 s-1", domain="15-25N, 120-150E, cosine-latitude weighted")
    ds.attrs.update(eof_domain="10-40N, 110-180E", eof_weighting=weighting,
                    eof_type="Covariance EOF; no grid-point variance normalization; temporal mean removed",
                    pc1_sign="Positive covariance with regional 850-hPa streamfunction (15-25N,120-150E)")
    return ds
