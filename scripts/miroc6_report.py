"""Rebuild the MIROC6 report figures from numerical climate data.

The historical method follows the supplied Task 3--4 notebook: linear
detrending in month-index space, monthly anomalies, sqrt(cos(lat)) EOF
weights, unit-standard-deviation PC1, and unweighted SST covariance / PC
variance. Rainfall phases use JJA means of the centered 108-month PC1 mean.
No image manipulation or generative image model is used.

Compute: python scripts/miroc6_report.py --data-root "E:/Sandro Dachuang" --compute
Redraw from the committed numerical results: python scripts/miroc6_report.py
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import warnings

import numpy as np
import xarray as xr
from eofs.xarray import Eof
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/miroc6_report"
FIGURES = ROOT / "figures/task1-4/report_miroc6"
SCENARIOS = ["ssp126", "ssp245", "ssp370", "ssp585"]
MEMBERS = ["r1", "r5", "r9"]
RAIN_KEYS = ("pos_mean", "neg_mean", "composite", "pval")
EOF_LEVELS = np.linspace(-1, 1, 21)
RAIN_LEVELS = np.linspace(-0.8, 0.8, 21)
CRS = ccrs.PlateCarree()


def runtime_versions():
    return {name: importlib.metadata.version(name) for name in
            ["numpy", "scipy", "xarray", "eofs", "matplotlib", "cartopy", "netCDF4", "pillow"]}


def normalize(da):
    rename = {k: v for k, v in [("latitude", "lat"), ("longitude", "lon")]
              if k in da.dims}
    da = da.rename(rename)
    if da.lat.values[0] > da.lat.values[-1]:
        da = da.sortby("lat")
    if float(da.lon.min()) < 0:
        da = da.assign_coords(lon=da.lon % 360).sortby("lon")
    return da


def read_region(path, variable, start, end, ocean=False):
    with xr.open_dataset(path) as ds:
        da = normalize(ds[variable])
        da = da.sel(time=slice(f"{start}-01-01", f"{end}-12-31"),
                    lat=slice(20, 60) if ocean else slice(15, 45),
                    lon=slice(110, 250) if ocean else slice(100, 130)).load()
    assert da.sizes["time"] == (end - start + 1) * 12, path
    return da.where(np.abs(da) < 1e10)


def pdo(sst):
    nt = sst.sizes["time"]
    flat = sst.values.reshape(nt, -1)
    valid = np.all(np.isfinite(flat), axis=0)
    detrended = flat.copy()
    design = np.vstack([np.arange(nt), np.ones(nt)]).T
    coeff = np.linalg.lstsq(design, flat[:, valid], rcond=None)[0]
    detrended[:, valid] = flat[:, valid] - design @ coeff
    da = xr.DataArray(detrended.reshape(sst.shape), coords=sst.coords, dims=sst.dims)
    anom = da.groupby("time.month") - da.groupby("time.month").mean("time")
    solver = Eof(anom * np.sqrt(np.cos(np.deg2rad(anom.lat))), weights=None)
    loading = solver.eofs(neofs=1).squeeze()
    pc = solver.pcs(npcs=1, pcscaling=1).squeeze()
    if float(loading.sel(lat=slice(35, 45), lon=slice(160, 200)).mean()) > 0:
        pc = -pc
    pc = (pc - pc.mean("time")) / pc.std("time")
    # Keep the original notebook's sample covariance / population variance convention.
    pattern = xr.cov(anom, pc, dim="time") / pc.var("time")
    if float(pattern.sel(lat=slice(35, 45), lon=slice(160, 200)).mean()) > 0:
        pattern, pc = -pattern, -pc
    smooth = pc.rolling(time=108, center=True).mean()
    summer = smooth.sel(time=smooth.time.dt.month.isin([6, 7, 8])).groupby("time.year").mean()
    years = pc.time.dt.year.values + (pc.time.dt.month.values - 0.5) / 12
    return dict(eof1=pattern.values, eof1_lats=pattern.lat.values,
                eof1_lons=pattern.lon.values, pc1=pc.values, pc1_years=years,
                pc1_smooth=smooth.values, pc1_smooth_years=years,
                pos_years=summer.year.values[summer.values > 0],
                neg_years=summer.year.values[summer.values < 0],
                var_frac=np.asarray(solver.varianceFraction(neigs=1).values[0]))


def rainfall(pr, result, observation=False):
    if observation:
        rate = pr / pr.time.dt.days_in_month
    else:
        units = pr.attrs.get("units", "").lower().replace(" ", "")
        if units in {"kgm-2s-1", "kg/m2/s", "kgm**-2s**-1", "kgm^-2s^-1"}:
            rate = pr * 86400.0
        elif units in {"mm/day", "mmd-1", "mmday-1"}:
            rate = pr
        else:
            raise ValueError(f"Unsupported precipitation units: {units}")
    summer = rate.sel(time=rate.time.dt.month.isin([6, 7, 8])).groupby("time.year").mean()
    anom = summer - summer.mean("year")
    pos = anom.sel(year=result["pos_years"])
    neg = anom.sel(year=result["neg_years"])
    pos_mean, neg_mean = pos.mean("year"), neg.mean("year")
    # Same all-finite Welch test as the source Task 2 notebook.
    valid = np.all(np.isfinite(pos.values), axis=0) & np.all(np.isfinite(neg.values), axis=0)
    pval = np.full(valid.shape, np.nan)
    pval[valid] = stats.ttest_ind(pos.values[:, valid], neg.values[:, valid],
                                 axis=0, equal_var=False).pvalue
    result.update(pos_mean=pos_mean.values, neg_mean=neg_mean.values,
                  composite=(pos_mean - neg_mean).values, pval=pval,
                  ec_lats=pr.lat.values, ec_lons=pr.lon.values)
    return result


def original_cache(folder, member):
    return {p.stem[len(member)+1:]: np.load(p, allow_pickle=False)
            for p in folder.glob(f"{member}_*.npy")}


def apply_land(result, land, lat, lon):
    assert np.array_equal(result["ec_lats"], lat)
    assert np.array_equal(result["ec_lons"], lon)
    result["land_mask"] = land
    for key in RAIN_KEYS:
        if key in result:
            result[key] = np.where(land, result[key], np.nan)
    return result


def save_case(name, data):
    np.savez_compressed(RESULTS / f"{name}.npz", **data)
    print(f"Saved {name}: EOF1 {float(data['var_frac'])*100:.2f}%; "
          f"PDO+ {len(data['pos_years'])}, PDO- {len(data['neg_years'])}", flush=True)


def load_case(name):
    with np.load(RESULTS / f"{name}.npz", allow_pickle=False) as ds:
        return {key: ds[key] for key in ds.files}


def source_record(path):
    stat = path.stat()
    return dict(path=str(path.resolve()), size_bytes=stat.st_size,
                mtime_ns=stat.st_mtime_ns)


def compute(data_root):
    RESULTS.mkdir(parents=True, exist_ok=True)
    records = []
    gpcc_path = data_root / "precip.mon.total.v2018.nc"
    gpcc = read_region(gpcc_path, "precip", 1900, 2014)
    # Gauge-data coverage defines the land grid; an ocean overlay additionally
    # clips interpolated contour edges precisely at the Natural Earth coastline.
    land = np.isfinite(gpcc.values).any(axis=0)
    lat, lon = gpcc.lat.values, gpcc.lon.values
    records.append(source_record(gpcc_path))
    obs_sst = data_root / "sst.mnmean.v4.nc"
    obs = rainfall(gpcc, pdo(read_region(obs_sst, "sst", 1900, 2014, True)), True)
    assert (len(obs["pos_years"]), len(obs["neg_years"])) == (53, 54)
    save_case("observation", apply_land(obs, land, lat, lon))
    records.append(source_record(obs_sst))
    for i in range(1, 11):
        member = f"r{i}"
        data = original_cache(data_root / "task3_npy_historical", member)
        assert len(data["pc1"]) == 1380
        if member == "r10":
            path = next((data_root / "MIROC6_Historical").glob(f"pr_*_{member}i1p1f1_*.nc"))
            fresh = rainfall(read_region(path, "pr", 1900, 2014), dict(data))
            for key in ["pos_mean", "neg_mean", "composite"]:
                np.testing.assert_array_equal(fresh[key], data[key])
            data["pval"] = fresh["pval"]
            records.append(source_record(path))
        save_case(f"historical_{member}", apply_land(data, land, lat, lon))
    # An independent existing member checks that the raw-file implementation
    # reproduces the audited notebook numerical outputs before new members run.
    check_root = data_root / "MIROC6_SSP126"
    tos = next(check_root.glob("tos_*_r3i1p1f1_*.nc"))
    pr = next(check_root.glob("pr_*_r3i1p1f1_*.nc"))
    check = rainfall(read_region(pr, "pr", 2015, 2100),
                     pdo(read_region(tos, "tos", 2015, 2100, True)))
    expected = original_cache(data_root / "task4_future_npy/ssp126", "r3")
    comparison = {}
    for key in ["eof1", "pc1", "pc1_smooth", "pos_mean", "neg_mean", "composite"]:
        np.testing.assert_allclose(check[key], expected[key], rtol=1e-6, atol=1e-7, equal_nan=True)
        comparison[key] = float(np.nanmax(np.abs(check[key]-expected[key])))
    for key in ["pos_years", "neg_years"]:
        np.testing.assert_array_equal(check[key], expected[key])
    print("Notebook reference case reproduced:", comparison, flush=True)
    for scenario in SCENARIOS:
        for member in MEMBERS:
            folder = data_root / f"MIROC6_{scenario.upper()}"
            tos = next(folder.glob(f"tos_*_{member}i1p1f1_*.nc"))
            pr = next(folder.glob(f"pr_*_{member}i1p1f1_*.nc"))
            print(f"Computing {scenario} {member}", flush=True)
            result = rainfall(read_region(pr, "pr", 2015, 2100),
                              pdo(read_region(tos, "tos", 2015, 2100, True)))
            save_case(f"{scenario}_{member}", apply_land(result, land, lat, lon))
            records.extend([source_record(tos), source_record(pr)])
    cache_sources = []
    for path in sorted((data_root / "task3_npy_historical").glob("*.npy")):
        cache_sources.append(dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    manifest = dict(historical_period=[1900, 2014], future_period=[2015, 2100],
                    future_members=MEMBERS, scenarios=SCENARIOS, observations=["ERSSTv4", "GPCC v2018"],
                    smooth_months=108, phase_definition="sign of JJA mean of centered smoothed PC1",
                    spread="sample variance across r1-r10 only, ddof=1; observation excluded",
                    rainfall_mask="finite GPCC coverage during 1900-2014; ocean covered at Natural Earth 1:50 million coastline",
                    selection="r1, r5, r9 explicitly requested; no skill ranking",
                    runtime_versions=runtime_versions(),
                    historical_cache_sha256=cache_sources, raw_sources=records,
                    notebook_reference_case=dict(case="ssp126_r3", max_absolute_errors=comparison))
    (RESULTS / "input_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")


def map_axis(ax, ocean=False):
    extent = [110, 250, 20, 60] if ocean else [100, 130, 15, 45]
    ax.set_extent(extent, crs=CRS)
    if ocean:
        ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#d4d4d4", zorder=4)
    else:
        ax.set_facecolor("white")
        ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", edgecolor="none", zorder=6)
        ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="white", edgecolor="none", zorder=6)
        ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=.35, linestyle=":", zorder=8)
    ax.coastlines(resolution="50m", linewidth=.55, zorder=8)
    gl = ax.gridlines(draw_labels=True, linewidth=.3, color="gray", alpha=.45, linestyle="--", zorder=9)
    gl.top_labels = gl.right_labels = False
    gl.xlocator = matplotlib.ticker.FixedLocator([120, 160, 200, 240] if ocean else [100, 110, 120, 130])
    gl.ylocator = matplotlib.ticker.FixedLocator([20, 40, 60] if ocean else [15, 25, 35, 45])
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = gl.ylabel_style = {"size": 8}


def field(ax, data, key, levels, ocean=False, dots=False, cmap=None):
    la, lo = ("eof1_lats", "eof1_lons") if ocean else ("ec_lats", "ec_lons")
    values = data[key] if ocean else np.where(data["land_mask"], data[key], np.nan)
    cf = ax.contourf(data[lo], data[la], values, levels=levels,
                     cmap=cmap or ("RdBu_r" if ocean else "BrBG"),
                     transform=CRS, extend="both" if levels[0] < 0 else "max")
    if dots and "pval" in data:
        xx, yy = np.meshgrid(data[lo], data[la])
        sig = (data["pval"] < .1) & data["land_mask"] & np.isfinite(values)
        ax.scatter(xx[sig], yy[sig], s=4, color="white", alpha=.85, transform=CRS, zorder=5)
    map_axis(ax, ocean)
    return cf


def figure_path(index):
    matches = list(FIGURES.glob(f"MIROC6_{index:02d}_*.png"))
    if len(matches) != 1:
        raise ValueError(f"Expected one existing filename for {index}: {matches}")
    return matches[0]


def save_figure(fig, index):
    path = figure_path(index)
    geo_axes = [ax for ax in fig.axes if isinstance(ax, GeoAxes)]
    pc_axes = [ax for ax in fig.axes if ax.get_ylabel() == "PC1"]
    entry = dict(index=index, path=path.relative_to(ROOT).as_posix(),
                 map_panels=len(geo_axes), pc_panels=len(pc_axes),
                 panel_titles=[ax.get_title() for ax in geo_axes+pc_axes])
    if index in [3,4,5,6,7,8,9,10,12,15,16,19,20,22]:
        entry["ocean_overlay_on_every_map"] = all(
            any(getattr(getattr(artist, "_feature", None), "name", None) == "ocean"
                for artist in ax.get_children()) for ax in geo_axes)
        assert entry["ocean_overlay_on_every_map"]
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    entry["image_size"] = list(Image.open(path).size)
    manifest_path = RESULTS / "figure_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest[str(index)] = entry
    manifest_path.write_text(json.dumps(manifest, indent=2)+"\n",encoding="utf-8")
    print(f"Rendered {path.name}", flush=True)


def positive_levels(array):
    values = array[np.isfinite(array)]
    return np.linspace(0, max(float(np.percentile(values, 98)), 1e-8), 16)


def variance(data, key):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanvar(np.stack([d[key] for d in data]), axis=0, ddof=1)


def colorbar(fig, cf, rect, label, ticks=None, scientific=False):
    bar = fig.colorbar(cf, cax=fig.add_axes(rect), orientation="horizontal", ticks=ticks)
    bar.set_label(label, fontsize=10)
    bar.ax.tick_params(labelsize=9)
    if scientific:
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 3))
        bar.formatter = formatter
        bar.update_ticks()
    return bar


def plot_historical(hist, obs):
    for key, index, ocean in [("eof1", 13, True), ("composite", 15, False)]:
        fig, axes = plt.subplots(4, 3, figsize=(16, 14 if ocean else 17),
                                 subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180) if ocean else CRS})
        flat = axes.ravel()
        levels = EOF_LEVELS if ocean else RAIN_LEVELS
        for i, d in enumerate(hist):
            cf = field(flat[i], d, key, levels, ocean)
            flat[i].set_title(f"r{i+1} | EOF1 {float(d['var_frac'])*100:.1f}%" if ocean else f"r{i+1} | PDO+ minus PDO-", fontsize=11)
        spread = dict(hist[0], **{key:variance(hist, key)})
        spread_levels = positive_levels(spread[key])
        sp = field(flat[10], spread, key, spread_levels, ocean, cmap="YlOrRd")
        flat[10].set_title("Spread: variance of r1-r10", fontsize=11)
        field(flat[11], obs, key, levels, ocean)
        flat[11].set_title("Observation: ERSSTv4 EOF1" if ocean else "Observation: ERSSTv4 + GPCC", fontsize=11, color="#174a70", fontweight="bold")
        fig.suptitle("Task 3 | Historical PDO EOF1 patterns, observation and member variance\n1900-2014" if ocean else "Task 3 | Historical JJA rainfall composites, observation and member variance\n1900-2014 | Land only", fontsize=16, y=.97)
        fig.subplots_adjust(left=.045, right=.97, top=.91, bottom=.16, wspace=.15, hspace=.3)
        colorbar(fig, cf, [.16,.095,.68,.018], "SST regression/covariance pattern (degC per normalized PC1)" if ocean else "JJA rainfall: PDO+ minus PDO- (mm/day)", np.linspace(levels[0],levels[-1],5))
        colorbar(fig, sp, [.16,.035,.68,.018], "Member variance of SST pattern (degC squared)" if ocean else "Member variance of rainfall composite ((mm/day) squared)", np.linspace(0,spread_levels[-1],5), True)
        save_figure(fig, index)
    plot_historical_pc(hist, obs)
    fig, axes = plt.subplots(1, 2, figsize=(12, 7), subplot_kw={"projection":CRS})
    arrays = [variance(hist,key) for key in ["pos_mean","neg_mean"]]
    levels = positive_levels(np.stack(arrays))
    for ax,key,arr,title in zip(axes,["pos_mean","neg_mean"],arrays,["PDO positive phase","PDO negative phase"]):
        cf = field(ax,dict(hist[0],**{key:arr}),key,levels,cmap="YlOrRd")
        ax.set_title(title, fontsize=12)
    fig.suptitle("Task 3 | Inter-member variance of JJA rainfall anomalies\nMIROC6 r1-r10 | 1900-2014 | Land only", fontsize=14)
    fig.subplots_adjust(left=.06,right=.97,top=.85,bottom=.18,wspace=.18)
    colorbar(fig,cf,[.2,.08,.6,.025],"Sample variance ((mm/day) squared)",scientific=True)
    save_figure(fig,16)


def plot_historical_pc(hist, obs):
    fig = plt.figure(figsize=(16, 16))
    gs = fig.add_gridspec(6, 2, left=.055, right=.98, bottom=.05, top=.93, hspace=.44, wspace=.12)
    for i, d in enumerate(hist):
        ax = fig.add_subplot(gs[i//2, i%2])
        pc_axis(ax, d, f"r{i+1}", 1900, 2014, xlabel=i>=8)
    ax = fig.add_subplot(gs[5, :])
    pc_axis(ax, obs, "Observation: ERSSTv4", 1900, 2014, xlabel=True)
    ax.set_title("Observation: ERSSTv4", color="#174a70", fontweight="bold", fontsize=12)
    fig.suptitle("Task 3 | Historical monthly PDO PC1 and 9-year running mean\nMIROC6 r1-r10 and ERSSTv4 observation | 1900-2014", fontsize=16, y=.98)
    save_figure(fig, 14)


def pc_axis(ax, data, title, start, end, xlabel=False):
    year, pc = data["pc1_years"], data["pc1"]
    ax.fill_between(year, pc, 0, where=pc>=0, color="#db3d3d", alpha=.65, linewidth=0)
    ax.fill_between(year, pc, 0, where=pc<0, color="#356ac3", alpha=.65, linewidth=0)
    ax.plot(year, data["pc1_smooth"], color="black", lw=1.3, label="9-year running mean")
    ax.axhline(0,color="gray",lw=.6)
    # Decimal-year coordinates include December of the final year.
    ax.set(xlim=(start,end+1),ylim=(-4,4),ylabel="PC1",title=title)
    ax.grid(alpha=.25, lw=.4)
    if xlabel:
        ax.set_xlabel("Year")
    ax.legend(loc="upper right",fontsize=7,frameon=False)


def plot_future():
    cases = {(s,m):load_case(f"{s}_{m}") for s in SCENARIOS for m in MEMBERS}
    for index,key,title in [(17,"eof1","PDO EOF1 patterns"),(19,"pos_mean","PDO positive-phase JJA rainfall anomalies"),(20,"neg_mean","PDO negative-phase JJA rainfall anomalies"),(22,"composite","JJA rainfall composites: PDO+ minus PDO-")]:
        ocean = key == "eof1"
        fig, axes = plt.subplots(4,3,figsize=(16,14 if ocean else 17),subplot_kw={"projection":ccrs.PlateCarree(central_longitude=180) if ocean else CRS})
        levels = EOF_LEVELS if ocean else RAIN_LEVELS
        for row,s in enumerate(SCENARIOS):
            for col,m in enumerate(MEMBERS):
                d=cases[s,m]
                cf=field(axes[row,col],d,key,levels,ocean)
                axes[row,col].set_title(f"{s.upper()} | {m}" + (f" | {float(d['var_frac'])*100:.1f}%" if ocean else ""),fontsize=11)
        fig.suptitle(f"Task 4 | Future {title}\nMIROC6 r1, r5, r9 | 2015-2100"+("" if ocean else " | Land only"),fontsize=16,y=.97)
        fig.subplots_adjust(left=.045,right=.97,top=.91,bottom=.10,wspace=.15,hspace=.28)
        colorbar(fig,cf,[.16,.045,.68,.02],"SST regression/covariance pattern (degC per normalized PC1)" if ocean else "JJA rainfall anomaly / composite (mm/day)",np.linspace(levels[0],levels[-1],5))
        save_figure(fig,index)
    plot_future_pc(cases)


def plot_future_pc(cases):
    fig, axes=plt.subplots(4,3,figsize=(16,12))
    for row,s in enumerate(SCENARIOS):
        for col,m in enumerate(MEMBERS):
            pc_axis(axes[row,col],cases[s,m],f"{s.upper()} | {m}",2015,2100,row==3)
    fig.suptitle("Task 4 | Future monthly PDO PC1 and 9-year running mean\nMIROC6 r1, r5, r9 | 2015-2100",fontsize=16)
    fig.tight_layout(rect=[0,0,1,.94],h_pad=2)
    save_figure(fig,18)


def plot_task2(obs, model):
    # Preserve the source notebook's color ranges while removing marine values.
    for d,start,label,phase_max,comp_max in [(obs,3,"Observation: ERSSTv4 + GPCC",.4,.7),(model,7,"MIROC6 historical r10",.5,1.)]:
        for offset,key,title,vmax in [(0,"pos_mean","PDO positive phase",phase_max),(1,"neg_mean","PDO negative phase",phase_max),(2,"composite","PDO+ minus PDO-",comp_max)]:
            fig,ax=plt.subplots(figsize=(7,7),subplot_kw={"projection":CRS})
            cf=field(ax,d,key,np.linspace(-vmax,vmax,21),dots=offset==2)
            ax.set_title(f"{label}\n{title} | JJA 1900-2014 | Land only",fontsize=12)
            cb=fig.colorbar(cf,ax=ax,shrink=.8,pad=.08)
            cb.set_label("mm/day")
            if offset==2:
                fig.text(.5,.025,"White dots: Welch p < 0.1",ha="center",fontsize=9)
            save_figure(fig,start+offset)
        rain_triptych([(label,d,phase_max,comp_max)],start+3)
    rain_triptych([("Observation: ERSSTv4 + GPCC",obs,.4,.7),("MIROC6 historical r10",model,.5,1.)],12)


def rain_triptych(rows,index):
    fig,axes=plt.subplots(len(rows),3,figsize=(15,6*len(rows)),subplot_kw={"projection":CRS},squeeze=False)
    for i,(label,d,pmax,cmax) in enumerate(rows):
        for j,(key,title,vmax) in enumerate([("pos_mean","PDO positive phase",pmax),("neg_mean","PDO negative phase",pmax),("composite","PDO+ minus PDO-",cmax)]):
            cf=field(axes[i,j],d,key,np.linspace(-vmax,vmax,21),dots=j==2)
            axes[i,j].set_title(f"{label}\n{title}",fontsize=10)
            cb=fig.colorbar(cf,ax=axes[i,j],orientation="horizontal",shrink=.9,pad=.08)
            cb.set_label("mm/day",fontsize=9)
    fig.suptitle("Task 2 | JJA rainfall anomalies and PDO phase composites\n1900-2014 | Land only | White dots: Welch p < 0.1",fontsize=14)
    fig.subplots_adjust(left=.05,right=.97,top=.85 if len(rows)==1 else .91,bottom=.06,wspace=.25,hspace=.35)
    save_figure(fig,index)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path)
    parser.add_argument("--compute",action="store_true")
    parser.add_argument("--compute-only",action="store_true")
    parser.add_argument("--section", choices=["all", "historical", "future"], default="all")
    args=parser.parse_args()
    if args.compute or args.compute_only:
        if args.data_root is None:
            parser.error("--data-root is required for computation")
        compute(args.data_root)
    if not args.compute_only:
        if args.section in ["all", "historical"]:
            hist=[load_case(f"historical_r{i}") for i in range(1,11)]
            obs=load_case("observation")
            plot_task2(obs,hist[9])
            plot_historical(hist,obs)
        if args.section in ["all", "future"]:
            plot_future()


if __name__ == "__main__":
    main()
