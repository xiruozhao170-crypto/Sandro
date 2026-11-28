"""Check scientific invariants and the MIROC6 figure/land-mask artifacts."""
import hashlib
import json
import warnings

import numpy as np
from PIL import Image
import xarray as xr

from miroc6_report import RESULTS, ROOT, MEMBERS, SCENARIOS, load_case, field, RAIN_LEVELS, CRS
import matplotlib.pyplot as plt


def main():
    names = ["observation"] + [f"historical_r{i}" for i in range(1, 11)]
    names += [f"{s}_{m}" for s in SCENARIOS for m in MEMBERS]
    checks = []
    obs = load_case("observation")
    for name in names:
        d = load_case(name)
        future = name.startswith("ssp")
        start, end = (2015, 2100) if future else (1900, 2014)
        pc = d["pc1"]
        years = d["pc1_years"]
        expected_years = np.repeat(np.arange(start, end+1), 12) + np.tile((np.arange(1, 13)-.5)/12, end-start+1)
        np.testing.assert_allclose(years, expected_years, atol=1e-10)
        assert len(pc) == (end-start+1)*12
        assert np.isfinite(pc).all()
        # ERSST inputs/EOF normalization retain float32 rounding before export.
        assert abs(pc.mean()) < 5e-7 and abs(pc.std()-1) < 5e-7
        expected_smooth = xr.DataArray(pc, dims="time").rolling(time=108, center=True).mean().values
        np.testing.assert_allclose(d["pc1_smooth"], expected_smooth, atol=1e-10, equal_nan=True)
        year = np.floor(years).astype(int)
        month = np.floor((years-year)*12).astype(int)+1
        yy = np.unique(year)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            jja = np.array([np.nanmean(expected_smooth[(year==y)&np.isin(month,[6,7,8])]) for y in yy])
        np.testing.assert_array_equal(d["pos_years"], yy[jja>0])
        np.testing.assert_array_equal(d["neg_years"], yy[jja<0])
        np.testing.assert_array_equal(d["land_mask"], obs["land_mask"])
        np.testing.assert_array_equal(d["ec_lats"], obs["ec_lats"])
        np.testing.assert_array_equal(d["ec_lons"], obs["ec_lons"])
        assert d["land_mask"].any() and (~d["land_mask"]).any()
        for key in ["pos_mean", "neg_mean", "composite"]:
            assert np.isnan(d[key][~d["land_mask"]]).all(), (name,key)
            assert np.isfinite(d[key][d["land_mask"]]).all(), (name,key)
        np.testing.assert_allclose(d["composite"], d["pos_mean"]-d["neg_mean"], atol=1e-7, equal_nan=True)
        core = d["eof1"][np.ix_((d["eof1_lats"]>=35)&(d["eof1_lats"]<=45),
                                 (d["eof1_lons"]>=160)&(d["eof1_lons"]<=200))]
        assert np.nanmean(core) < 0
        assert 0 < float(d["var_frac"]) < 1
        checks.append(dict(case=name, months=len(pc), positive_years=len(d["pos_years"]),
                           negative_years=len(d["neg_years"]),
                           eof1_variance_percent=float(d["var_frac"])*100,
                           numerical_sha256=hashlib.sha256((RESULTS/f"{name}.npz").read_bytes()).hexdigest()))
    assert (len(obs["pos_years"]),len(obs["neg_years"])) == (53,54)
    ref = json.loads((RESULTS/"input_manifest.json").read_text())
    assert ref["future_members"] == MEMBERS
    assert max(ref["notebook_reference_case"]["max_absolute_errors"].values()) < 1e-7
    figures = json.loads((RESULTS/"figure_manifest.json").read_text())
    expected_indices = list(range(3,11))+list(range(12,21))+[22]
    assert set(map(int,figures)) == set(expected_indices)
    for index in expected_indices:
        entry=figures[str(index)]
        path=ROOT/entry["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        with Image.open(path) as im:
            assert list(im.size)==entry["image_size"]
            assert min(im.size)>1000
            im.verify()
        if index in [13,15]:
            assert entry["map_panels"]==12
            assert any("Observation" in title for title in entry["panel_titles"])
        if index==14:
            assert entry["pc_panels"]==11
            assert any("Observation" in title for title in entry["panel_titles"])
        if index in [17,18,19,20,22]:
            assert entry["map_panels"]+entry["pc_panels"]==12
            assert all(any(f"| {m}" in title for m in MEMBERS) for title in entry["panel_titles"])
            assert all(any(s.upper() in title for s in SCENARIOS) for title in entry["panel_titles"])
        if index in [3,4,5,6,7,8,9,10,12,15,16,19,20,22]:
            assert entry["ocean_overlay_on_every_map"]
    # Check actual rendered colors at three independent offshore locations,
    # including a plot with significance dots, rather than only mask flags.
    fig, ax = plt.subplots(figsize=(6,6), subplot_kw={"projection":CRS})
    field(ax, obs, "composite", RAIN_LEVELS, dots=True)
    fig.canvas.draw()
    pixels = np.asarray(fig.canvas.buffer_rgba())
    ocean_samples = []
    for lon, lat in [(124.3,19.7),(125.3,28.2),(123.3,33.2)]:
        x, y = ax.transData.transform((lon,lat))
        rgb = pixels[pixels.shape[0]-1-int(round(y)),int(round(x)),:3]
        assert np.all(rgb == 255), (lon,lat,rgb)
        ocean_samples.append(dict(lon=lon,lat=lat,rgb=rgb.tolist()))
    plt.close(fig)
    result=dict(cases=checks, updated_figures=len(figures),
                monthly_time_coverage_verified=True, phase_years_rederived=True,
                all_precipitation_ocean_grid_values_nan=True,
                every_precipitation_map_has_coastline_ocean_overlay=True,
                rendered_ocean_samples=ocean_samples,
                notebook_reference_case_reproduced=True,
                figure_hashes_and_required_panels_verified=True)
    (RESULTS/"validation.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
