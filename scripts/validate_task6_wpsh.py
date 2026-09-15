"""Independent numerical checks for the submitted WPSH composites."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import netCDF4
import numpy as np
import pandas as pd
from PIL import Image
import xarray as xr

REPO = Path(__file__).resolve().parents[1]


def validate_figures():
    import contourpy
    from shapely.geometry import LineString, MultiLineString
    from task6_wpsh import draw_panel, EXTENT, plt, ccrs

    images = []
    expected = {"overview", "selected_members", "common_period", "ssp245_phase_separated",
                "ssp585_phase_separated", "observation_historical_phase_separated"}
    paths = sorted((REPO / "figures/task6_wpsh").glob("*.png"))
    assert {p.stem.removeprefix("Figure4_JJA_WPSH_") for p in paths} == expected
    for path in paths:
        with Image.open(path) as img:
            assert img.width >= 4000 and img.height >= 2500
            images.append(dict(file=path.name, width=img.width, height=img.height))
        svg = ET.parse(path.with_suffix(".svg"))
        assert not svg.findall(".//{http://www.w3.org/2000/svg}image"), "SVG must contain vector geometry"
    audit = pd.read_csv(REPO / "results/task6_wpsh/contour_display_audit.csv")
    assert len(audit) == 32 and (audit.missing_values == 0).all()
    assert (audit.new_visible_length_degrees + 1e-8 >= audit.old_visible_length_degrees).all()
    checked = 0
    for path in sorted((REPO / "results/task6_wpsh").glob("*.nc")):
        with xr.open_dataset(path) as ds:
            for phase in ("pos", "neg"):
                f = ds[f"z500_{phase}"]
                expected_segments = contourpy.contour_generator(x=f.lon.values, y=f.lat.values,
                                                                z=f.values).lines(5870)
                fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree(central_longitude=140)})
                draw_panel(ax, ds, "Rendering check", phases=(phase,))
                np.testing.assert_allclose(ax.get_extent(ccrs.PlateCarree(central_longitude=140)),
                                           [EXTENT[0]-140, EXTENT[1]-140, EXTENT[2], EXTENT[3]], atol=1e-7)
                actual = [c for c in ax.collections if hasattr(c, "levels")]
                if expected_segments:
                    assert len(actual) == 1
                    np.testing.assert_array_equal(actual[0].levels, [5870])
                    expected_lines = MultiLineString([LineString(s) for s in expected_segments])
                    drawn_lines = MultiLineString([LineString(s) for s in actual[0].allsegs[0] if len(s)>1])
                    assert expected_lines.hausdorff_distance(drawn_lines) < 1e-8
                else:
                    assert not actual and any("no 5870-gpm contour" in t.get_text() for t in ax.texts)
                plt.close(fig)
                checked += 1
    return images, dict(phase_panels_checked=checked, contour_geometry="Matches native-grid 5870-gpm contours",
                        map_extent_degrees=list(EXTENT), missing_composite_values=0,
                        clipping="Expanded map shows at least as much of every curve as the old map",
                        absent_contours="Explicitly labelled; no replacement threshold or displaced lines")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--figures-only", action="store_true")
    args = parser.parse_args()
    root = REPO / "results/task6_wpsh"
    if args.figures_only:
        result = json.loads((root / "validation.json").read_text())
        assert result["status"] == "passed"
        result["figures"], result["rendering_recheck"] = validate_figures()
        result["svg"] = "All six SVG figures contain no embedded raster images"
        result["numerical_checks"] = "Previous raw-season and composite checks retained; this rerun checks rendering only"
        (root / "validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(json.dumps(result["rendering_recheck"], indent=2))
        return
    if args.data_root is None:
        parser.error("--data-root is required unless --figures-only is used")
    manifest = json.loads((root / "input_manifest.json").read_text())
    assert manifest["selected_members"] == ["r1", "r2", "r5"]
    assert len(manifest["model_files"]) == 14
    assert sum(r["used_for_composite"] for r in manifest["model_files"]) == 12
    table = pd.read_csv(root / "phase_years.csv")
    summary = json.loads((root / "summary.json").read_text())
    geojson = json.loads((root / "contours_5870gpm.geojson").read_text())
    checks = []
    for path in sorted(root.glob("*.nc")):
        with xr.open_dataset(path) as ds:
            assert ds.attrs["pressure_hPa"] == 500
            assert ds.attrs["wpsh_contour_gpm"] == 5870
            for phase in ("pos", "neg"):
                field = ds[f"z500_{phase}"]
                exists = bool(float(field.min()) < 5870 < float(field.max()))
                assert summary[path.stem]["contour_exists_in_computed_domain"][phase] == exists
                emitted = any(f["properties"]["dataset"] == path.stem and f["properties"]["phase"] == phase
                              for f in geojson["features"])
                assert emitted == exists
            np.testing.assert_allclose(ds.z500_diff, ds.z500_pos - ds.z500_neg, atol=1e-10)
            if "year" not in ds.dims:
                continue
            for phase, code in [("pos", 1), ("neg", -1)]:
                mask = ds.phase.values == code
                years = ds.year.values[mask]
                assert len(years) == ds.attrs[f"n_{phase}"]
                expected = ds.z500_jja.values[mask].sum(axis=0) / len(years)
                np.testing.assert_allclose(ds[f"z500_{phase}"], expected, atol=1e-10)
                rows = table[(table.dataset == path.stem) & (table.pdo_phase == code)]
                np.testing.assert_array_equal(rows.year.values, years)
            if path.stem.startswith(("historical_r", "ssp245_r", "ssp585_r")):
                experiment, member = path.stem.split("_")
                # Verify a full source-grid season from EACH phase directly against raw files.
                for code in [1, -1]:
                    year = int(ds.year.values[np.flatnonzero(ds.phase.values == code)[0]])
                    raw_path = next(args.data_root / r["file"] for r in manifest["model_files"]
                                    if r["experiment"] == experiment and r["member"].startswith(member + "i")
                                    and int(r["first_month"][:4]) <= year <= int(r["last_month"][:4]))
                    with netCDF4.Dataset(raw_path) as raw:
                        t = raw["time"]
                        dates = netCDF4.num2date(t[:], t.units, t.calendar)
                        it = [i for i, d in enumerate(dates) if d.year == year and d.month in (6, 7, 8)]
                        assert len(it) == 3
                        iz = np.flatnonzero(np.isclose(raw["plev"][:], 50000, rtol=0, atol=.001)).item()
                        la = np.searchsorted(raw["lat"][:], ds.lat.values)
                        lo = np.searchsorted(raw["lon"][:], ds.lon.values)
                        direct = np.stack([raw["zg"][i, iz, la, lo].astype("float64") for i in it]).mean(axis=0)
                    error = float(np.max(abs(direct - ds.z500_jja.sel(year=year).values)))
                    assert error < 1e-9
                    checks.append(dict(dataset=path.stem, year=year, phase=code, max_raw_JJA_error_gpm=error))
    for experiment in ("historical", "ssp245", "ssp585"):
        with xr.open_dataset(root / f"{experiment}_member_mean.nc") as ensemble:
            for phase in ("pos", "neg"):
                fields = []
                for member in ("r1", "r2", "r5"):
                    with xr.open_dataset(root / f"{experiment}_{member}.nc") as member_ds:
                        fields.append(member_ds[f"z500_{phase}"].values)
                np.testing.assert_allclose(ensemble[f"z500_{phase}"], np.stack(fields).mean(axis=0), atol=1e-10)
    with xr.open_dataset(REPO / "data/observations/ncep_hgt500_monthly_1948_2014.nc") as raw:
        with xr.open_dataset(root / "obs.nc") as saved:
            direct = raw.hgt.sel(level=500, time=slice("1948-06-01", "1948-08-31")).sortby("lat").astype("float64").mean("time")
            np.testing.assert_allclose(direct, saved.z500_jja.sel(year=1948), atol=1e-10)
    images, rendering = validate_figures()
    assert len(images) == 6 and len(checks) == 18
    result = dict(status="passed", raw_model_seasons_checked=checks,
                  composite_reconstruction="Passed for all 13 annual datasets and both PDO phases",
                  member_means="Passed for historical, SSP245, SSP585",
                  observation_raw_season="1948 JJA agrees with downloaded NCEP data",
                  csv_year_lists="Agree with NetCDF phase lists",
                  contour_availability="GeoJSON and summary agree with each composite's range; below-threshold fields have no invented contours",
                  inputs="14 inventoried, 12 used, correct members/500 hPa/5870 gpm",
                  figures=images, rendering_recheck=rendering,
                  svg="All six SVG figures contain no embedded raster images")
    (root / "validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
