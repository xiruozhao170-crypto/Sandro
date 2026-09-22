"""Task 6: WPSH position (5870-gpm contour of JJA 500-hPa height) by PDO phase.

Reproduces the design of Dong (2016, Atmos. Sci. Lett. 17:115-120,
doi:10.1002/asl.634) Figure 4 with this project's data: the JJA-mean
500-hPa geopotential height is composited over PDO-positive and
PDO-negative years and the western Pacific subtropical high (WPSH)
position is drawn as the 5870-gpm contour. JJA only.

Phase years are reused from the already-saved Task 2 / Task 4 results
(attrs pos_years / neg_years), so the composites here are phase-consistent
with the rainfall composites:
  observation  ERSSTv4 phases + NCEP/NCAR R1 hgt500 (1948-2014 only)
  historical   ACCESS-CM2 r1/r2/r5 zg500, 1900-2014
  ssp245/585   ACCESS-CM2 r1/r2/r5 zg500, 2015-2100

Figures:
  Task6_Figure1_WPSH_5870_JJA.png        2x2 obs / historical / ssp245 / ssp585
  Task6_Figure2_WPSH_5870_by_member.png  3x3 member detail with climatology
  Task6_Figure3_WPSH_bias_adjusted.png   2x2 as Fig 1 but model contours at
                                         5870+bias (bias = member 1948-2014
                                         climatology minus NCEP, box mean
                                         10-40N/110-180E) - ACCESS-CM2 zg500
                                         sits ~a few tens of gpm too low, so
                                         the raw 5870 line barely exists in
                                         the historical run
Outputs:
  results/task6_summary.json  (n years used, 5870 western edge per dataset)
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"
RES = REPO / "results"
FIG = REPO / "figures"

MEMBERS = ["r1", "r2", "r5"]
SCENARIOS = ["ssp245", "ssp585"]
THRESH = 5870.0
JJA = (6, 7, 8)

C_POS, C_NEG, C_CLIM = "#e15759", "#4e79a7", "0.45"
MEMBER_LS = {"r1": "-", "r2": "--", "r5": ":"}
EXTENT = [100, 200, 0, 55]


def phase_years(nc_path):
    ds = xr.open_dataset(nc_path)
    pos = np.array([int(y) for y in ds.attrs["pos_years"].split(",")])
    neg = np.array([int(y) for y in ds.attrs["neg_years"].split(",")])
    return pos, neg


def jja_by_year(paths, y0, y1):
    """Yearly JJA-mean zg500 [gpm] from one or more monthly files."""
    das = [xr.open_dataset(p)["zg500"] for p in paths]
    da = das[0] if len(das) == 1 else xr.concat(das, dim="time")
    da = da.sel(time=da.time.dt.month.isin(JJA))
    da = da.sel(time=(da.time.dt.year >= y0) & (da.time.dt.year <= y1))
    return da.groupby("time.year").mean("time").load()


def composites(yearly, pos, neg):
    pos = np.intersect1d(pos, yearly.year.values)
    neg = np.intersect1d(neg, yearly.year.values)
    return {
        "pos": yearly.sel(year=pos).mean("year"),
        "neg": yearly.sel(year=neg).mean("year"),
        "clim": yearly.mean("year"),
        "n_pos": len(pos), "n_neg": len(neg),
    }


def west_edge(field, thresh=THRESH):
    """Westernmost longitude reached by the >=5870 gpm area in 10-40N,
    100-180E (the classic WPSH west-extension index). None if absent."""
    sub = field.sel(lat=slice(10, 40), lon=slice(100, 180))
    hit = (sub >= thresh).any("lat")
    lons = sub.lon.values[hit.values]
    return round(float(lons.min()), 1) if lons.size else None


def wpsh_axis(ax):
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.coastlines(lw=0.8)
    ax.add_feature(cfeature.LAND, facecolor="0.93", zorder=0)
    ax.set_xticks(np.arange(100, 201, 20), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(0, 56, 10), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(labelsize=8)
    ax.grid(lw=0.3, ls="--", alpha=0.6)


def draw(ax, comp, color_pos=C_POS, color_neg=C_NEG, ls="-", lw=2.0,
         clim=False, thresh=THRESH):
    for key, color in [("pos", color_pos), ("neg", color_neg)]:
        f = comp[key]
        ax.contour(f.lon, f.lat, f, levels=[thresh], colors=[color],
                   linestyles=ls, linewidths=lw, transform=ccrs.PlateCarree())
    if clim:
        f = comp["clim"]
        ax.contour(f.lon, f.lat, f, levels=[thresh], colors=[C_CLIM],
                   linestyles=":", linewidths=1.4,
                   transform=ccrs.PlateCarree())


def missing_members(comps, exp, thresh=THRESH):
    """Members whose pos AND neg composites never reach `thresh` on the map."""
    out = []
    for m in MEMBERS:
        c = comps[(exp, m)]
        mx = max(float(c["pos"].max()), float(c["neg"].max()))
        if mx < thresh:
            out.append(f"{m} (max {mx:.0f})")
    return out


def main():
    # ---- gather data -----------------------------------------------------
    comps = {}  # key -> composites dict

    def box_mean(field):
        sub = field.sel(lat=slice(10, 40), lon=slice(110, 180))
        w = np.cos(np.deg2rad(sub.lat))
        return float(sub.weighted(w).mean())

    obs_pos, obs_neg = phase_years(RES / "task2_obs.nc")
    obs_yearly = jja_by_year([PROC / "zg500_ncep_194801-latest.nc"], 1948, 2014)
    comps["obs"] = composites(obs_yearly, obs_pos, obs_neg)
    obs_box = box_mean(obs_yearly.mean("year"))

    bias = {}  # member -> zg500 climatological bias vs NCEP [gpm], 1948-2014
    for m in MEMBERS:
        pos, neg = phase_years(RES / f"task2_{m}.nc")
        files = sorted(PROC.glob(f"zg500_historical_{m}i1p1f1_*.nc"))
        hist_yearly = jja_by_year(files, 1900, 2014)
        comps[("historical", m)] = composites(hist_yearly, pos, neg)
        bias[m] = round(box_mean(
            hist_yearly.sel(year=slice(1948, 2014)).mean("year")) - obs_box, 1)
        for scen in SCENARIOS:
            pos, neg = phase_years(RES / f"task4_{scen}_{m}.nc")
            files = sorted(PROC.glob(f"zg500_{scen}_{m}i1p1f1_*.nc"))
            comps[(scen, m)] = composites(
                jja_by_year(files, 2015, 2100), pos, neg)

    # ---- Figure 1: 2x2 overview ------------------------------------------
    proj = ccrs.PlateCarree(central_longitude=150)
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.6), constrained_layout=True,
                             subplot_kw={"projection": proj})
    panels = [("obs", "(a) Observation  NCEP/NCAR hgt500, ERSST phases "
               f"(1948–2014, pos={comps['obs']['n_pos']}yr "
               f"neg={comps['obs']['n_neg']}yr)"),
              ("historical", "(b) ACCESS-CM2 historical 1900–2014 "
               "(r1/r2/r5)"),
              ("ssp245", "(c) ACCESS-CM2 ssp245 2015–2100 (r1/r2/r5)"),
              ("ssp585", "(d) ACCESS-CM2 ssp585 2015–2100 (r1/r2/r5)")]
    for ax, (key, title) in zip(axes.flat, panels):
        wpsh_axis(ax)
        if key == "obs":
            draw(ax, comps["obs"], clim=True)
        else:
            for m in MEMBERS:
                draw(ax, comps[(key, m)], ls=MEMBER_LS[m], lw=1.8)
            miss = missing_members(comps, key)
            if miss:
                ax.text(0.02, 0.04, "no 5870-gpm area: " + ", ".join(miss),
                        transform=ax.transAxes, fontsize=8, color="0.35")
        ax.set_title(title, fontsize=10)
    handles = [Line2D([], [], color=C_POS, lw=2, label="PDO positive"),
               Line2D([], [], color=C_NEG, lw=2, label="PDO negative"),
               Line2D([], [], color=C_CLIM, lw=1.4, ls=":",
                      label="climatology (obs panel)")]
    handles += [Line2D([], [], color="k", lw=1.4, ls=MEMBER_LS[m],
                       label=f"member {m}") for m in MEMBERS]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=9,
               frameon=False)
    fig.suptitle("Task 6: WPSH position by PDO phase — 5870-gpm contour "
                 "of JJA-mean 500-hPa geopotential height\n"
                 "(after Dong 2016 ASL Fig. 4; phase years identical to the "
                 "Task 2 / Task 4 rainfall composites)", fontsize=12)
    fig.savefig(FIG / "Task6_Figure1_WPSH_5870_JJA.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 2: member detail 3x3 -------------------------------------
    rows = ["historical"] + SCENARIOS
    fig, axes = plt.subplots(3, 3, figsize=(15.5, 10.5),
                             constrained_layout=True,
                             subplot_kw={"projection": proj})
    for i, exp in enumerate(rows):
        for j, m in enumerate(MEMBERS):
            ax = axes[i, j]
            wpsh_axis(ax)
            c = comps[(exp, m)]
            draw(ax, c, clim=True)
            mx = max(float(c["pos"].max()), float(c["neg"].max()),
                     float(c["clim"].max()))
            if mx < THRESH:
                ax.text(0.03, 0.06, f"no 5870-gpm area (max {mx:.0f} gpm)",
                        transform=ax.transAxes, fontsize=9, color="0.35")
            ax.set_title(f"{exp} {m}  (pos={c['n_pos']}yr neg={c['n_neg']}yr)",
                         fontsize=10)
    handles = [Line2D([], [], color=C_POS, lw=2, label="PDO positive"),
               Line2D([], [], color=C_NEG, lw=2, label="PDO negative"),
               Line2D([], [], color=C_CLIM, lw=1.4, ls=":",
                      label="all-year JJA climatology")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=10,
               frameon=False)
    fig.suptitle("Task 6: WPSH 5870-gpm contour (JJA zg500) by member, "
                 "ACCESS-CM2\nhistorical 1900–2014 vs ssp245 / ssp585 "
                 "2015–2100", fontsize=13)
    fig.savefig(FIG / "Task6_Figure2_WPSH_5870_by_member.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 3: bias-adjusted thresholds for the model panels ----------
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.6), constrained_layout=True,
                             subplot_kw={"projection": proj})
    for ax, (key, title) in zip(axes.flat, panels):
        wpsh_axis(ax)
        if key == "obs":
            draw(ax, comps["obs"], clim=True)
            ax.set_title(title, fontsize=10)
        else:
            for m in MEMBERS:
                draw(ax, comps[(key, m)], ls=MEMBER_LS[m], lw=1.8,
                     thresh=THRESH + bias[m])
            ax.set_title(title.split(" (r1")[0]
                         + "  —  contour at 5870+$\\Delta_m$", fontsize=10)
    handles = [Line2D([], [], color=C_POS, lw=2, label="PDO positive"),
               Line2D([], [], color=C_NEG, lw=2, label="PDO negative"),
               Line2D([], [], color=C_CLIM, lw=1.4, ls=":",
                      label="climatology (obs panel)")]
    handles += [Line2D([], [], color="k", lw=1.4, ls=MEMBER_LS[m],
                       label=f"member {m}") for m in MEMBERS]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=9,
               frameon=False)
    fig.suptitle("Task 6 (bias-adjusted): model contours drawn at "
                 "5870+$\\Delta_m$ gpm, where $\\Delta_m$ = member JJA zg500 "
                 "climatology (1948–2014) minus NCEP, 10–40°N/110–180°E mean\n"
                 + ", ".join(f"$\\Delta_{{{m}}}$={bias[m]:+.1f} gpm"
                             for m in MEMBERS), fontsize=11)
    fig.savefig(FIG / "Task6_Figure3_WPSH_bias_adjusted.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # ---- summary ----------------------------------------------------------
    def entry(c, thresh=THRESH):
        return {"n_pos": c["n_pos"], "n_neg": c["n_neg"],
                "west_edge_pos": west_edge(c["pos"], thresh),
                "west_edge_neg": west_edge(c["neg"], thresh),
                "west_edge_clim": west_edge(c["clim"], thresh)}

    summary = {"threshold_gpm": THRESH, "season": "JJA",
               "west_edge_domain": "10-40N, 100-180E",
               "model_bias_gpm_vs_ncep_1948_2014": bias,
               "obs": entry(comps["obs"])}
    for exp in rows:
        summary[exp] = {m: entry(comps[(exp, m)]) for m in MEMBERS}
        summary[exp + "_bias_adjusted"] = {
            m: entry(comps[(exp, m)], THRESH + bias[m]) for m in MEMBERS}
    (RES / "task6_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
