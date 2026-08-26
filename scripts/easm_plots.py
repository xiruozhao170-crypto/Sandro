"""Shared figures for Tasks 2-4 (style follows the course example outputs)."""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def _sym_levels(*fields, n=21, q=99):
    lim = max(np.nanpercentile(np.abs(f), q) for f in fields)
    return np.linspace(-lim, lim, n)


def _map_axis(ax):
    # NB: gridlines(draw_labels=True) breaks constrained_layout with colorbars
    ax.coastlines(lw=1.0)
    ax.add_feature(cfeature.BORDERS, lw=0.4, ls=":")
    ax.set_xticks(np.arange(100, 131, 10), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(20, 45, 5), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.grid(lw=0.3, ls="--")


def rainfall_three_panels(comp, title, path):
    """(a) PDO+ anomaly, (b) PDO- anomaly, (c) composite diff w/ p<0.1 dots."""
    lat, lon = comp["pos"].lat, comp["pos"].lon
    lev_ab = _sym_levels(comp["pos"].values, comp["neg"].values)
    lev_c = _sym_levels(comp["diff"].values)

    fig, axes = plt.subplots(1, 3, figsize=(19, 6.2), constrained_layout=True,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    panels = [
        ("(a)  PDO Positive Phase\nJJA Rainfall Anomaly", comp["pos"], lev_ab),
        ("(b)  PDO Negative Phase\nJJA Rainfall Anomaly", comp["neg"], lev_ab),
        ("(c)  Composite (Pos $-$ Neg)\nwhite dots: p<0.1", comp["diff"], lev_c),
    ]
    for i, (ax, (ttl, fld, lev)) in enumerate(zip(axes, panels)):
        cf = ax.contourf(lon, lat, fld, levels=lev, cmap="BrBG",
                         transform=ccrs.PlateCarree(), extend="both")
        if i == 2:
            sig = comp["pval"].values < 0.1
            LO, LA = np.meshgrid(lon, lat)
            ax.plot(LO[sig], LA[sig], "o", color="white", ms=2.4, mew=0,
                    transform=ccrs.PlateCarree())
        _map_axis(ax)
        ax.set_title(ttl, fontsize=11)
        fig.colorbar(cf, ax=ax, label="mm/day", shrink=0.8, pad=0.03)
    fig.suptitle(title, fontsize=13)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def pdo_timeseries_panels(entries, path):
    """Stacked PDO-index panels: entries = [(label, pdo_dict), ...]."""
    n = len(entries)
    fig, axes = plt.subplots(n, 1, figsize=(13, 3.2 * n), sharex=True,
                             squeeze=False)
    for ax, (label, p) in zip(axes[:, 0], entries):
        raw, sm = p["pc_raw"], p["pc_smooth"]
        ax.fill_between(raw.index, raw.values, 0, where=raw.values > 0,
                        color="#e15759", alpha=0.75, lw=0)
        ax.fill_between(raw.index, raw.values, 0, where=raw.values <= 0,
                        color="#4e79a7", alpha=0.75, lw=0)
        ax.plot(sm.index, sm.values, "k", lw=2, label="9-yr running mean")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel("PDO index (s.d.)")
        ax.set_title(f"{label}   (pos={len(p['pos_years'])}yr  "
                     f"neg={len(p['neg_years'])}yr)")
        ax.legend(frameon=False, loc="upper right")
    axes[-1, 0].set_xlabel("Year")
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
