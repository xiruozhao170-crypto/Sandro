"""Task 6: Fourier power-spectrum analysis of the PDO PC1 time series.

Follows Figure 1c of the advisor's paper (Atmosphere 2020, 11, 3,
doi:10.3390/atmos11010003): the periodogram (Fourier transform) of the
monthly normalized PC1 is plotted against period, together with the
theoretical red-noise (AR1) spectrum fitted from the lag-1
autocorrelation (Gilman et al. 1963) and its 90%/95% chi-square
confidence levels.

The period of the dominant cycle is the period of the strongest spectral
peak that stands above the confidence levels.  The "energy" of that
cycle is the variance it carries: the integral of the power spectral
density over the contiguous band around the peak where the spectrum
exceeds the red-noise background (by Parseval the full integral equals
the total variance of PC1, = 1 for the normalized index).

Series analysed
  historical (1900-2014): observed ERSSTv4 PC1 + ACCESS-CM2 r1-r5
  future (2015-2100):     ssp245 and ssp585, the best-3 members of Task 3

Figures:
  Task6_Figure1_Historical_PC1_Spectrum.png   2x3: obs + r1-r5
  Task6_Figure2_Future_PC1_Spectrum.png       2x3: ssp245/ssp585 x best3
  Task6_Figure3_Cycle_Energy_Summary.png      dominant period + energy bars
Outputs: results/task6_summary.json
"""
import json
from pathlib import Path

import numpy as np
import xarray as xr
from scipy import signal, stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"
FIG = REPO / "figures"

MEMBERS = ["r1", "r2", "r3", "r4", "r5"]
SCENARIOS = ["ssp245", "ssp585"]
DT_YEARS = 1.0 / 12.0  # monthly data


def analyze_pc(pc):
    """Periodogram + AR1 red-noise background and cycle diagnostics.

    Returns a dict with frequencies (cycles/yr), PSD (sigma^2 * yr), the
    red-noise spectrum, 90/95% confidence spectra, and the dominant-cycle
    period/energy.
    """
    x = np.asarray(pc, dtype=float)
    x = x - x.mean()
    n = x.size

    freq, psd = signal.periodogram(x, fs=1.0 / DT_YEARS)
    freq, psd = freq[1:], psd[1:]  # drop the zero-frequency (mean) bin
    # modest Daniell (5-point running mean) smoothing: the raw periodogram
    # has only 2 dof per bin and single-bin noise spikes cross any
    # confidence level; smoothing gives ~10 dof and stable peaks
    NSMOOTH = 5
    kernel = np.ones(NSMOOTH) / NSMOOTH
    psd = np.convolve(np.pad(psd, NSMOOTH // 2, mode="reflect"),
                      kernel, mode="valid")
    df = freq[1] - freq[0]
    total_var = psd.sum() * df

    # Gilman et al. (1963) theoretical red-noise spectrum from the lag-1
    # autocorrelation, scaled to carry the same total variance
    r1 = float(np.corrcoef(x[:-1], x[1:])[0, 1])
    red = (1.0 - r1**2) / (1.0 + r1**2
                           - 2.0 * r1 * np.cos(np.pi * freq / freq[-1]))
    red *= psd.mean() / red.mean()

    # 2 dof per raw estimate x 5-point Daniell smoothing
    dof = 2 * NSMOOTH
    conf90 = red * stats.chi2.ppf(0.90, dof) / dof
    conf95 = red * stats.chi2.ppf(0.95, dof) / dof

    # dominant cycle: strongest peak above the 95% level (fall back to
    # 90%, then to the strongest peak above the red-noise background).
    # Search window: periods between 2 yr and 1/3 of the record length --
    # longer periods complete <3 cycles in the record and cannot be
    # separated from the trend, and the 1-yr line is residual seasonality.
    period = 1.0 / freq
    window = (period >= 2.0) & (period <= n * DT_YEARS / 3.0)
    for level, name in [(conf95, "95%"), (conf90, "90%"), (red, "red noise")]:
        above = (psd > level) & window
        if above.any():
            ipk = int(np.argmax(np.where(above, psd, -np.inf)))
            sig_level = name
            break
    else:
        ipk = int(np.argmax(np.where(window, psd, -np.inf)))
        sig_level = "not significant"

    # energy of the cycle: integrate the PSD over the contiguous band
    # around the peak where it exceeds the red-noise background, keeping
    # the band inside the resolvable range (period <= half the record)
    resolvable = period <= n * DT_YEARS / 2.0
    i0 = ipk
    while i0 > 0 and psd[i0 - 1] > red[i0 - 1] and resolvable[i0 - 1]:
        i0 -= 1
    i1 = ipk
    while i1 < psd.size - 1 and psd[i1 + 1] > red[i1 + 1]:
        i1 += 1
    band = slice(i0, i1 + 1)
    energy = float(psd[band].sum() * df)

    return {
        "freq": freq, "psd": psd, "red": red,
        "conf90": conf90, "conf95": conf95, "r1": r1,
        "peak_period": float(1.0 / freq[ipk]),
        "peak_psd": float(psd[ipk]),
        "sig_level": sig_level,
        "band_period": (float(1.0 / freq[i1]), float(1.0 / freq[i0])),
        "energy": energy,
        "energy_frac": float(energy / total_var),
        "total_var": float(total_var),
    }


def spectrum_axis(ax, sp, title):
    per = 1.0 / sp["freq"]
    ax.plot(per, sp["psd"], "k", lw=1.0, label="PC1 spectrum")
    ax.plot(per, sp["red"], "r", lw=1.6, label="Red noise")
    ax.plot(per, sp["conf90"], color="0.55", ls="--", lw=1.4,
            label="90% confidence")
    ax.plot(per, sp["conf95"], color="0.55", ls="-", lw=1.4,
            label="95% confidence")
    ax.axvline(sp["peak_period"], color="#e15759", ls=":", lw=1.2)
    ax.set_xscale("log", base=2)
    ax.set_xlim(64, 0.5)
    ax.set_xticks([64, 32, 16, 8, 4, 2, 1])
    ax.set_xticklabels(["64", "32", "16", "8", "4", "2", "1"])
    ax.set_title(title, fontsize=11)
    ax.grid(lw=0.3, ls="--", alpha=0.6)
    txt = (f"cycle ≈ {sp['peak_period']:.1f} yr ({sp['sig_level']})\n"
           f"energy = {sp['energy']:.2f} σ² "
           f"({sp['energy_frac'] * 100:.0f}% of var)")
    ax.text(0.03, 0.95, txt, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(fc="white", ec="0.7", alpha=0.85))


def spectrum_figure(entries, path, suptitle):
    """entries = [(label, sp_dict), ...] laid out 2x3."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), sharex=True)
    for ax, (label, sp) in zip(axes.flat, entries):
        spectrum_axis(ax, sp, label)
    for ax in axes.flat[len(entries):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("PSD (σ² · yr)")
    for ax in axes[-1, :]:
        ax.set_xlabel("Period (years)")
    axes.flat[0].legend(frameon=False, fontsize=8, loc="center left")
    fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def summary_figure(order, specs, groups, path):
    colors = {"obs": "#4d4d4d", "historical": "#4e79a7",
              "ssp245": "#f28e2b", "ssp585": "#e15759"}
    y = np.arange(len(order))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=True)
    for ax, key, xlabel in [
            (axes[0], "peak_period", "Dominant cycle period (years)"),
            (axes[1], "energy_frac", "Cycle energy (fraction of variance)")]:
        vals = [specs[k][key] for k in order]
        ax.barh(y, vals, height=0.62,
                color=[colors[groups[k]] for k in order])
        for yi, v in zip(y, vals):
            ax.text(v, yi, f" {v:.2f}", va="center", fontsize=9)
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", lw=0.3, ls="--", alpha=0.6)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(order)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors.values()]
    axes[1].legend(handles, colors.keys(), frameon=False, fontsize=9,
                   loc="lower right")
    fig.suptitle("Task 6: PDO PC1 dominant cycle and its energy\n"
                 "(periodogram peak above red-noise background)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    best3 = json.loads((RES / "task3_summary.json").read_text())["best3"]

    specs, groups = {}, {}

    # historical: observations + all five members
    for tag, label in [("obs", "obs (ERSSTv4)")] + [(m, m) for m in MEMBERS]:
        ds = xr.open_dataset(RES / f"task2_{tag}.nc")
        specs[label] = analyze_pc(ds["pdo_pc"].values)
        groups[label] = "obs" if tag == "obs" else "historical"

    # future: the best-3 members per scenario (Task 4 outputs)
    for scen in SCENARIOS:
        for m in best3:
            ds = xr.open_dataset(RES / f"task4_{scen}_{m}.nc")
            label = f"{scen} {m}"
            specs[label] = analyze_pc(ds["pdo_pc"].values)
            groups[label] = scen

    hist_labels = ["obs (ERSSTv4)"] + MEMBERS
    fut_labels = [f"{scen} {m}" for scen in SCENARIOS for m in best3]

    spectrum_figure(
        [(f"{l}   (r1={specs[l]['r1']:.2f})", specs[l]) for l in hist_labels],
        FIG / "Task6_Figure1_Historical_PC1_Spectrum.png",
        "Task 6 Figure 1: Power spectrum of the historical PDO PC1 "
        "(1900–2014, monthly)\nblack: periodogram | red: AR1 red noise | "
        "gray: 90% (dashed) and 95% (solid) confidence levels")

    spectrum_figure(
        [(f"{l}   (r1={specs[l]['r1']:.2f})", specs[l]) for l in fut_labels],
        FIG / "Task6_Figure2_Future_PC1_Spectrum.png",
        "Task 6 Figure 2: Power spectrum of the future PDO PC1 "
        "(2015–2100, monthly, best-3 members)\nblack: periodogram | red: AR1 "
        "red noise | gray: 90% (dashed) and 95% (solid) confidence levels")

    order = hist_labels + fut_labels
    summary_figure(order, specs, groups,
                   FIG / "Task6_Figure3_Cycle_Energy_Summary.png")

    summary = {}
    for l in order:
        s = specs[l]
        summary[l] = {
            "group": groups[l],
            "lag1_autocorr": round(s["r1"], 3),
            "dominant_period_years": round(s["peak_period"], 2),
            "significance": s["sig_level"],
            "band_period_years": [round(s["band_period"][0], 2),
                                  round(s["band_period"][1], 2)],
            "cycle_energy_sigma2": round(s["energy"], 4),
            "cycle_energy_frac_of_variance": round(s["energy_frac"], 4),
        }
        print(f"{l:18s} cycle={s['peak_period']:6.1f} yr "
              f"({s['sig_level']:>9s})  energy={s['energy']:.3f} σ² "
              f"({s['energy_frac'] * 100:4.1f}% of var)")

    for grp in ["historical", "ssp245", "ssp585"]:
        keys = [l for l in order if groups[l] == grp]
        summary[f"_mean_{grp}"] = {
            "dominant_period_years": round(float(np.mean(
                [specs[l]["peak_period"] for l in keys])), 2),
            "cycle_energy_frac_of_variance": round(float(np.mean(
                [specs[l]["energy_frac"] for l in keys])), 4),
        }

    (RES / "task6_summary.json").write_text(json.dumps(summary, indent=2))
    print("saved results/task6_summary.json")


if __name__ == "__main__":
    main()
