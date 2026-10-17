"""Summarize all six ACCESS WPSH runs and both below-ground extensions."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/task7_wpsh_eof"
FIGURES = ROOT / "figures/task7_wpsh_eof"


def save_figure(fig, name):
    for ext in ("png", "svg"):
        target = FIGURES / f"{name}.{ext}"
        fig.savefig(target, dpi=200)
        if ext == "svg":
            target.write_text("\n".join(line.rstrip() for line in target.read_text(encoding="utf-8").splitlines())+"\n", encoding="utf-8")


def main():
    rows = []
    fig, axes = plt.subplots(3, 2, figsize=(14, 9), sharex=True, sharey=True, layout="constrained")
    for column, scenario in enumerate(("ssp245", "ssp585")):
        for row, member in enumerate(("r1", "r2", "r5")):
            tag = f"access_{scenario}_{member}"
            primary = pd.read_csv(RESULTS / f"{tag}_index.csv")
            other = pd.read_csv(RESULTS / f"{tag}_nearest_index.csv")
            summary = json.loads((RESULTS / f"{tag}_summary.json").read_text())
            alternate = json.loads((RESULTS / f"{tag}_nearest_summary.json").read_text())
            if not np.array_equal(primary.year, other.year):
                raise ValueError("Sensitivity years differ")
            corr = float(np.corrcoef(primary.pc1, other.pc1)[0, 1])
            disagreement = primary.year[(primary.pc1_reference > 1.5) != (other.pc1_reference > 1.5)].tolist()
            centuries = summary["century_comparison"]
            item = dict(scenario=scenario, member=member,
                eof1_variance_percent=summary["eof1_variance_percent"],
                eof2_variance_percent=summary["eof2_variance_percent"],
                pc1_regional_sf_correlation=summary["pc1_regional_sf_correlation"],
                future_to_present_sd_ratio=centuries["sd_future"],
                sd_change_percent=centuries["sd_change_percent"],
                strong_events_present=centuries["strong_events_present"],
                strong_events_future=centuries["strong_events_future"],
                harmonic_nearest_pc1_correlation=corr,
                harmonic_nearest_pc1_reference_rmse=float(np.sqrt(np.mean((primary.pc1_reference-other.pc1_reference)**2))),
                fill_sensitive_event_years=disagreement,
                nearest_sd_change_percent=alternate["century_comparison"]["sd_change_percent"],
                nearest_events_present=alternate["century_comparison"]["strong_events_present"],
                nearest_events_future=alternate["century_comparison"]["strong_events_future"])
            rows.append(item)
            ax = axes[row, column]
            ax.axvspan(1900, 1999.5, color="#eeeeee", zorder=0)
            ax.plot(primary.year, primary.pc1_reference, color="#183c62", lw=1, label="Harmonic extension")
            ax.plot(other.year, other.pc1_reference, color="#d67728", lw=.7, alpha=.65, label="Nearest extension")
            ax.axhline(0, color="#777777", lw=.5)
            ax.axhline(1.5, color="#73468b", lw=.8, ls="--", label="Strong-event threshold")
            ax.set_title(f"{scenario.upper()} / {member} | EOF1 {summary['eof1_variance_percent']:.1f}% | fill correlation {corr:.6f}", fontsize=10, loc="left")
            ax.set_xlim(1900, 2099)
            ax.set_ylabel("PC1 (1900-1999 SD)", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(axis="y", alpha=.16)
    axes[0, 0].legend(loc="upper left", fontsize=7, ncol=1)
    fig.suptitle("ACCESS-CM2 | JJA 850-hPa WPSH from streamfunction EOF1\nMonthly quadratic detrending, 1900-2099; shaded years define the reference period", fontsize=14)
    save_figure(fig, "access_overview")
    plt.close(fig)
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "access_summary.csv", index=False)
    report = dict(reference_period="1900-1999", comparison_period="2000-2099",
                  strong_threshold="PC1 > 1.5 reference-period sample SD",
                  note="Six runs of one model, not six independent models; extensions do not recover real below-ground winds",
                  runs=rows)
    (RESULTS / "access_summary.json").write_text(json.dumps(report, indent=2)+"\n")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
    x = np.arange(len(rows))
    labels = [f"{r['scenario'].upper()}\n{r['member']}" for r in rows]
    axes[0].bar(x, table.sd_change_percent, color=["#3477a1"]*3+["#bc653b"]*3)
    axes[0].scatter(x, table.nearest_sd_change_percent, color="black", s=17, marker="x", label="Nearest-extension sensitivity")
    axes[0].axhline(0, color="#555555", lw=.7)
    axes[0].set_ylabel("Change in PC1 standard deviation (%)")
    axes[0].set_title("(a) 2000-2099 relative to 1900-1999", loc="left", fontsize=11)
    axes[0].legend(fontsize=8)
    axes[1].bar(x-.18, table.strong_events_present, width=.36, label="1900-1999", color="#3477a1")
    axes[1].bar(x+.18, table.strong_events_future, width=.36, label="2000-2099", color="#bc653b")
    axes[1].set_ylabel("Events per 100 years")
    axes[1].set_title("(b) PC1 > 1.5 reference-period SD", loc="left", fontsize=11)
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xticks(x, labels, fontsize=8)
        ax.grid(axis="y", alpha=.15)
        ax.set_axisbelow(True)
    fig.suptitle("ACCESS-CM2 WPSH | century comparisons conditional on wind extension", fontsize=13)
    save_figure(fig, "access_century_comparison")
    plt.close(fig)
    lines = ["# ACCESS-CM2 WPSH results", "",
             "Computed from global 850-hPa winds with harmonic below-ground extension.",
             "JJA EOF domain: 10–40°N, 110–180°E. Monthly quadratic detrending: 1900–2099.",
             "Reference: 1900–1999; future comparison: 2000–2099 (historical through 2014).", "",
             "| Scenario | Member | EOF1 variance | SD change | Strong events: reference → future | Harmonic/nearest PC1 correlation |",
             "|---|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['scenario'].upper()} | {r['member']} | {r['eof1_variance_percent']:.2f}% | {r['sd_change_percent']:+.2f}% | {r['strong_events_present']} → {r['strong_events_future']} | {r['harmonic_nearest_pc1_correlation']:.6f} |")
    lines += ["", "Strong events use PC1 > 1.5 reference-period sample SD. The separate", "full-period-threshold counts are retained in each run's summary JSON.", "",
              "These six runs are from one model with three members; they do not establish",
              "a multimodel consensus or statistical significance of a climate response.",
              "Below-ground winds are a numerical extension. Results remain conditional",
              "on this choice even when both extensions yield similar indices.", "",
              "## Threshold sensitivity", "",
              "Years whose strong-event classification changes between the two extensions:", ""]
    for r in rows:
        lines.append(f"- {r['scenario'].upper()} / {r['member']}: {', '.join(map(str, r['fill_sensitive_event_years'])) or 'none'}.")
    lines += ["", "[Overview](../figures/task7_wpsh_eof/access_overview.png) ·",
              "[Century comparison](../figures/task7_wpsh_eof/access_century_comparison.png) ·",
              "[Numerical summary](../results/task7_wpsh_eof/access_summary.csv) ·",
              "[Full method and reproduction](Task7_WPSH_streamfunction_EOF.md)", ""]
    (ROOT / "docs/Task7_ACCESS_results.md").write_text("\n".join(lines), encoding="utf-8")
    print(table.drop(columns="fill_sensitive_event_years").to_string(index=False))


if __name__ == "__main__":
    main()
