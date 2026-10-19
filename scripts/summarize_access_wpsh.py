"""Summarize all six ACCESS WPSH runs and both below-ground extensions."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from wpsh_streamfunction import FUTURE_PERIOD, HISTORICAL_PERIOD

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/task7_wpsh_eof"
FIGURES = ROOT / "figures/task7_wpsh_eof"
HISTORICAL_LABEL = "-".join(map(str, HISTORICAL_PERIOD))
FUTURE_LABEL = "-".join(map(str, FUTURE_PERIOD))


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
            periods = summary["period_comparison"]
            nearest = alternate["period_comparison"]
            for run in (summary, alternate):
                if run["baseline"] != list(HISTORICAL_PERIOD):
                    raise ValueError(f"Stale historical baseline: {run['tag']}")
                if (run["period_comparison"]["historical"], run["period_comparison"]["future"]) != (HISTORICAL_LABEL, FUTURE_LABEL):
                    raise ValueError(f"Stale comparison periods: {run['tag']}")
            item = dict(scenario=scenario, member=member,
                historical_period=periods["historical"], future_period=periods["future"],
                n_years_historical=periods["n_years_historical"], n_years_future=periods["n_years_future"],
                eof1_variance_percent=summary["eof1_variance_percent"],
                eof2_variance_percent=summary["eof2_variance_percent"],
                pc1_regional_sf_correlation=summary["pc1_regional_sf_correlation"],
                future_to_historical_sd_ratio=periods["future_to_historical_sd_ratio"],
                sd_change_percent=periods["sd_change_percent"],
                strong_events_historical=periods["strong_events_historical"],
                strong_events_future=periods["strong_events_future"],
                strong_event_rate_historical_percent=periods["strong_event_rate_historical_percent"],
                strong_event_rate_future_percent=periods["strong_event_rate_future_percent"],
                harmonic_nearest_pc1_correlation=corr,
                harmonic_nearest_pc1_reference_rmse=float(np.sqrt(np.mean((primary.pc1_reference-other.pc1_reference)**2))),
                fill_sensitive_event_years=disagreement,
                nearest_sd_change_percent=nearest["sd_change_percent"],
                nearest_events_historical=nearest["strong_events_historical"],
                nearest_events_future=nearest["strong_events_future"],
                nearest_event_rate_historical_percent=nearest["strong_event_rate_historical_percent"],
                nearest_event_rate_future_percent=nearest["strong_event_rate_future_percent"])
            rows.append(item)
            ax = axes[row, column]
            ax.axvspan(HISTORICAL_PERIOD[0], HISTORICAL_PERIOD[1]+.5, color="#eeeeee", zorder=0)
            ax.plot(primary.year, primary.pc1_reference, color="#183c62", lw=1, label="Harmonic extension")
            ax.plot(other.year, other.pc1_reference, color="#d67728", lw=.7, alpha=.65, label="Nearest extension")
            ax.axhline(0, color="#777777", lw=.5)
            ax.axhline(1.5, color="#73468b", lw=.8, ls="--", label="Strong-event threshold")
            ax.set_title(f"{scenario.upper()} / {member} | EOF1 {summary['eof1_variance_percent']:.1f}% | fill correlation {corr:.6f}", fontsize=10, loc="left")
            ax.set_xlim(1900, 2099)
            ax.set_ylabel(f"PC1 ({HISTORICAL_LABEL} SD)", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(axis="y", alpha=.16)
    axes[0, 0].legend(loc="upper left", fontsize=7, ncol=1)
    fig.suptitle(f"ACCESS-CM2 | JJA 850-hPa WPSH from streamfunction EOF1\nMonthly quadratic detrending, 1900-2099; shaded historical reference: {HISTORICAL_LABEL}", fontsize=14)
    save_figure(fig, "access_overview")
    plt.close(fig)
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "access_summary.csv", index=False)
    report = dict(reference_period=HISTORICAL_LABEL, comparison_period=FUTURE_LABEL,
                  historical_years=HISTORICAL_PERIOD[1]-HISTORICAL_PERIOD[0]+1,
                  future_years=FUTURE_PERIOD[1]-FUTURE_PERIOD[0]+1,
                  strong_threshold="PC1 > 1.5 historical-period sample SD",
                  event_rate_units="Percent of summers; 100 * strong-event count / years in each period",
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
    axes[0].set_title(f"(a) {FUTURE_LABEL} relative to {HISTORICAL_LABEL}", loc="left", fontsize=11)
    axes[0].legend(fontsize=8)
    axes[1].bar(x-.18, table.strong_event_rate_historical_percent, width=.36, label=f"{HISTORICAL_LABEL} (115 years)", color="#3477a1")
    axes[1].bar(x+.18, table.strong_event_rate_future_percent, width=.36, label=f"{FUTURE_LABEL} (85 years)", color="#bc653b")
    axes[1].set_ylabel("Strong-event frequency (% of summers)")
    axes[1].set_title("(b) PC1 > 1.5 historical-period SD", loc="left", fontsize=11)
    axes[1].margins(y=.25)
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xticks(x, labels, fontsize=8)
        ax.grid(axis="y", alpha=.15)
        ax.set_axisbelow(True)
    fig.suptitle("ACCESS-CM2 WPSH | historical and future comparison", fontsize=13)
    save_figure(fig, "access_period_comparison")
    plt.close(fig)
    lines = ["# ACCESS-CM2 WPSH results", "",
             "Computed from global 850-hPa winds with harmonic below-ground extension.",
             "JJA EOF domain: 10–40°N, 110–180°E. Monthly quadratic detrending: 1900–2099.",
             "Historical reference: 1900–2014 (115 years); future comparison: 2015–2099 (85 years).", "",
             "Both climatological anomalies and index normalization use the historical period.",
             "Each member/scenario retains one EOF fit over 1900–2099; ERA5/NCEP remain 1979–2020.", "",
             "| Scenario | Member | EOF1 variance | SD change | Strong events: historical → future | Frequency (% of summers): historical → future | Harmonic/nearest PC1 correlation |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['scenario'].upper()} | {r['member']} | {r['eof1_variance_percent']:.2f}% | {r['sd_change_percent']:+.2f}% | {r['strong_events_historical']} → {r['strong_events_future']} | {r['strong_event_rate_historical_percent']:.2f} → {r['strong_event_rate_future_percent']:.2f} | {r['harmonic_nearest_pc1_correlation']:.6f} |")
    lines += ["", "Strong events use PC1 > 1.5 historical-period sample SD after historical-period centering.",
              "Because the periods contain 115 and 85 summers, compare annual frequencies, not raw counts.",
              "Frequency (%) = 100 × event count / number of summers. Counts are retained for transparency.",
              "The separate full-period-threshold statistics are retained in each run's summary JSON.", "",
              "These results supersede the earlier 1900–1999 versus 2000–2099 century comparison.",
              "The two scenario fits can yield slightly different historical indices for the same member",
              "because detrending and EOF fitting use each complete historical-plus-scenario record.", "",
              "These six runs are from one model with three members; they do not establish",
              "a multimodel consensus or statistical significance of a climate response.",
              "Below-ground winds are a numerical extension. Results remain conditional",
              "on this choice even when both extensions yield similar indices.", "",
              "## Threshold sensitivity", "",
              "Years whose strong-event classification changes between the two extensions:", ""]
    for r in rows:
        lines.append(f"- {r['scenario'].upper()} / {r['member']}: {', '.join(map(str, r['fill_sensitive_event_years'])) or 'none'}.")
    lines += ["", "[Overview](../figures/task7_wpsh_eof/access_overview.png) ·",
              "[Historical/future comparison](../figures/task7_wpsh_eof/access_period_comparison.png) ·",
              "[Numerical summary](../results/task7_wpsh_eof/access_summary.csv) ·",
              "[Full method and reproduction](Task7_WPSH_streamfunction_EOF.md)", ""]
    (ROOT / "docs/Task7_ACCESS_results.md").write_text("\n".join(lines), encoding="utf-8")
    print(table.drop(columns="fill_sensitive_event_years").to_string(index=False))


if __name__ == "__main__":
    main()
