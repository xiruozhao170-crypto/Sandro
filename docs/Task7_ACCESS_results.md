# ACCESS-CM2 WPSH results

Computed from global 850-hPa winds with harmonic below-ground extension.
JJA EOF domain: 10–40°N, 110–180°E. Monthly quadratic detrending: 1900-2100.
Historical reference: 1900-2014 (115 years); future comparison: 2015-2100 (86 years).

Both climatological anomalies and index normalization use the historical period.
Each member/scenario retains one EOF fit over 1900-2100; ERA5/NCEP remain 1979–2020.

| Scenario | Member | EOF1 variance | SD change | Strong events: historical → future | Frequency (% of summers): historical → future | Harmonic/nearest PC1 correlation |
|---|---|---:|---:|---:|---:|---:|
| SSP245 | r1 | 76.02% | +19.98% | 5 → 12 | 4.35 → 13.95 | 0.999996 |
| SSP245 | r2 | 73.15% | +2.41% | 8 → 8 | 6.96 → 9.30 | 0.999996 |
| SSP245 | r5 | 72.50% | +4.39% | 9 → 10 | 7.83 → 11.63 | 0.999996 |
| SSP585 | r1 | 74.30% | +16.65% | 5 → 11 | 4.35 → 12.79 | 0.999995 |
| SSP585 | r2 | 73.55% | +15.83% | 9 → 10 | 7.83 → 11.63 | 0.999996 |
| SSP585 | r5 | 74.82% | +27.34% | 9 → 12 | 7.83 → 13.95 | 0.999996 |

Strong events use PC1 > 1.5 historical-period sample SD after historical-period centering.
Because the periods contain 115 and 86 summers, compare annual frequencies, not raw counts.
Frequency (%) = 100 × event count / number of summers. Counts are retained for transparency.
The separate full-period-threshold statistics are retained in each run's summary JSON.

These results supersede the earlier 1900–1999 versus 2000–2099 century comparison.
The future period now includes 2100, using additional winds from the same pinned source versions.
The two scenario fits can yield slightly different historical indices for the same member
because detrending and EOF fitting use each complete historical-plus-scenario record.

These six runs are from one model with three members; they do not establish
a multimodel consensus or statistical significance of a climate response.
Below-ground winds are a numerical extension. Results remain conditional
on this choice even when both extensions yield similar indices.

## Threshold sensitivity

Years whose strong-event classification changes between the two extensions:

- SSP245 / r1: none.
- SSP245 / r2: none.
- SSP245 / r5: none.
- SSP585 / r1: none.
- SSP585 / r2: none.
- SSP585 / r5: none.

[Overview](../figures/task7_wpsh_eof/access_overview.png) ·
[Historical/future comparison](../figures/task7_wpsh_eof/access_period_comparison.png) ·
[Numerical summary](../results/task7_wpsh_eof/access_summary.csv) ·
[Full method and reproduction](Task7_WPSH_streamfunction_EOF.md)
