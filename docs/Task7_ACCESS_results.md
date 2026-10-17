# ACCESS-CM2 WPSH results

Computed from global 850-hPa winds with harmonic below-ground extension.
JJA EOF domain: 10–40°N, 110–180°E. Monthly quadratic detrending: 1900–2099.
Reference: 1900–1999; future comparison: 2000–2099 (historical through 2014).

| Scenario | Member | EOF1 variance | SD change | Strong events: reference → future | Harmonic/nearest PC1 correlation |
|---|---|---:|---:|---:|---:|
| SSP245 | r1 | 75.96% | +16.11% | 5 → 12 | 0.999996 |
| SSP245 | r2 | 73.43% | +6.51% | 6 → 11 | 0.999996 |
| SSP245 | r5 | 72.55% | +1.45% | 8 → 10 | 0.999996 |
| SSP585 | r1 | 74.40% | +13.87% | 4 → 10 | 0.999995 |
| SSP585 | r2 | 73.66% | +18.67% | 6 → 16 | 0.999996 |
| SSP585 | r5 | 74.84% | +20.00% | 8 → 12 | 0.999996 |

Strong events use PC1 > 1.5 reference-period sample SD. The separate
full-period-threshold counts are retained in each run's summary JSON.

These six runs are from one model with three members; they do not establish
a multimodel consensus or statistical significance of a climate response.
Below-ground winds are a numerical extension. Results remain conditional
on this choice even when both extensions yield similar indices.

## Threshold sensitivity

Years whose strong-event classification changes between the two extensions:

- SSP245 / r1: none.
- SSP245 / r2: 2085.
- SSP245 / r5: none.
- SSP585 / r1: 1973.
- SSP585 / r2: none.
- SSP585 / r5: none.

[Overview](../figures/task7_wpsh_eof/access_overview.png) ·
[Century comparison](../figures/task7_wpsh_eof/access_century_comparison.png) ·
[Numerical summary](../results/task7_wpsh_eof/access_summary.csv) ·
[Full method and reproduction](Task7_WPSH_streamfunction_EOF.md)
