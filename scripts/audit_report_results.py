"""Check report evidence against committed results; no raw climate rerun."""
from pathlib import Path
import json
import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / 'results'

def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m], b[m])[0, 1])

def main():
    report = {'scope': 'Audit of committed outputs, not rerun from raw data', 'historical': {}, 'future': {}}
    obs = xr.open_dataset(RES / 'task2_obs.nc')
    summary3 = json.loads((RES / 'task3_summary.json').read_text())
    summary4 = json.loads((RES / 'task4_summary.json').read_text())
    assert 'pattern_corr_vs_ensemble_mean_composite' not in summary3
    assert 'ssp585_vs_ssp245_ens_diff_corr' not in summary4
    for s in ['ssp245', 'ssp585']:
        assert 'ens_diff_corr_vs_hist' not in summary4['scenarios'][s]
    assert not (ROOT / 'figures/task1-4/Task4_T4-5_Composite_Difference_MemberMean.png').exists()
    rainfall = {phase: [] for phase in ['pos', 'neg', 'diff']}
    for m in ['r1','r2','r3','r4','r5']:
        h = xr.open_dataset(RES / f'task2_{m}.nc')
        np.testing.assert_allclose(h.rain_pos-h.rain_neg, h.rain_diff, atol=1e-12)
        for phase in rainfall:
            rainfall[phase].append(h[f'rain_{phase}'].values)
        rc = corr(h.rain_diff.values, obs.rain_diff.values)
        assert abs(rc-summary3['pattern_corr_vs_obs_composite'][m]) < 1e-10
        report['historical'][m] = {'rain_pattern_correlation_with_observation': rc,
            'n_pos':int(h.attrs['n_pos']), 'n_neg':int(h.attrs['n_neg'])}
        for s in ['ssp245','ssp585']:
            f = xr.open_dataset(RES / f'task4_{s}_{m}.nc')
            np.testing.assert_allclose(f.rain_pos-f.rain_neg, f.rain_diff, atol=1e-12)
            for key in ['n_pos','n_neg']:
                assert f.attrs[key] == summary4['scenarios'][s][m][key]
            report['future'][f'{s}_{m}'] = {
                'rain_pattern_correlation_with_own_historical': corr(f.rain_diff.values,h.rain_diff.values),
                'eof_correlation_with_observation_unweighted': corr(f.pdo_eof.values, obs.pdo_eof.values),
                'n_pos':int(f.attrs['n_pos']), 'n_neg':int(f.attrs['n_neg'])}
            assert round(report['future'][f'{s}_{m}']['rain_pattern_correlation_with_own_historical'], 3) == summary4['scenarios'][s][m]['rain_diff_corr_vs_own_hist']
            f.close()
        with xr.open_dataset(RES / f'task4_ssp245_{m}.nc') as low, xr.open_dataset(RES / f'task4_ssp585_{m}.nc') as high:
            assert round(corr(high.rain_diff.values, low.rain_diff.values), 3) == summary4['ssp585_vs_ssp245_diff_corr_by_member'][m]
        h.close()
    obs.close()
    with xr.open_dataset(RES / 'task3_spread.nc') as spread:
        for phase, fields in rainfall.items():
            assert f'rain_{phase}_mean' not in spread
            values = np.stack(fields)
            land = np.isfinite(values).all(axis=0)
            np.testing.assert_allclose(spread[f'rain_{phase}_var'].values[land], values[:, land].var(axis=0, ddof=1))
            assert np.isnan(spread[f'rain_{phase}_var'].values[~land]).all()
    report['rainfall_output_policy'] = 'Per-member composites and comparisons; member variance retained; rainfall member means and their diagnostics removed.'
    s5=json.loads((RES/'task5_summary.json').read_text())
    report['task5_matched_members_peak_band_fraction']={g: float(np.mean([
        s5[m if g=='historical' else f'{g} {m}']['cycle_energy_frac_of_variance']
        for m in ['r2','r3','r5']])) for g in ['historical','ssp245','ssp585']}
    report['task5_note']='Different peak bands and record lengths: descriptive fractions, not a fixed-band significance test.'
    # Independently reconstruct current equivalent thresholds from the earlier
    # committed annual height fields, whose raw inputs have a local manifest.
    summary6=json.loads((RES/'task6_summary.json').read_text())
    box=dict(lat=slice(10,40),lon=slice(110,180))
    with xr.open_dataset(RES/'task6_wpsh/obs.nc') as d:
        oc=d.z500_jja.mean('year').sortby('lat').sel(**box)
        weights=np.cos(np.deg2rad(oc.lat)).broadcast_like(oc)
        frac=float(weights.where(oc>=5870).sum()/weights.sum())
    report['task6_equivalent_threshold_reconstruction']={'obs_fraction':frac,'members':{}}
    for s in ['historical','ssp245','ssp585']:
        for m in ['r1','r2','r5']:
            with xr.open_dataset(RES/f'task6_wpsh/{s}_{m}.nc') as d:
                clim=d.z500_jja.mean('year').sel(**box)
                v=clim.values.ravel()
                w=np.cos(np.deg2rad(clim.lat)).broadcast_like(clim).values.ravel()
                order=np.argsort(v)[::-1]
                threshold=float(v[order][np.searchsorted(np.cumsum(w[order])/w.sum(),frac)])
                recorded=summary6[f'{s}_equivalent_threshold'][m]['threshold_gpm']
                assert abs(threshold-recorded)<=0.051, (s,m,threshold,recorded)
                edges={}
                for phase in ['pos','neg']:
                    field=d[f'z500_{phase}'].sel(lat=slice(10,40),lon=slice(100,180))
                    hit=(field>=threshold).any('lat')
                    edge=round(float(field.lon.values[hit.values].min()),1)
                    assert edge==summary6[f'{s}_equivalent_threshold'][m][f'west_edge_{phase}']
                    edges[phase]=edge
                report['task6_equivalent_threshold_reconstruction']['members'][f'{s}_{m}']={
                    'threshold_gpm_unrounded':threshold,'recorded_rounded':recorded,'west_edges':edges}
    dest=RES/'report_audit_2026-11-23.json'
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
