"""Read local manifests/NetCDF headers and sample paired ACCESS zg downloads.

This does not rehash the 18.9 GB raw archive or compare every time step.
Usage: python scripts/audit_report_data_provenance.py --root "E:/Sandro Dachuang"
"""
from pathlib import Path
import argparse
import hashlib
import json
from datetime import datetime, timezone
import netCDF4
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    args = ap.parse_args()
    root = args.root
    repo = Path(__file__).resolve().parents[1]
    subset_root = root / 'ACCESS-CM2_zg_script_20260911'
    first = json.loads((subset_root / 'completion_manifest.json').read_text(encoding='utf-8-sig'))
    second = json.loads((root / 'download_job_20260912/original_download_status.json').read_text(encoding='utf-8-sig'))
    later = {f['name']: f for f in second['files']}
    inputs = json.loads((repo / 'results/task6_wpsh/input_manifest.json').read_text(encoding='utf-8'))
    input_map = {f['file']: f for f in inputs['model_files']}
    out = {'audited_utc': datetime.now(timezone.utc).isoformat(),
           'scope': 'Manifest checksum agreement, current file headers and three whole-grid time slices per paired file; not a new full-file checksum or all-month numerical comparison.',
           'subset_completed_utc': first['completed_utc'],
           'original_completed_utc': second['updated_utc'], 'paired_zg': [], 'tos_pr_headers': [], 'early_task6_figure_matches': []}
    for entry in first['results']:
        rec = entry['original_script_record']
        name = rec['name']
        raw_path = root / name
        sub_path = subset_root / entry['file']
        item = {'file': name, 'publication_version': rec['version'], 'source_url': rec['url'],
                'original_sha256_recorded': rec['sha256'],
                'original_sha256_agrees_between_manifests': rec['sha256'] == later[name]['expected_sha256'],
                'subset_relative_path': str(sub_path.relative_to(root)), 'subset_sha256_recorded': entry['subset_sha256'],
                'raw_bytes': raw_path.stat().st_size,
                'raw_bytes_match_download_record': raw_path.stat().st_size == later[name]['total_bytes']}
        with netCDF4.Dataset(raw_path) as raw, netCDF4.Dataset(sub_path) as sub:
            k = int(np.argmin(abs(raw['plev'][:] - 50000)))
            n = len(raw.dimensions['time'])
            times = sorted(set([0, n // 2, n - 1]))
            a = np.ma.filled(raw['zg'][times, k, :, :], np.nan)
            b = np.ma.filled(sub['zg'][times, 0, :, :], np.nan)
            item.update(tracking_id=getattr(raw, 'tracking_id', None),
                        source_header_version=getattr(raw, 'version', None),
                        pressure_Pa=float(raw['plev'][k]), n_months=n,
                        sample_time_indices=times, sample_grid_values=int(a.size),
                        sample_equal=bool(np.array_equal(a, b, equal_nan=True)),
                        sample_max_abs_diff=float(np.nanmax(abs(a-b))),
                        encoded_coordinates_equal=all(np.array_equal(raw[v][:], sub[v][:]) for v in ['time','lat','lon']))
            def dates(ds):
                v = ds['time']
                return [str(t) for t in netCDF4.num2date(v[:], v.units, calendar=getattr(v, 'calendar', 'standard'))]
            item['coordinates_equal'] = (dates(raw) == dates(sub) and
                all(np.array_equal(raw[v][:], sub[v][:]) for v in ['lat','lon']))
            item['matches_committed_task6_input_manifest'] = (
                item['tracking_id'] == input_map[name]['tracking_id'] and item['raw_bytes'] == input_map[name]['bytes'])
        out['paired_zg'].append(item)
        print(name, item['sample_equal'], flush=True)
    for p in sorted((root / 'ACCESS-CM2_data').rglob('*.nc')):
        item = {'file': str(p.relative_to(root)), 'bytes':p.stat().st_size,
                'filesystem_modified_utc': datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat(),
                'validation_scope': 'Header only; does not certify compressed payload integrity or identity to server replacement files.'}
        try:
            with netCDF4.Dataset(p) as ds:
                for key in ['source_id','experiment_id','variant_label','grid_label','creation_date','tracking_id','version']:
                    item[key] = getattr(ds, key, None)
                item['dimensions'] = {k:len(v) for k,v in ds.dimensions.items()}
        except Exception as exc:
            item['header_error'] = str(exc)
        out['tos_pr_headers'].append(item)
    for p in sorted((repo / 'figures/task6_wpsh').glob('*.png')):
        matches = list((root / 'task6_wpsh_outputs').rglob(p.name))
        for local in matches:
            sha = lambda f: hashlib.sha256(f.read_bytes()).hexdigest()
            out['early_task6_figure_matches'].append({'repo_file':str(p.relative_to(repo)),
                'local_file':str(local.relative_to(root)), 'sha256':sha(p), 'identical':sha(p)==sha(local)})
    out['summary'] = {key:all(x[key] for x in out['paired_zg']) for key in [
        'original_sha256_agrees_between_manifests','raw_bytes_match_download_record',
        'sample_equal','coordinates_equal','matches_committed_task6_input_manifest']}
    out['summary'].update(paired_files=len(out['paired_zg']), tos_pr_files=len(out['tos_pr_headers']),
        identical_early_task6_figures=sum(x['identical'] for x in out['early_task6_figure_matches']))
    dest = repo / 'results/report_data_provenance_2026-11-23.json'
    dest.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(out['summary'], indent=2))

if __name__ == '__main__':
    main()
