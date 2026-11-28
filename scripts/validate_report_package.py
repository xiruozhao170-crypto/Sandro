"""Validate report links, image provenance and inclusion/exclusion rules."""
from pathlib import Path
from urllib.parse import unquote
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]

def main():
    documents = list((ROOT / 'docs').glob('Report*2026-11-23.md'))
    assert len(documents) == 3
    images = []
    links_checked = 0
    for path in documents:
        text = path.read_text(encoding='utf-8')
        assert '\ufffd' not in text and '{{FIG:' not in text
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', text):
            if target.startswith(('https://', 'http://', '#')):
                continue
            target = unquote(target.split('#', 1)[0])
            assert (path.parent / target).is_file(), (path.name, target)
            links_checked += 1
        for target in re.findall(r'!\[[^\]]*\]\(([^)]+)\)', text):
            images.append((path.parent / unquote(target)).resolve())
    manifest = json.loads((ROOT / 'results/report_figure_manifest_2026-11-23.json').read_text(encoding='utf-8'))
    assert manifest['included_research_figures'] == 65
    assert manifest['main_figures'] == 17
    assert manifest['appendix_research_figures'] == 48
    assert len(images) == 66 and len(set(images)) == 66
    labels = set()
    for entry in manifest['figures']:
        p = ROOT / entry['repo_path']
        assert p.resolve() in images, p
        assert hashlib.sha256(p.read_bytes()).hexdigest() == entry['source_sha256'], p
        assert entry['figure_label'] not in labels
        labels.add(entry['figure_label'])
        assert not any(s in p.name.lower() for s in ['regioncheck', 'region_validation', 'membermean', '3membermean'])
    assert {x['source_index'] for x in manifest['excluded']} == {0,21,23,47}
    assert len(list((ROOT/'figures/task1-4/report_miroc6').glob('*.png'))) == 21
    main_text = (ROOT/'docs/Report_2026-11-23.md').read_text(encoding='utf-8')
    assert len(re.findall(r'^\d+\. ',main_text.split('## 参考文献')[1],re.M)) == 22
    result = {'report_documents':len(documents),'research_figures':65,'reference_figures':1,
              'main_figures':17,'appendix_research_figures':48,'references':22,
              'local_links_checked':links_checked,'all_figure_hashes_match':True,
              'inclusion_exclusion_rules_pass':True,'scope':'Artifact integrity; scientific evidence reviewed separately.'}
    (ROOT/'results/report_package_validation_2026-11-23.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
