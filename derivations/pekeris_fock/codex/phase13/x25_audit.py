import hashlib, json
from pathlib import Path
out=Path(__file__).resolve().parent; root=out.parents[1]
baseline=json.loads((out/'x25_baseline.json').read_text())
checks={}; details={}
for name,record in baseline.items():
    data=(root/name).read_bytes()
    if name in ['COLLAB.md','THEORY.md']:
        prefix=(out/('x25_prefix_'+name.replace('.md','.txt'))).read_bytes()
        checks[name+'_append_only']=data.startswith(prefix)
        details[name+'_appended_bytes']=len(data)-len(prefix)
    else:
        checks[name+'_unchanged']=hashlib.sha256(data).hexdigest()==record['sha256']
for name in ['exact','periods']:
    data=json.loads((out/('x25_'+name+'.json')).read_text())
    checks[name+'_all_scientific_checks']=all(data['checks'].values())
    checks[name+'_under_2GB']=data['peak_rss_bytes']<2_000_000_000
    details[name+'_check_count']=len(data['checks'])
    details[name+'_peak_rss_bytes']=data['peak_rss_bytes']
checks['all_phase13_artifacts_authorized']=all(p.name.startswith('x25_') for p in out.iterdir())
available=sorted(p.name for p in (root/'claude11').iterdir())
details['claude11_final_inventory']=available
details['missing_audit_inputs']=[n for n in ['NOTES.md','claim_*.md']
                               if not list((root/'claude11').glob(n))]
details['artifact_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(out.glob('x25_*'))
                          if p.is_file() and p.name!='x25_audit.json'}
checks['report_exists']=(out/'x25_report.md').is_file()
checks['summary_log_exists']=(out/'x25_log.md').is_file()
assert all(checks.values()),[k for k,v in checks.items() if not v]
result={'checks':checks,'details':details,
        'scope':'Input hash and append-prefix audit; no git repository available, no OS syscall trace.'}
(out/'x25_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':True,'audit_checks':len(checks),
                  'scientific_checks':details['exact_check_count']+details['periods_check_count'],
                  'max_scientific_rss':max(details['exact_peak_rss_bytes'],details['periods_peak_rss_bytes']),
                  'claude11_available':available,'missing':details['missing_audit_inputs']},indent=2))
