from pathlib import Path
import json,hashlib
root=Path.cwd(); out=root/'codex/phase13'
base=json.loads((out/'x29_audit.json').read_text()); checks={}
for name,rec in base.items():
    data=(root/name).read_bytes()
    if name in ['COLLAB.md','THEORY.md']:
        data=data[:rec['bytes']]
    checks[name]=hashlib.sha256(data).hexdigest()==rec['sha256']
for name in ['x29_exact.json','x29_guess.json','x29_controls.json','x29_samples.json']:
    z=json.loads((out/name).read_text()); checks[name]=bool(z)
e=json.loads((out/'x29_exact.json').read_text()); g=json.loads((out/'x29_guess.json').read_text()); c=json.loads((out/'x29_controls.json').read_text())
checks['minimal_order_degree']=e['order']==5 and e['degree']==22 and e['exact_cyclic_rank']==5
checks['exact_reconstruction']=g['records'][-1]['rational_reconstruction_matches'] is True
checks['210_digits']=g['dps']==210
checks['automatic_left_E']=c['automatic_factor_orders']==[2,1,1,1] and c['automatic_hits'][0][0:2]==[0,'E']
checks['resource_recorded_runs']=max(e['rss_bytes'],g['peak_rss_bytes'],c['peak_rss_bytes'])<1500000000
checks['final_summary_line']=(out/'x29_log.md').read_text().splitlines()[-1].startswith('FINAL SUMMARY:')
assert all(checks.values()),checks
artifacts={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(out.glob('x29_*')) if p.is_file() and p.name!='x29_audit_result.json'}
(out/'x29_audit_result.json').write_text(json.dumps({'checks':checks,'artifacts':artifacts},indent=2)+'\n')
print('AUDIT PASS',len(checks),'checks;',len(artifacts),'artifacts; shared-log prefixes intact')
