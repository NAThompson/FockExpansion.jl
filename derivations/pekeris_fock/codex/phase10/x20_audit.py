"""Read-only input/output inventory, shared prefix verification, and result audit."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,resource,time
root=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
checks={};prefixes={}
for name in ['COLLAB.md','THEORY.md']:
 pre=(out/('x20_prefix_'+name.replace('.md','.txt'))).read_bytes();now=(root/name).read_bytes()
 checks['append_only_'+name]=now.startswith(pre)
 prefixes[name]={'initial_length':len(pre),'current_length':len(now),'initial_sha256':hashlib.sha256(pre).hexdigest(),'current_sha256':hashlib.sha256(now).hexdigest()}
base=json.loads((out/'x20_baseline.json').read_text());changed=[];missing=[]
for rel,old in base.items():
 p=root/rel
 if not p.is_file():missing.append(rel)
 elif sha(p)!=old['sha256']:
  changed.append({'path':rel,'old_sha256':old['sha256'],'new_sha256':sha(p),'ownership':'authorized shared append' if rel in ['COLLAB.md','THEORY.md'] else 'read-only/concurrent artifact; not modified by x20 commands'})
checks['no_baseline_file_removed']=not missing
verify=json.loads((out/'x20_verify.json').read_text());checks['all_scientific_checks']=all(verify['checks'].values())
checks['claimed_ranks']=verify['facts']['summary']['e2']['Picard']==17 and verify['facts']['summary']['e1']['Picard']==18
metrics=[]
for nm in ['x20_counts.json','x20_counts_line.json']:
 x=json.loads((out/nm).read_text());checks[nm+'_under_2GB']=max(x['peak_child_rss_bytes'],x['peak_python_rss_bytes'])<2_000_000_000
 metrics.append({'job':nm,'elapsed_seconds':x['elapsed_seconds'],'peak_child_rss_bytes':x['peak_child_rss_bytes'],'peak_python_rss_bytes':x['peak_python_rss_bytes']})
checks['symbolic_under_2GB']=verify['peak_rss_bytes']<2_000_000_000
checks['only_x20_files']=all(p.name.startswith('x20_') for p in out.iterdir())
report=(out/'x20_report.md').read_text()
checks['research_limits_stated']=all(word in ' '.join(report.split()) for word in ['**OPEN:** nonzero survival','stricter integral definition','Artin–Tate','ψ52','no completed'])
checks['shared_entries_present']=all('TASK_codex_20' in (root/n).read_text() for n in ['COLLAB.md','THEORY.md'])
inputs=['TASK_codex_20.md','codex/phase9/x19_report.md','codex/phase9/x19_geometry.py','codex/phase9/x18_report.md']
inputs += [str(p.relative_to(root)) for p in sorted((root/'claude10').glob('*')) if p.is_file()]
manifest={str(p.relative_to(root)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(out.iterdir()) if p.is_file() and p.name not in ['x20_audit.json','x20_audit.out']}
result={'checks':checks,'scientific_check_count':len(verify['checks']),'input_hashes':{n:sha(root/n) for n in inputs},'shared_prefixes':prefixes,'baseline_changes':changed,'baseline_missing':missing,'outputs':manifest,'jobs':metrics,'symbolic_job':{'elapsed_seconds':verify['elapsed_seconds'],'peak_rss_bytes':verify['peak_rss_bytes']},'write_scope':'x20_* artifacts and append-only COLLAB/THEORY. Claude files read only; concurrent changes attributed separately. Initial compiler TMPDIR was not explicit, so transient external compiler writes are not audited. Subsequent compilation uses in-directory TMPDIR.','agents_used':False,'serial_jobs':True}
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
(out/'x20_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'checks':checks,'scientific_check_count':len(verify['checks']),'baseline_changes':changed,'max_measured_rss_bytes':max([verify['peak_rss_bytes']]+[max(j['peak_child_rss_bytes'],j['peak_python_rss_bytes']) for j in metrics]),'outputs':list(manifest)},indent=2))
