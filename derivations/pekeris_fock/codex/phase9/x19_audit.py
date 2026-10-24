"""Final saved-evidence, input-stability and append-only audit; only x19_* writes."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,re
O=Path(__file__).resolve().parent;R=O.parents[1]
sha=lambda b:hashlib.sha256(b).hexdigest()
checks={}
for name in ['x19_check41.json','x19_geometry.json','x19_panelcheck.json']:
 d=json.loads((O/name).read_text())
 if 'checks' in d:assert all(d['checks'].values())
 checks[name]=True
pde=json.loads((O/'x19_check41.json').read_text());assert set(pde['exact_PDE_residuals'])=={'0'} and pde['integer_relations_match']
for n,h in pde['inputs'].items():assert sha((R/n).read_bytes())==h,n
baseline=json.loads((O/'x19_baseline.json').read_text())
inputs=['claude9/closed41.pkl','claude9/coef5_ident.txt','claude9/coef5.txt','claude9/spec40.txt','claude9/spec41.txt','claude9/data2.txt','codex/phase2/order2.py','codex/phase8/x16_summary.json','codex/phase8/x16_clearing.txt','codex/phase8/verify.py','codex/phase8/rational_panel.py']
for n in inputs:
 st=(R/n).stat();assert [st.st_size,st.st_mtime_ns]==baseline[n],('scientific input changed during task',n)
outputs=['x19_check41.out','x19_purity.out','x19_fit40.out','x19_fit40_big.out','x19_fit40_singlet.out','x19_geometry.out','x19_panelcheck.out']
rss=[pde['peak_rss_bytes'],json.loads((O/'x19_geometry.json').read_text())['peak_rss_bytes'],json.loads((O/'x19_panelcheck.json').read_text())['peak_rss_bytes']]
for n in outputs:
 s=(O/n).read_text();assert 'ERROR:' not in s and 'Traceback' not in s,n
 vals=re.findall(r'(?:peak_rss_bytes=|maxrss=)(\d+)',s);rss.extend(map(int,vals))
assert max(rss)<2_000_000_000,rss
for n,field in [('x19_fit40_results.txt',7),('x19_fit40_big_results.txt',6)]:
 rows=[l.split() for l in (O/n).read_text().splitlines() if l and not l.startswith('#')]
 control=[r for r in rows if r[0].startswith('control41')];assert control
 assert max(float(r[field]) for r in control)<(1e-65 if 'big' in n else 1e-9)
sing=[l.split() for l in (O/'x19_fit40_singlet_results.txt').read_text().splitlines() if l and not l.startswith('#')];assert len(sing)==13 and float(sing[-1][5])<1e-35
for p in O.glob('x19_*.py'):compile(p.read_text(),str(p),'exec')
for n,v in json.loads((O/'x19_shared_prefix.json').read_text()).items():assert sha((R/n).read_bytes()[:v['size']])==v['sha256'],('not append-only',n)
current={str(p.relative_to(R)):[p.stat().st_size,p.stat().st_mtime_ns] for p in R.rglob('*') if p.is_file() and not str(p.relative_to(R)).startswith(('codex/phase9/x19_','venv/','.git/'))}
changed=[n for n,v in baseline.items() if n not in ['COLLAB.md','THEORY.md'] and current.get(n)!=v]
new=[n for n in current if n not in baseline]
# Other collaborator activity is observed and attributed only as concurrent, never overwritten.
result=dict(checks=checks,exact_PDE=True,scientific_inputs_stable=True,shared_logs_append_only=True,max_peak_rss_bytes=max(rss),observed_concurrent_existing_changes=changed,observed_concurrent_new_files=new,input_sha256={n:sha((R/n).read_bytes()) for n in inputs},artifact_sha256={p.name:sha(p.read_bytes()) for p in sorted(O.glob('x19_*')) if p.is_file() and p.name!='x19_audit.json'},qualification='Recorded concurrent file changes are observations, not writes made by this task; all task mutations are x19_* or the two authorized shared append-only logs.')
(O/'x19_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['artifact_sha256','input_sha256']},indent=2))
