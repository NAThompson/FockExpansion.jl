"""Audit scientific evidence, read-only inputs, append-only prefixes and scope."""
import sys
sys.dont_write_bytecode=True
import ast,hashlib,json,resource,time
from pathlib import Path
START=time.time();OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while chunk:=f.read(1024*1024):h.update(chunk)
 return h.hexdigest()
base=json.loads((OUT/'x24_baseline.json').read_text());checks={};changes=[]
for rel,old in base['files'].items():
 p=ROOT/rel
 if not p.is_file():changes.append({'path':rel,'kind':'missing'});continue
 st=p.stat()
 if 'sha256' in old:
  now=sha(p)
  if now!=old['sha256']:
   changes.append({'path':rel,'kind':'changed','initial_sha256':old['sha256'],'current_sha256':now})
 elif st.st_size!=old['bytes'] or st.st_mtime_ns!=old['mtime_ns']:
  changes.append({'path':rel,'kind':'metadata_changed_large_unhashed_file'})
shared=['COLLAB.md','THEORY.md','NOTES_codex.md']
def protected(rel):
 first=Path(rel).parts[0];name=Path(rel).name
 return first.startswith('claude') or first in ('sym','julia') or name in ('ANSATZ.md','NOTES_claude.md','RESPONSES.md') or name.startswith('REQUEST_')
checks['Claude_owned_snapshots_unchanged']=not any(protected(c['path']) for c in changes)
prefixinfo={}
for name in shared:
 before=(OUT/('x24_prefix_'+name.replace('.md','.txt'))).read_bytes();now=(ROOT/name).read_bytes()
 checks[name+'_append_only']=now.startswith(before)
 checks[name+'_x24_entry_present']=b'TASK_codex_24' in now[len(before):]
 prefixinfo[name]={'initial_bytes':len(before),'current_bytes':len(now),
                   'initial_sha256':hashlib.sha256(before).hexdigest(),'current_sha256':sha(ROOT/name)}
static=['TASK_codex_24.md','codex/phase3/structure.md','codex/phase9/x18_report.md',
 'codex/phase9/x19_report.md','codex/phase10/x20_report.md',
 'codex/phase11/x21_report.md','codex/phase11/x22_report.md',
 'codex/phase11/x22_verify.py','codex/phase11/x22_verify.json']
checks['all_used_static_inputs_unchanged']=all(sha(ROOT/p)==base['files'][p]['sha256'] for p in static)
checks['no_baseline_input_missing']=not any(x['kind']=='missing' for x in changes)
# Record independently changing/harness files without claiming they were ours.
for entry in changes:
 entry['classification']='authorized append-only log' if entry['path'] in shared else 'independent concurrent/harness change; no x24 command writes this path'
known={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and not any(x in p.relative_to(ROOT).parts for x in ('venv','.git','__pycache__'))}
new=sorted(known-set(base['files']))
ours=[p for p in new if p.startswith('codex/phase12/')]
checks['all_new_task_outputs_x24_namespace']=all(Path(p).name.startswith('x24_') for p in ours)
other_new=[p for p in new if not p.startswith('codex/phase12/')]
results={};total=0
for script,data in [('algebra','exact'),('cut','cut'),('numeric','numeric'),('transport','transport'),('q1','q1'),('modular','modular'),('cap_pilot','cap_pilot')]:
 j=json.loads((OUT/('x24_'+data+'.json')).read_text());r=json.loads((OUT/('x24_'+script+'_resource.json')).read_text())
 scientific=j['checks'];total+=len(scientific)
 checks[script+'_all_checks_pass']=bool(scientific) and all(scientific.values())
 checks[script+'_exit_zero']=r['exit_code']==0 and r['reason'] is None
 checks[script+'_stderr_empty']=(OUT/('x24_'+script+'.err')).stat().st_size==0
 checks[script+'_peak_below_2GB']=max(j['peak_rss_bytes'],r['peak_rss_bytes'])<2_000_000_000
 results[script]={'checks':len(scientific),'elapsed_seconds':j['elapsed_seconds'],
                  'peak_rss_bytes':max(j['peak_rss_bytes'],r['peak_rss_bytes'])}
checks['208_recorded_checks']=total==208
checks['207_checks_plus_only_one_cap_finiteness_screen']=total-len(json.loads((OUT/'x24_cap_pilot.json').read_text())['checks'])==207
for p in sorted(OUT.glob('x24_*.py')):
 ast.parse(p.read_text(),filename=str(p));checks[p.name+'_syntax']=True
checks['no_task_pycache_created']=not (OUT/'__pycache__').exists()
report=(OUT/'x24_report.md').read_text()
checks['report_preserves_research_gap']=all(x in report for x in ['**OPEN: the physical T5 pairing','Neither quotient nor middle positions','No new five-state connection or R psi40 witness','no finite engineering bound'])
checks['plain_text_handoffs_present']=all((OUT/n).stat().st_size>500 for n in ['x24_exact.txt','x24_cut.txt','x24_boundaries.txt','x24_q1.txt'])
checks['ledger_records_open_work']='Full limiting physical boundary/nearby-cycle/Gysin map remains OPEN' in (OUT/'x24_ledger.md').read_text()
# All rational export expressions have only documented variables and Rational coefficients.
import sympy as S
j=json.loads((OUT/'x24_cut.json').read_text());allowed=set(S.symbols('n z R1 R2 u v q'))
checks['cut_rows_documented_rational_variables']=all(S.sympify(e).free_symbols<=allowed and S.sympify(e).is_rational_function(*allowed) for item in j['terms'] for e in item['cut_row'])
maxrss=max(v['peak_rss_bytes'] for v in results.values())
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;rss=rss if sys.platform=='darwin' else rss*1024
checks['audit_peak_under_2GB']=rss<2_000_000_000
outputs={p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(OUT.glob('x24_*')) if p.is_file() and not p.name.startswith('x24_audit')}
audit={'checks':checks,'scientific_runs':results,'recorded_check_count':total,'scientific_peak_rss_bytes':maxrss,
 'audit_peak_rss_bytes':rss,'audit_elapsed_seconds':time.time()-START,'shared_prefixes':prefixinfo,
 'baseline_changes':changes,'new_files_outside_task_namespace':other_new,'outputs':outputs,
 'static_inputs':{p:{'sha256':sha(ROOT/p),'bytes':(ROOT/p).stat().st_size} for p in static},
 'scope':'Task commands write codex/phase12/x24_* and append COLLAB.md, THEORY.md, NOTES_codex.md only. Claude-owned files read only. No outside-directory writes requested by any command.',
 'limits':'Hash/command-scope audit, not syscall tracing. Large pre-existing files used metadata only. Independent concurrent/harness changes are recorded, not attributed to this task.',
 'execution':'No agents, concurrent scientific jobs, installs, commits, energies or stencils. Successful scientific runs use one in-process runner with exclusive lock and RSS guard.',
 'mathematical_result':'Exact six-to-four reduction, rational local cap forms, interior physical-germ transport and q=1 S_e2 boundary chart. Full limiting physical T5 map/pairing, rational R psi40 and quotient/middle survival OPEN.'}
(OUT/'x24_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
summary={'audit_checks':len(checks),'all_pass':all(checks.values()),'recorded_scientific_checks':total,'scientific_peak_rss_bytes':maxrss,
 'baseline_changes':changes,'new_files_outside_task_namespace':other_new,'output_count':len(outputs)}
print(json.dumps(summary,indent=2))
assert all(checks.values()),[k for k,v in checks.items() if not v]
