"""Read-only input/prefix checks, output manifest and measured-resource audit."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib, json, resource, time

out=Path(__file__).resolve().parent
root=out.parents[1]
start=time.time()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while chunk:=f.read(1024*1024): h.update(chunk)
    return h.hexdigest()

checks={}
base=json.loads((out/'x21_baseline.json').read_text())
changes=[]
for rel,old in base.items():
    p=root/rel
    if not p.is_file():
        changes.append({'path':rel,'change':'missing'})
    elif sha(p) != old['sha256']:
        changes.append({'path':rel,'initial_sha256':old['sha256'],'current_sha256':sha(p),
                        'classification':'authorized append-only shared file' if rel in ('COLLAB.md','THEORY.md') else
                        'observed concurrent/input change; no x21 tool command writes this path'})
checks['baseline_files_present']=all(c.get('change')!='missing' for c in changes)
protected_changes=[c for c in changes if c['path'].startswith(('claude','sym/')) or
    c['path'] in ('ANSATZ.md','NOTES_claude.md','RESPONSES.md') or c['path'].startswith('REQUEST_')]
prefixes={}
for n in ('COLLAB.md','THEORY.md'):
    prefix=(out/('x21_prefix_'+n.replace('.md','.txt'))).read_bytes()
    current=(root/n).read_bytes()
    checks['append_only_'+n]=current.startswith(prefix)
    checks['entry_present_'+n]='TASK_codex_21' in current[len(prefix):].decode()
    prefixes[n]={'initial_bytes':len(prefix),'current_bytes':len(current),
                 'initial_sha256':hashlib.sha256(prefix).hexdigest(),'current_sha256':sha(root/n)}
verify=json.loads((out/'x21_verify.json').read_text())
checks['all_25_scientific_checks']=len(verify['checks'])==25 and all(verify['checks'].values())
checks['scientific_peak_under_2GB']=verify['peak_rss_bytes']<2_000_000_000
checks['checker_stderr_empty']=(out/'x21_verify.err').stat().st_size==0
checks['only_x21_outputs']=all(p.name.startswith('x21_') for p in out.iterdir())
report=(out/'x21_report.md').read_text()
checks['research_limits_explicit']=all(t in report for t in (
    '**OPEN: Conjecture A','complete weight-graded relative','not yet the requested transcendental-cycle pairing',
    'Pm=0 AND I m=0','3+3-5=1','currently uncosted part'))
checks['all_plan_steps_closed']=all(line.startswith('- [x]') for line in (out/'x21_plan.md').read_text().splitlines() if line.startswith('- ['))
inputs=['TASK_codex_21.md','codex/phase10/x20_report.md','codex/phase9/x19_report.md',
        'codex/phase9/x18_report.md','codex/phase3/structure.md','codex/phase7/global_and_domain.md',
        'claude9/spec40.txt','claude9/psi40_big.jl','claude9/fock_spectral_copy.jl',
        'claude10/op_e2_Q.sobj','claude10/gd.sage','claude10/analyze.sage']
checks['protected_scientific_inputs_unchanged']=not any(c['path'] in inputs for c in protected_changes)
outputs={p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(out.iterdir())
         if p.is_file() and p.name not in ('x21_audit.json','x21_audit.out')}
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if sys.platform!='darwin':rss*=1024
checks['audit_peak_under_2GB']=rss<2_000_000_000
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
result={'checks':checks,'scientific_check_count':len(verify['checks']),
        'scientific_elapsed_seconds':verify['elapsed_seconds'],'scientific_peak_rss_bytes':verify['peak_rss_bytes'],
        'audit_peak_rss_bytes':rss,'audit_elapsed_seconds':time.time()-start,
        'shared_prefixes':prefixes,'baseline_changes':changes,
        'observed_concurrent_protected_output_changes':protected_changes,
        'inputs':{n:{'bytes':(root/n).stat().st_size,'sha256':sha(root/n)} for n in inputs},
        'outputs':outputs,'serial_jobs':True,'agents_used':False,
        'write_scope':'x21_* in codex/phase11; append-only COLLAB.md and THEORY.md Log. No outside-directory or Claude-owned writes requested by any x21 command.',
        'limits':'Hash audit covers root files and claude9/claude10/sym inputs, not every multi-GB historical Claude output. It does not claim an OS syscall trace. The first optional OS time wrapper returned a sandbox sysctl error; the final direct checker exits successfully with empty stderr.',
        'conjecture_A':'OPEN','full_relative_GrW':'OPEN'}
(out/'x21_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'checks':checks,'scientific_check_count':len(verify['checks']),
                  'max_measured_rss_bytes':max(rss,verify['peak_rss_bytes']),
                  'baseline_changes':changes,'outputs':list(outputs),
                  'conjecture_A':'OPEN','full_relative_GrW':'OPEN'},indent=2))
