"""Input hashes, append-only prefixes, final verification, output inventory."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, resource, time

start=time.time()
out=Path(__file__).resolve().parent
root=out.parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while chunk:=f.read(1024*1024): h.update(chunk)
    return h.hexdigest()

base=json.loads((out/'x22_baseline.json').read_text())
changes=[]
for rel,old in base.items():
    p=root/rel
    if not p.is_file():
        changes.append({'path':rel,'status':'missing'})
    elif sha(p)!=old['sha256']:
        changes.append({'path':rel,'initial_sha256':old['sha256'],'current_sha256':sha(p),
                        'classification':'authorized append-only shared file' if rel in ['COLLAB.md','THEORY.md'] else
                        'observed concurrent or harness change; no x22 command writes this path'})

inputs=['TASK_codex_22.md','codex/phase3/structure.md','codex/phase9/x18_report.md',
        'codex/phase9/x19_report.md','codex/phase10/x20_report.md',
        'codex/phase11/x21_report.md']
checks={}
checks['no_baseline_files_missing']=not any(c.get('status')=='missing' for c in changes)
checks['task_and_static_scientific_inputs_unchanged']=not any(c['path'] in inputs for c in changes)
prefixes={}
for n in ['COLLAB.md','THEORY.md']:
    before=(out/('x22_prefix_'+n.replace('.md','.txt'))).read_bytes()
    now=(root/n).read_bytes()
    checks[n+'_append_only']=now.startswith(before)
    checks[n+'_x22_entry_present']='TASK_codex_22' in now[len(before):].decode()
    prefixes[n]={'initial_bytes':len(before),'current_bytes':len(now),
                 'initial_sha256':hashlib.sha256(before).hexdigest(),'current_sha256':sha(root/n)}
verify=json.loads((out/'x22_verify.json').read_text())
checks['all_49_scientific_checks']=len(verify['checks'])==49 and all(verify['checks'].values())
checks['scientific_peak_under_2GB']=verify['peak_rss_bytes']<2_000_000_000
checks['final_scientific_stderr_empty']=(out/'x22_verify_final.err').stat().st_size==0
checks['final_log_reports_49_passes']=(out/'x22_verify_final.out').read_text().count('PASS ')==49
checks['all_outputs_x22_namespace']=all(p.name.startswith('x21_') or p.name.startswith('x22_') for p in out.iterdir())
report=(out/'x22_report.md').read_text()
checks['positive_and_open_claims_distinguished']=all(x in report for x in [
    'does **not** contain T5','Step 1 of x21 §5 was attempted',
    'no internal density residue is relabelled','full psi40',
    'This computational filtration need NOT','underlying DIFFERENTIAL',
    'SIX-parameter','Engineering estimates, not benchmarks'])
checks['ledger_completion_present']='Final: all requested fallback deliverables complete' in (out/'x22_plan.md').read_text()

repo=Path('/Users/nathompson7/FockExpansion.jl')
external=['src/feynman_k4.jl','derivations/k4/README.md']
external_hashes={str(repo/n):{'bytes':(repo/n).stat().st_size,'sha256':sha(repo/n)} for n in external}
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if sys.platform!='darwin':rss*=1024
checks['audit_peak_under_2GB']=rss<2_000_000_000
outputs={p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(out.glob('x22_*'))
         if p.is_file() and p.name not in ['x22_audit.json','x22_audit.out','x22_audit.err']}
result={'checks':checks,'baseline_changes':changes,'shared_prefixes':prefixes,
        'scientific_check_count':len(verify['checks']),
        'scientific_elapsed_seconds':verify['elapsed_seconds'],
        'scientific_peak_rss_bytes':verify['peak_rss_bytes'],
        'audit_peak_rss_bytes':rss,'audit_elapsed_seconds':time.time()-start,
        'static_scientific_inputs':{n:{'bytes':(root/n).stat().st_size,'sha256':sha(root/n)} for n in inputs},
        'read_only_repository_inputs_final_hashes':external_hashes,
        'outputs':outputs,'agents_used':False,'concurrent_scientific_jobs':False,
        'initial_hashing_overlap':'Baseline input hashing briefly overlapped an abandoned symbolic primitive exploration; no two scientific jobs ran together.',
        'write_scope':'x22_* in codex/phase11; authorized append-only COLLAB.md and THEORY.md; no Claude-owned or outside-directory writes in any x22 command.',
        'audit_limits':'Hash audit and command inventory, not an OS syscall trace. Concurrent protected outputs can change; used static inputs are checked separately. The abandoned exploration and initial failed JSON-export run did not save their own peak RSS; all successful measured checks are about 70 MB.',
        'source_independent_C3_T5':'PROVED excluded under the cited Hodge inputs, exact algebra VERIFIED',
        'full_density_physical_T5_cut':'OPEN; step 1 not certified, task-authorized fallback used'}
(out/'x22_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'checks':checks,'scientific_check_count':len(verify['checks']),
                  'max_measured_rss_bytes':max(rss,verify['peak_rss_bytes']),
                  'baseline_changes':changes,'output_count':len(outputs)},indent=2))
assert all(checks.values()),[n for n,v in checks.items() if not v]
