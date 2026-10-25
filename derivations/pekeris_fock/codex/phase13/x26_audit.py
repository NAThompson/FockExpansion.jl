"""Final ownership, arithmetic and reported-exception audit."""
import json,hashlib,re
from decimal import Decimal as D,getcontext
from pathlib import Path
getcontext().prec=90
p=Path(__file__).resolve().parent;r=p.parents[1]
baseline=json.loads((p/'x26_baseline.json').read_text())
checks={};detail={};failures=[]
for name,b in baseline.items():
 data=(r/name).read_bytes()
 if name in ('COLLAB.md','THEORY.md'):
  checks[name+'_append_only']=hashlib.sha256(data[:b['length']]).hexdigest()==b['sha256']
  detail[name+'_appended_bytes']=len(data)-b['length']
 else:
  if hashlib.sha256(data).hexdigest()!=b:failures.append(name)
detail['concurrent_changed_original_files']=failures
checks['protected_original_files_unchanged_except_observed_concurrent_activity']=set(failures)<= {'codex_26.log','claude11/chain.jl','claude11/guess_ode.jl'}
new=[str(f.relative_to(r)) for f in r.rglob('*') if f.is_file() and str(f.relative_to(r)) not in baseline]
detail['concurrent_new_files_outside_x26']=[s for s in new if not s.startswith('codex/phase13/x26_')]
checks['new_files_outside_x26_only_in_active_claude11']=all(s.startswith('codex/phase13/x26_') or s.startswith('claude11/') for s in new)
detail['new_files']=new
exact=json.loads((p/'x26_exact.json').read_text())
checks['exact_all_pass']=all(exact['checks'].values());detail['exact_check_count']=len(exact['checks'])
metrics={k:D(v) for k,v in (line.split('\t') for line in (p/'x26_metrics.tsv').read_text().splitlines())}
checks['corrected_double_fit_below_1e55']=metrics['corrected_double_fit_max']<D('1e-55')
checks['corrected_double_other_below_1e55']=metrics['corrected_double_other_sectors']<D('1e-55')
checks['triple_stored_zero']=metrics['corrected_triple_all']==0
checks['corrected_matches_extend']=max(v for k,v in metrics.items() if k.startswith('corrected_vs_extend'))<D('1e-65')
checks['duplicate_psi41_demonstrated']=metrics['psi41_original_minus_twice_corrected']<D('1e-44')
checks['duplicate_psi40_effect_nonzero']=metrics['psi40_first_log_original_minus_corrected']>D('.009')
checks['corrected_sampled_convergence']=max(v for k,v in metrics.items() if k.startswith('n64_n80_step0'))<D('1e-34')
mode=json.loads((p/'x26_mode.json').read_text())
positive=D(mode['sin_positive_log_coefficient'])
legacy=(p/'x26_saved.out').read_text()
saved=D(re.search(r'psi20 Z1 code A1 projection onto sin\(Y/2\): ([^\n]+)',legacy).group(1))
checks['independent_flux_magnitude_and_sign']=abs(positive+saved)<D('1e-42')
mon=json.loads((p/'x26_monodromy.json').read_text())
checks['reflection_involution']=D(mon['R_u_minus']['involution_residual'])<D('1e-17')
checks['v1_unipotent_square']=D(mon['U_v_plus']['nilpotent_square_residual'])<D('1e-17')
checks['actual_T5_shared_line_control']=D(mon['conjugation_R']['stacked_images'][2])>1 and D(mon['conjugation_R']['stacked_images'][3])<D('1e-18') and D(mon['conjugation_R']['product'][0])>100 and D(mon['conjugation_R']['triple'][0])<D('1e-15')
checks['actual_T5_v1_nonblind_control']=D(mon['conjugation_U']['stacked_images'][3])>D('.3') and D(mon['conjugation_U']['triple'][1])>100
TA=json.loads((p/'x26_TA.json').read_text())
checks['opposite_circle_nonblind_control']=D(TA['images_concatenated_singular_values'][3])>D('.1') and D(TA['Nminus_Nplus_Nminus_singular_values'][1])>20
for stem in ('x26_chain80_bounded','x26_verify_saved'):
 text=(p/(stem+'.out')).read_text();rss=int(re.search(r'x26 peak_rss_bytes (\d+)',text).group(1))
 checks[stem+'_under2GB']=rss<2_000_000_000;detail[stem+'_peak_rss_bytes']=rss
checks['bounded_solver_no_stderr']=not (p/'x26_chain80_bounded.err').read_text()
checks['verifier_no_stderr']=not (p/'x26_verify_saved.err').read_text()
report=(p/'x26_report.md').read_text();log=(p/'x26_log.md').read_text()
checks['resource_exception_disclosed']='2,518,417,408' in report and '2,518,417,408' in log
checks['remaining_open_questions_disclosed']='actual W_a versus transported W_b' in report and 'attribution' in report
checks['summary_log_present']='**Final summary:**' in log
exceptions=[{'constraint':'all runs <2 GB','status':'VIOLATED','recorded_peak_bytes':2518417408,'earlier_peaks':'unknown','mitigation':'interrupted pending repeat; bounded n80 rerun and final verifier below 2 GB'}]
result={'checks':checks,'details':detail,'known_exceptions':exceptions,'workspace_change_attribution':'Observed concurrent Claude11 activity and harness log updates; ownership is not inferred from a global unchanged-files assertion. Our explicit write commands target only x26_* and append-only shared logs.','full_task_solved':False,'research_gaps':['actual based T5 strip transport','term-by-term original-source derivation of rational coefficients'],'artifact_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in p.glob('x26_*') if f.is_file() and f.name not in ('x26_audit.json','x26_audit.out')}}
(p/'x26_audit.json').write_text(json.dumps(result,indent=2))
assert all(checks.values()),[k for k,v in checks.items() if not v]
print(json.dumps({'audit_checks_passed':len(checks),'exact_checks_passed':len(exact['checks']),'known_exceptions':exceptions,'full_task_solved':False,'details':{k:v for k,v in detail.items() if 'rss' in k}},indent=2))
