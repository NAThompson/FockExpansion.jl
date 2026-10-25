import hashlib,json,pathlib
out=pathlib.Path(__file__).resolve().parent
root=out.parent.parent
base=json.loads((out/'x28_inputs.json').read_text())
checks={}
for name,record in base.items():
 data=(root/name).read_bytes()
 if name in ('COLLAB.md','THEORY.md'):
  prefix=(out/('x28_prefix_'+name.replace('.md','.txt'))).read_bytes()
  checks[name+'_prefix_unchanged']=data.startswith(prefix)
  checks[name+'_entry_present']=b'TASK_codex_28' in data[len(prefix):] or b'TASK 28' in data[len(prefix):]
 else:
  checks[name+'_unchanged']=hashlib.sha256(data).hexdigest()==record['sha256']
stdout=(out/'x28_check.out').read_text()
checks['scientific_checks_pass']='CHECKS PASS' in stdout
checks['checker_stderr_empty']=(out/'x28_check.err').stat().st_size==0
rss=int(stdout.split('peak_rss_bytes=')[1].split()[0])
checks['rss_below_2GB']=rss<2_000_000_000
checks['final_summary_in_log']='FINAL SUMMARY:' in (out/'x28_log.md').read_text()
result={'checks':checks,'all_pass':all(checks.values()),'peak_rss_bytes':rss,'scope':'Task commands write x28_* files and append shared logs only; no full-tree before/after baseline was captured.'}
(out/'x28_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
assert result['all_pass']
