#!/usr/bin/env python3
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,hashlib,re
P=Path(__file__).resolve().parent
checks={}
for name,token in [('sync4_direct.out','PASS 48 actual'),('sync4_connection.out','PASS 144 exact'),('sync4_weakform.out','PASS manufactured'),('sync4_ritz24.out','N=11'),('sync4_ritz32.out','N=11'),('sync4_a21_refined.out','N=11 a21=0.476748379513083')]:
 assert token in (P/name).read_text(),name
 checks[name]='PASS';print('PASS evidence',name)
s=(P/'sync4_a21.out').read_text()
for N in range(2,12):
 for q in (24,48,96):assert f'bits=192 nq={q} N={N} ' in s,(N,q)
for N in (10,11):assert f'bits=288 nq=48 N={N} ' in s
assert 'bits=192 nq=96 N=11 a21=0.476748379513075' in s
checks['exported_state_checks']=32
print('PASS all 30 angular-order and 2 higher-precision exported-state checks')
for original,digest in json.loads((P/'sync4_input_state_hashes.json').read_text()).items():
 assert hashlib.sha256((P/'input_states'/Path(original).name).read_bytes()).hexdigest()==digest
checks['input_snapshots']='10 SHA256 matches'
print('PASS captured input SHA256 values')
old=json.loads((P/'sync4_protected_hashes.json').read_text())
changed=[k for k,v in old.items() if not (P.parents[1]/k).is_file() or hashlib.sha256((P.parents[1]/k).read_bytes()).hexdigest()!=v]
checks['concurrent_protected_changes']=changed
print('Concurrent protected-file changes (not Codex writes):',changed)
assert set(changed)<=set(['claude2/run_quadcheck.log','claude2/run_quadcheck.out','claude2/states/atoms2_N11.txt']),changed
for f in ['codex_crosscheck2.md','codex_minimality2.md','sync4_counts.md']:assert (P/f).is_file()
for f in ['REVIEW_4_codex.txt','NOTES_codex.md','COLLAB.md']:assert (P.parents[1]/f).is_file()
assert 'Codex position on phase-2 final claims' in (P.parents[1]/'COLLAB.md').read_text()
checks['required_deliverables']='PASS'
(P/'sync4_final_summary.json').write_text(json.dumps(checks,indent=2)+'\n')
print('PASS all required deliverables and recorded verification evidence')
