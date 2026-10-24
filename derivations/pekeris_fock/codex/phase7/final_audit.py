from pathlib import Path
import ast,json,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
for f in P.glob('*.py'):ast.parse(f.read_text(),filename=str(f))
for name,rec in json.loads((P/'append_prefixes.json').read_text()).items():
 raw=(ROOT/name).read_bytes()[:rec['bytes']]
 assert hashlib.sha256(raw).hexdigest()==rec['sha256'],name
for name,h in json.loads((P/'input_hashes.json').read_text()).items():
 assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
for log,token in [('connections.out',"classical [('dt', 28, [2, 3])"),('angle_recurrence.out','PASS: 9 exact Laguerre actions'),('verify.out','PASS all serial verification checks'),('fresh_line_check.out','PASS 988 fresh coefficients'),('fixed_endpoint_check.out','PASS 24 pointwise involution checks'),('line_matrices.out','PASS four complete')]:
 assert token in (P/log).read_text(),log
line=json.loads((P/'line_audit.json').read_text())
assert line['training']==128 and line['holdouts']==72
for name,b in line['blocks'].items():assert b['failed_holdouts']==(72 if name=='classical' else 0)
assert sum(b['nonzero'] for name,b in line['blocks'].items() if name!='classical')==438
assert json.loads((P/'fresh_line_check.json').read_text())['coefficients_checked']==988
conn=json.loads((P/'connections.json').read_text());assert len(conn['atoms'])==55
assert conn['cover']['plus_h']['max_degrees']==[6,7]
assert conn['cover']['plus_z']['max_degrees']==[4,5]
assert conn['classical']['dt']['max_degrees']==[2,3]
manifest={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(P.iterdir()) if f.is_file() and f.name not in ['manifest.json','final_audit.out']}
(P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('PASS script syntax; append-only prefixes; unchanged imported inputs; all saved verification logs; degree/count/holdout audits')
print('PASS artifact manifest:',len(manifest),'files; no agents or parallel jobs used')
print('OPEN: minimal/generic ranks; minimal radicals; bivariate reduced rational matrices; sharper unweighted panel estimates; complete elliptic Hamiltonian recurrence and future Claude comparison')
