"""Final small serial audit; no writes except explicit x14 outputs."""
import sys,json,hashlib,ast,runpy
sys.dont_write_bytecode=True
from pathlib import Path
import sympy as S
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
# Fresh full small-module comparison, not a stale-success check.
runpy.run_path(str(OUT/'x14_check.py'),run_name='__main__')
check=json.loads((OUT/'x14_check.json').read_text())
for p,digest in check['inputs'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==digest
raw=json.loads((ROOT/'claude7/classical_ops.json').read_text())
counts={n:sum(ast.literal_eval(k)[0]==n for k in raw['ops']) for n in ['sigma2P','C1','C2']}
assert counts=={'sigma2P':87,'C1':76,'C2':21},counts
assert len(raw['atoms'])==21
print('PASS saved full operator inventory:',counts, '21 sectors',flush=True)
# Independent exact sufficient denominator spot checks, including difficult E/FP rows.
h,z=S.symbols('h z',real=True);loc={'h':h,'z':z}
dat=json.loads((OUT/'connections.json').read_text());size=json.loads((OUT/'x14_size.json').read_text())
rows={r['row']:r for r in size['cover_rows']};n=0
for ax in ['plus_h','plus_z']:
 for i,j,c,deg in dat['cover'][ax]['entries']:
  if i not in [0,5,8,48,49,50,53,54]:continue
  D=S.sympify(rows[i]['D'],locals=loc);f=S.sympify(c,locals=loc)
  S.Poly(S.cancel(D*f),h,z)
  # Row second-order denominator must cover differentiated own coefficients.
  D2=S.sympify(rows[i]['D2'],locals=loc)
  for t in [h,z]:S.Poly(S.cancel(D2*S.diff(f,t)),h,z)
  n+=1
assert size['max_first_order_row_degree']==12 and size['max_second_order_row_degree']==30
print('PASS cover denominator and differentiated-coefficient checks:',n,flush=True)
assert 2*55*4+2*55*8==1320
assert [int((2*d+1)**3) for d in [30,60,120]]==[226981,1771561,13997521]
print('PASS allocation/cube arithmetic; no expanded distance-degree claim',flush=True)
paths=[OUT/'x14_check.py',OUT/'x14_check.json',OUT/'x14_size.py',OUT/'x14_size.json',OUT/'x14_report.md']
(OUT/'x14_manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2)+'\n')
print('PASS final audit; exact stencil comparison conditional on shard inventory',flush=True)
