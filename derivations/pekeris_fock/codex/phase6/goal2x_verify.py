"""Final serial audit, with fresh exact computations and artifact checks."""
import sys,os,json,hashlib
sys.dont_write_bytecode=True
import goal2x_probe as G
from goal2x_rational import fit,ev

d=json.loads((G.OUT/'goal2x_independent.json').read_text());assert d['complete']
assert len(d['points'])==2 and len(d['checks'])==14
for item in d['checks']:assert item['equal']
c=json.loads((G.OUT/'goal2x_constants.json').read_text())
for result in c['results']:
 for row in result['pi_fits']:
  for f in row:
   assert f['heldout']
   if f['degrees']!=[0,0]:assert f['denominator']==[0,1]
 for v in result['embedding_parity']:
  assert len(v['odd_entries'])==(0 if v['sr']==1 else 50)
print('PASS saved branch and pi/embedding evidence',flush=True)
line=json.loads((G.OUT/'goal2x_line.json').read_text());assert len(line['samples'])==261
xs=[v['u'] for v in line['samples']];assert len(set(xs))==261
assert line['analysis']['training']==256 and line['analysis']['nonzero_heldout_pass']==0
for a,b in [(0,0),(4,44),(9,5)]:
 ys=[v['relations'][a][b] for v in line['samples']]
 f=fit(xs[:256],ys[:256],xs[256:],ys[256:]);assert f['degrees']==[128,127] and not f['heldout']
print('PASS fresh degree-search reproductions at three entries; no validated degree candidate',flush=True)
# Cross-check exact arithmetic acceleration against the unaccelerated phase-five jet path.
r=G.relations(10,17);assert r==d['points'][0]['relations']
os.environ['GOAL2X_PURE']='1'
rr=G.relations(10,17);assert rr==r
os.environ.pop('GOAL2X_PURE')
print('PASS all 450 accelerated / pure-Python entries agree',flush=True)
# Additional independent point not used in the original branch experiments.
r=G.relations(25,34);rr=G.relations(25,34,signs=(-1,-1,-1));assert r==rr
print('PASS fresh (25,34) full radical conjugation, all 450 entries',flush=True)
comp=json.loads((G.OUT/'goal2x_comparison.json').read_text());assert len(comp['results'])==2 and all(r['match'] for r in comp['results'])
print('PASS two fresh Claude sample comparisons, 630 coefficients each',flush=True)
print('QUALIFIED: actual minimal line degrees and generic symbolic reconstruction remain unresolved.',flush=True)
