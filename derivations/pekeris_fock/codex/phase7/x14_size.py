"""Small serial sufficient cover row denominator measurement, no matrix expansion.
Use rational norms to avoid expensive Gaussian-field factorization.
"""
import json
from pathlib import Path
from functools import lru_cache
import sympy as S
OUT=Path(__file__).resolve().parent
h,z=S.symbols('h z',real=True)
@lru_cache(None)
def factors(q):
 p=S.Poly(q,h,z); q=S.expand(q/p.LC())
 try:S.Poly(q,h,z,domain=S.QQ)
 except S.polys.polyerrors.CoercionFailed:q=S.expand(q*S.conjugate(q))
 return dict(S.factor_list(q,h,z)[1])
data=json.loads((OUT/'connections.json').read_text());rows={i:{} for i in range(55)};neighbors={i:set() for i in rows}
for key in ['plus_h','plus_z']:
 for i,j,c,deg in data['cover'][key]['entries']:
  neighbors[i].add(j)
  den=S.fraction(S.sympify(c,locals={'h':h,'z':z}))[1]
  for f,k in factors(den).items():rows[i][f]=max(rows[i].get(f,0),k)
res=[];allf=set()
for i,own in rows.items():
 more=dict(own)
 for j in neighbors[i]:
  for f,k in rows[j].items():more[f]=max(more.get(f,0),k)
 second={f:own.get(f,0)+more.get(f,0) for f in set(own)|set(more)}
 degree=lambda a:int(sum(S.Poly(f,h,z).total_degree()*k for f,k in a.items()))
 res.append({'row':i,'D':str(S.prod(f**k for f,k in own.items())),'degree':degree(own),'D2':str(S.prod(f**k for f,k in second.items())),'second_order_sufficient_degree':degree(second)})
 allf.update(map(str,second))
r={'note':'Sufficient QQ denominators; complex denominator norms can overclear. Neighbor rows included in second-order clearing. Cover only, before physical pullback.', 'cover_rows':res,'max_first_order_row_degree':max(a['degree'] for a in res),'max_second_order_row_degree':max(a['second_order_sufficient_degree'] for a in res),'factors':sorted(allf)}
(OUT/'x14_size.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:v for k,v in r.items() if k!='cover_rows'},indent=2))
