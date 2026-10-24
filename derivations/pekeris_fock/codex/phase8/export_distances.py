"""Serial pullback and sharp common row denominators over Q[x,y,c]."""
import json
from pathlib import Path
import sympy as S
from rational_panel import H,k,C,encode
OUT=Path(__file__).resolve().parent
x,y,c=S.symbols('x y c');Q=x*x+y*y;e=2*Q-c*c
hv=c*c/e;kv=2*(x*x-y*y)/e
D=json.loads((OUT/'rational_connection.json').read_text());mats=[]
for axis in ['H','k']:
 mats.append({(i,j):S.sympify(v,locals={'H':H,'k':k}) for i,j,v,_ in D['matrices'][axis]['entries']})
keys=set(mats[0])|set(mats[1]);result={}
# Substitute once per distinct coefficient, not once per entry.
cache={v:C(v.subs({H:hv,k:kv},simultaneous=True)) for M in mats for v in M.values()}
for r in [x,y,c]:
 A=S.diff(hv,r);B=S.diff(kv,r)
 M={ij:C(A*cache.get(mats[0].get(ij,0),0)+B*cache.get(mats[1].get(ij,0),0)) for ij in keys}
 result[str(r)]=encode(M,[x,y,c]);print(r,result[str(r)]['nnz'],result[str(r)]['max_degrees'],flush=True)
(OUT/'distance_connection.json').write_text(json.dumps(dict(atoms=D['atoms'],panel=3,field='Q(x,y,c)',matrices=result,panel4='exchange x and y, including the differential coordinate'),indent=2)+'\n')
