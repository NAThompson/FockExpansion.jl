"""Check reduced fractions, sharp row LCMs and exchange, serially."""
import json,hashlib
from pathlib import Path
import sympy as S
O=Path(__file__).resolve().parent;R=O.parents[1]
x,y,c=S.symbols('x y c');loc=dict(x=x,y=y,c=c)
d=json.loads((O/'distance_connection.json').read_text());rows={};entries=0
for ax,mat in d['matrices'].items():
 lcms={}
 for i,j,v,dg in mat['entries']:
  a,b=S.fraction(S.sympify(v,locals=loc));assert S.Poly(S.gcd(a,b),x,y,c).total_degree()==0
  assert [S.Poly(z,x,y,c).total_degree() for z in (a,b)]==dg
  lcms[i]=S.lcm(lcms.get(i,1),b);rows[i]=S.lcm(rows.get(i,1),b);entries+=1
 for i,v in lcms.items():assert S.cancel(v/S.sympify(mat['row_clearing'][str(i)],locals=loc)).free_symbols==set()
# Record joint smallest polynomial denominator up to a nonzero rational unit.
joint={str(i):dict(atom=d['atoms'][i],factor=str(S.factor(rows.get(i,S.Integer(1)))),degree=int(S.Poly(rows.get(i,1),x,y,c).total_degree())) for i in range(55)}
(O/'row_clearing.json').write_text(json.dumps(joint,indent=2)+'\n')
print('PASS',entries,'reduced distance entries, degrees, axis LCMs; 55 sharp joint row LCMs')
print('joint maximum row degree',max(v['degree'] for v in joint.values()))
print('joint clearing factors',sorted(set(v['factor'] for v in joint.values())))
# Compare generated pullback against symbolic chain rule at rational fresh points.
H,k=S.symbols('H k');src=json.loads((O/'rational_connection.json').read_text());hv=c*c/(2*(x*x+y*y)-c*c);kv=2*(x*x-y*y)/(2*(x*x+y*y)-c*c)
source=[{(i,j):S.sympify(v,locals={'H':H,'k':k}) for i,j,v,_ in src['matrices'][a]['entries']} for a in ['H','k']]
count=0
for xx,yy,cc in [(2,3,4),(3,4,2)]:
 sub={x:S.Rational(xx),y:S.Rational(yy),c:S.Rational(cc)};hk={H:hv.subs(sub),k:kv.subs(sub)}
 for r in [x,y,c]:
  got={(i,j):S.sympify(v,locals=loc).subs(sub) for i,j,v,_ in d['matrices'][str(r)]['entries']}
  for ij in set(source[0])|set(source[1])|set(got):
   expected=S.diff(hv,r).subs(sub)*source[0].get(ij,S.Integer(0)).subs(hk)+S.diff(kv,r).subs(sub)*source[1].get(ij,S.Integer(0)).subs(hk)
   assert S.cancel(got.get(ij,0)-expected)==0;count+=1
print('PASS',count,'fresh exact distance pullback identities (panel4 by electron exchange)')
files=sorted(p for p in (R/'claude8').glob('panel1_*') if p.is_file())
(O/'claude_inventory.json').write_text(json.dumps([dict(path=str(p.relative_to(R)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files],indent=2)+'\n')
print('Claude panel1 artifacts:',len(files),'— conditional cross-check pending' if not files else '— inspect and cross-check required')

outs=json.loads((O/'outputs.json').read_text())
for i,v in outs['panel3_R3psi30']:
 a,b=S.fraction(S.sympify(v,locals=loc));assert S.Poly(a,x,y,c).is_homogeneous and S.Poly(b,x,y,c).is_homogeneous
 assert S.Poly(a,x,y,c).total_degree()-S.Poly(b,x,y,c).total_degree()==3
print('PASS all 22 panel output coefficients have homogeneous degree three')
