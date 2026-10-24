"""Audit actual integral exports, exact raw pullbacks, and sharp reduced row LCMs."""
import ast,json,hashlib
from pathlib import Path
import sympy as S
import mpmath as m
import rational_panel as R
O=Path(__file__).resolve().parent;m.mp.dps=40
node=next(n for n in ast.parse((O/'verify.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='actual');ns={'R':R,'m':m};exec(compile(ast.Module(body=[node],type_ignores=[]),'actual','exec'),ns);actual=ns['actual']
x,y,c=S.symbols('r1 r2 r12');syms=(x,y,c);D=json.loads((O/'rational_connection.json').read_text());summ=json.loads((O/'x16_summary.json').read_text())
base={ax:{(R.names[i],R.names[j]):S.sympify(v) for i,j,v,_ in D['matrices'][ax]['entries']} for ax in ['H','k']}
def readmat(file):
 out={}
 for line in (O/file).read_text().splitlines():
  p,ax,i,j,a,b=line.split('|');out[int(p),ax,i,j]=S.sympify(a)/S.sympify(b)
 return out
raw=readmat('x16_raw55_connection.txt');red=readmat('x16_connection.txt');LC={}
for (p,ax,i,j),v in red.items():LC[p,i]=S.lcm(LC.get((p,i),S.Integer(1)),S.denom(v))
for line in (O/'x16_clearing.txt').read_text().splitlines():
 p,_,i,v=line.split('|');assert S.cancel(S.sympify(v)/LC.get((int(p),i),1)).free_symbols==set()
print('PASS sharp common row LCMs',len(LC),flush=True)
outs={}
for line in (O/'x16_output.txt').read_text().splitlines():
 p,_,i,v=line.split('|');outs[int(p),i]=S.sympify(v)
count=0;worst=m.mpf(0)
for p in range(1,5):
 hv=S.sympify(summ[str(p)]['H']);kv=S.sympify(summ[str(p)]['k'])
 for point in [(S.Rational(13,10),S.Rational(11,10),S.Rational(7,5)),(S.Rational(9,10),S.Rational(7,5),S.Rational(6,5))]:
  sub=dict(zip(syms,point));hh=hv.subs(sub);kk=kv.subs(sub)
  X=actual(m.mpf(str(S.N(hh,42))),m.mpf(str(S.N(kk,42))));vals=dict(zip(R.names,X))
  def num(v):return m.mpf(str(S.N(v.subs(sub),42)))
  # Every raw exported entry checked against an independent exact chain rule.
  for ax in map(str,syms):
   keys={(i,j) for pp,aa,i,j in raw if pp==p and aa==ax};expected=set(base['H'])|set(base['k']);r=S.Symbol(ax)
   for i,j in keys|expected:
    ref=(S.diff(hv,r).subs(sub)*base['H'].get((i,j),S.Integer(0)).subs({R.H:hh,R.k:kk})+S.diff(kv,r).subs(sub)*base['k'].get((i,j),S.Integer(0)).subs({R.H:hh,R.k:kk}))
    assert S.cancel(raw.get((p,ax,i,j),S.Integer(0)).subs(sub)-ref)==0;count+=1
   for i in summ[str(p)]['reachable']:
    a=sum(num(v)*vals[j] for (pp,aa,ii,j),v in raw.items() if (pp,aa,ii)==(p,ax,i))
    b=sum(num(v)*vals[j] for (pp,aa,ii,j),v in red.items() if (pp,aa,ii)==(p,ax,i))
    err=abs(a-b);worst=max(worst,err);assert err<m.mpf('1e-33'),(p,ax,i,err)
  # Direct original physical moment integral; no output-row construction reuse.
  h=m.sqrt(m.mpf(str(S.N(hh,42))));k=m.mpf(str(S.N(kk,42)));b=(h*h-1)/(h*h+1)
  yy=(k/(h*h+1)+1)/2
  def moment(n):
   def fun(s):
    A=m.pi/2-m.atan(s);B=m.atan(s/h);T=(1-b)*s**4+(4*yy-2)*s*s+1+b
    return 2/m.pi*((A*A-B*B) if n==0 else ((h*B/s)**2-(s*A)**2))/m.sqrt(T)
   return m.quad(fun,[0,1])
  xx,zz,cc=map(lambda v:m.mpf(str(v)),point);Q=xx*xx+zz*zz
  if p<3:
   co0=-m.sqrt(2)*2*(18*b*b-36-8*b*b+8)/576
   co2=-m.sqrt(2)*2*(1-b)*(180*yy-90-18*b-32*yy+16)/288
  else:co0=4*(1-b*b)/18;co2=-8*(1-b)*(2*yy-1)/9
  ref=Q**m.mpf('1.5')*(co0*moment(0)+co2*moment(2))
  got=sum(num(v)*vals[i] for (pp,i),v in outs.items() if pp==p)
  err=abs(got-ref);worst=max(worst,err);assert err<m.mpf('1e-33'),(p,err)
 print('PASS panel',p,'two actual points, output and reduced actions',flush=True)
print('PASS',count,'exact raw entry pullbacks; worst actual error',m.nstr(worst,6),flush=True)
files=sorted(O.glob('x16_*'));manifest={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files if f.suffix not in ['.out'] and f.name!='x16_manifest.json'}
(O/'x16_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
