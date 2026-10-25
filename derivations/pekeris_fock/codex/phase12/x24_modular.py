"""Independent finite-field Taylor check of every exported numerator cut row."""
import sys
sys.dont_write_bytecode=True
import resource,time,json,math
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
from pathlib import Path
import sympy as S
OUT=Path(__file__).resolve().parent;START=time.time()
cut=json.loads((OUT/'x24_cut.json').read_text());checks={};records=[]
monoms=[(a,d-a) for d in range(4) for a in range(d+1)]
for p in (1000033,1000037):
 assert S.isprime(p)
 def fq(x):
  x=S.Rational(x);return int(x.p)%p*pow(int(x.q)%p,-1,p)%p
 class Jet:
  def __init__(self,x=0):self.c={k:v%p for k,v in x.items() if v%p} if isinstance(x,dict) else ({(0,0):int(x)%p} if int(x)%p else {})
  def __add__(self,o):
   o=o if isinstance(o,Jet) else Jet(o);z=self.c.copy()
   for k,v in o.c.items():z[k]=(z.get(k,0)+v)%p
   return Jet(z)
  __radd__=__add__
  def __neg__(self):return Jet({k:-v for k,v in self.c.items()})
  def __sub__(self,o):return self+-asjet(o)
  def __rsub__(self,o):return asjet(o)+-self
  def __mul__(self,o):
   o=asjet(o);z={}
   for (a,b),v in self.c.items():
    for (c,d),w in o.c.items():
     if a+b+c+d<=3:z[a+c,b+d]=(z.get((a+c,b+d),0)+v*w)%p
   return Jet(z)
  __rmul__=__mul__
  def __pow__(self,n):
   if n<0:return self.inv()**(-n)
   z=Jet(1);b=self
   while n:
    if n&1:z=z*b
    b=b*b;n//=2
   return z
  def inv(self):
   c=self.c.get((0,0),0);assert c
   ci=pow(c,-1,p);e=self*ci-1
   return ci*(1-e+e*e-e*e*e)
  def __truediv__(self,o):return self*asjet(o).inv()
  def __rtruediv__(self,o):return asjet(o)*self.inv()
 def asjet(x):return x if isinstance(x,Jet) else Jet(x)
 r1=Jet({(0,0):fq('1/10'),(1,0):1});r2=Jet({(0,0):fq('9/10'),(0,1):1})
 uv,vv,qq=fq('1/5'),fq('6/5'),fq('1/2')
 n=r1*r1+r2*r2;z=r1*uv+r2*vv
 n0=n.c[0,0];z0=z.c[0,0]
 ratio=(1-n)/(z*z-n);ratio0=ratio.c[0,0];e=ratio*pow(ratio0,-1,p)-1
 # sqrt(ratio/ratio0) = 1+e/2-e^2/8+e^3/16 at total degree 3.
 hh=1+fq('1/2')*e-fq('1/8')*e*e+fq('1/16')*e*e*e
 R3=z*(4*z*z-3*n)
 a1=(4*z*z-n)*(3*n-4)/(12*n**3);c1=-(n-4)*R3/(12*n**3)
 a2=-(4*z*z-n)*(n-4)/(4*n**3);c2=-R3*(3*n-4)/(4*n**3*(n-1))
 al=-(3*n*n*z*z+3*n*n-12*n*z*z-2*n+8*z*z)/(24*n**3);cl=(1-n)*R3/(12*n**3)
 ax=-(-6*n**3*z+15*n*n*z*z+5*n*n-20*n*z*z-2*n+8*z*z)/(120*n**3);cx=(1-n)**2*R3/(60*n**3)
 prim={0:(a2,c2*hh),1:(a1,c1*hh),2:(-al,-cl*hh),3:(ax,cx*hh)}
 syms=S.symbols('n z R1 R2 u v q');loc={str(x):x for x in syms}
 point=dict(zip(syms,[S.Rational(41,50),S.Rational(11,10),S.Rational(1,10),S.Rational(9,10),S.Rational(1,5),S.Rational(6,5),S.Rational(1,2)]))
 def ev(expr):return fq(S.sympify(expr,locals=loc).subs(point))
 derived={}
 for a,b in monoms:
  fac=math.factorial(a)*math.factorial(b)
  got=[val.c.get((a,b),0)*fac%p for val in prim[a+b]]
  expected=[ev(x) for x in cut['monomial_rows'][str(a)+','+str(b)]]
  ok=got==expected;checks[f'p{p}_monomial_{a}_{b}']=ok;assert ok
  derived[a,b]=got
 Y1,Y2=S.symbols('Y1 Y2')
 for idx,item in enumerate(cut['terms']):
  pol=S.Poly(S.sympify(item['four_numerator'],locals={**loc,'Y1':Y1,'Y2':Y2}).subs(syms[-1],S.Rational(1,2)),Y1,Y2)
  got=[sum(fq(co)*derived[powers][j] for powers,co in pol.terms())%p for j in (0,1)]
  expected=[ev(x) for x in item['cut_row']]
  ok=got==expected;checks[f'p{p}_actual_term_{idx}']=ok;assert ok
 records.append({'prime':p,'point':{str(k):str(v) for k,v in point.items()},'monomial_rows':10,'actual_term_rows':8,'mismatches':0})
 print(json.dumps(records[-1]),flush=True)
facts={'checks':checks,'records':records,'method':'degree-three bivariate truncated polynomial jets over Fp, independent of symbolic differentiation; sqrt-state normalized by its constant character',
       'status':'verifies exported local rational rows, not a T5 reduction or cycle pairing','elapsed_seconds':time.time()-START}
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;facts['peak_rss_bytes']=rss if sys.platform=='darwin' else rss*1024
(OUT/'x24_modular.json').write_text(json.dumps(facts,indent=2)+'\n')
print(json.dumps({'checks':len(checks),'elapsed_seconds':facts['elapsed_seconds'],'peak_rss_bytes':facts['peak_rss_bytes']},indent=2))
