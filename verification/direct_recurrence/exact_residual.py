"""Direct recurrence certificate. All zero tests are exact rational identities."""
import sys,pickle,time
from pathlib import Path
import sympy as s
INPUTS=Path(__file__).with_name('inputs')
sys.path.insert(0,str(INPUTS/'round2'))
import sym_classical as sc
T,U,Z,AL,PS=sc.T,sc.U,sc.Z,sc.AL,sc.PS
pi=s.pi; rt=s.sqrt(2); G=sc.Gc
# Keep Clausen atoms at half-angle scale; only logarithms use duplication.
def cl(m,n,k):
 m,n,k=map(int,(m,n,k));sg=1
 if m==n==0:return s.S(0) if k%8==0 else s.Symbol('Clconst'+str(k%16))
 if (m,n)<(0,0):m,n,k,sg=-m,-n,-k,-1
 k%=16;v=s.Symbol(f'C_{m}_{n}_{k}');sc.REG[v]=('cl',(m,n,k));return sg*v
sc.cl=cl
import sym_H as sh
sh.cl=cl
lam=sc.lam
r1,r2=sc.ca,sc.sa;xi,eta=rt*sc.sp_,rt*sc.cp
sg,dl=r1+r2,r1-r2
a=r1*r1-r2*r2;v=sc.cp**2-sc.sp_**2;sa=2*r1*r2;sb=2*sc.cp*sc.sp_
# Build actual compact formula through the exact adapter.
sys.path.insert(0,str(INPUTS/'round8_classical'))
import compact_classical as cc
class MP:
 pi=pi;catalan=G
 def mpf(self,x):return s.sympify(x)
 def sqrt(self,x):return s.sqrt(x)
 def sin(self,x):
  if x==PS/2:return sc.sp_
  if x==AL/2:return r2
  raise ValueError(x)
 def cos(self,x):
  if x==PS/2:return sc.cp
  if x==AL/2:return r1
  raise ValueError(x)
 def log(self,x):assert x==2;return s.log(2)
cc.mp=MP()
def indices(q):
 q=s.expand(q);return (int(4*q.coeff(AL)),int(4*q.coeff(PS)),int(8*q.subs({AL:0,PS:0})/pi))
def lf(q):return lam(*indices(q))
def cf(q):return cl(*indices(q))
cc.L=lambda q:q*lf(q)+cf(2*q)/2
cc.T=lambda q:-q*lf(q+pi/2)+cf(pi-2*q)/2
cc.lc=lambda q:lf(q+pi/2)
# Rational derivative of log|2sin(q)|, q=m AL/4+n PS/4+k pi/8.
def cotq(m,n,k):
 # exp(i AL/2) = (1+iT)^2/(1+T²), similarly PS.
 # cot q = i (exp(2iq)+1)/(exp(2iq)-1).
 I=s.I;phase=s.expand_complex(s.exp(I*s.pi*s.Rational(k,4)))
 z=phase*((1+I*T)/(1-I*T))**m*((1+I*U)/(1-I*U))**n
 return s.cancel(I*(z+1)/(z-1),extension=[rt,I])
from functools import lru_cache
cotq=lru_cache(None)(cotq)
def scalar_der(e,j):
 x=(T,U)[j];ang=(AL,PS)[j]
 return s.diff(e,x)*(1+x*x)/4+s.diff(e,ang)
def derivative(e,j):
 out=scalar_der(e,j)
 for atom in e.free_symbols.intersection(sc.REG):
  typ,(m,n,k)=sc.REG[atom];mult=s.Rational((m,n)[j],4)
  if not mult:continue
  d=mult*cotq(m,n,k) if typ=='lam' else -mult*lam(m//2,n//2,k//2)
  if typ=='cl':assert m%2==n%2==k%2==0,(m,n,k)
  out+=s.diff(e,atom)*d
 return out
# Angular operator in (AL,PS), with a=cos AL, v=cos PS.
# Lambda²=-4[d_AL²+d_PS²-2 av/(sinAL sinPS)d_AL d_PS+2cotAL d_AL+2cotPS d_PS].
def op(e):
 ea=derivative(e,0);ep=derivative(e,1)
 return -4*(derivative(ea,0)+derivative(ep,1)-2*a*v/(sa*sb)*derivative(ea,1)+2*a/sa*ea+2*v/sb*ep)-21*e
# universal integral action obtained in panel_operator.py; conic evaluations already proved.
from panel_operator import b,y,t,P,h
paneldata=pickle.load(open(Path(__file__).with_name('panels.pkl'),'rb'))
def image(name,bv,yv,ab,be,sinab,sy,s1mb,tanab2):
 Q,R,W,S=paneldata[name];H=s.cancel(S/P);h0=s.expand(H).coeff(t,1);h2=s.expand(H).coeff(t,3)
 rat=s.factor(W-H*h[t]/P)
 D=1+b+(1-b)*t*t
 c1,c2,c3=s.symbols('c1 c2 c3')
 eq=s.Poly(s.cancel((rat-c1*t/(1+t*t)-c2*t/D-c3*t/D**2)*(1+t*t)*D**2),t)
 sol=s.solve(eq.all_coeffs(),(c1,c2,c3),dict=True)[0]
 sub={b:bv,y:yv}; abx=sc.ang(*ab);bex=sc.ang(*be)
 uf=sc.universal(ab,be,s1mb,sy,-yv,tanab2,bv)
 # integral L t/D² = ab/(4 sinab).
 elementary=sol[c1]*pi**2/16+sol[c2]*abx**2/(4*(1-b))+sol[c3]*abx/(4*sinab)
 return (2*h2).subs(sub)*uf['Kinf']+(h0+h2).subs(sub)*uf['K1']+(h0+h2).subs(sub)*(1+2*bex/pi)/(2*sy)*s.log(2)+elementary.subs(sub)/pi
F=s.Rational
panels=[(a,-sc.sp_,(1,0,0),(0,-F(1,2),0),sa,sc.cp,rt*r2,r2/r1),(-a,-sc.sp_,(-1,0,1),(0,-F(1,2),0),sa,sc.cp,rt*r1,r1/r2)]
rot=[(-v,-r1,(0,-1,1),(F(1,2),0,-F(1,2)),sb,r2,rt*sc.cp,sc.cp/sc.sp_),(-v,-r2,(0,-1,1),(-F(1,2),0,0),sb,r1,rt*sc.cp,sc.cp/sc.sp_)]
# classical coefficients separated in Z and state parameters to keep expressions small.
E,A21=s.symbols('E a21')
C=cc.CL_compact(AL,PS,Z)
state=E*(Z*sg*(2+r1*r2)/18-(6-xi*xi)*xi/72)-A21*(Z*sg*v/2-(6-5*xi*xi)*xi/12)
chi=sh.chi()
p10=xi/2-Z*sg
p20=(1-2*E)/12+Z*(chi-sg*xi/3)+Z**2*(F(1,3)+sa/2)+(A21-Z*(pi+4)/(9*pi))*v
p31=Z*(pi-2)/(36*pi)*xi*(5*xi*xi-6)+Z**2*(pi-2)/(6*pi)*sg*v
src=10*p31-2*(1/xi-2*Z*sg/sa)*p20+2*E*p10
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--sector',default='all');ap.add_argument('--reuse',action='store_true');ap.add_argument('--emit-only',action='store_true');args=ap.parse_args()
 for tag,zp,ep,apow in [('Z0',0,0,0),('Z1',1,0,0),('Z2',2,0,0),('Z3',3,0,0),('E',0,1,0),('EZ',1,1,0),('a21',0,0,1),('a21Z',1,0,1)]:
  if args.sector not in ('all',tag):continue
  def coeff(e):return s.expand(e).coeff(Z,zp).coeff(E,ep).coeff(A21,apow)
  raw=Path(__file__).with_name('raw_'+tag+'.pkl')
  if args.reuse and raw.exists():residual=pickle.load(open(raw,'rb'))
  else:
   c=coeff(C+state);print(tag,'differentiate',flush=True)
   residual=op(c)-coeff(src)
   if tag in ('Z1','Z2'):
    for pa in panels:residual+=rt*image('angZ' if zp==1 else 'angZ2',*pa)
    if tag=='Z2':
     for pa in rot:residual+=image('rotZ2',*pa)
  # The table uses twice the K_1 of the conic proof.
  # Clausen duplication: Cl2(2q)=2Cl2(q)+2Cl2(q+pi).
  repl={}
  for atom in residual.free_symbols:
   if str(atom).startswith('C_'):
    m,n,k=map(int,str(atom).split('_')[1:])
    if max(abs(m),abs(n))>=4 and m%2==n%2==k%2==0:
     repl[atom]=2*cl(m//2,n//2,k//2)+2*cl(m//2,n//2,k//2+8)
  residual=residual.xreplace(repl)
  print(tag,'collect',flush=True)
  atoms=sorted(residual.free_symbols.intersection(sc.REG)|{AL,PS},key=str)
  pickle.dump(residual,open(Path(__file__).with_name('raw_'+tag+'.pkl'),'wb'))
  if args.emit_only:
   print(tag,'exact residual generated',flush=True);continue
  aset=set(atoms);zero=(0,)*len(atoms)
  def mul(p,q):
   out={}
   for m,c in p.items():
    for n,d in q.items():
     key=tuple(i+j for i,j in zip(m,n));out[key]=out.get(key,0)+c*d
   return out
  @lru_cache(None)
  def collect(e):
   if not e.free_symbols.intersection(aset):return {zero:e}
   if e in aset:
    m=list(zero);m[atoms.index(e)]=1;return {tuple(m):s.S.One}
   if e.is_Add:
    out={}
    for arg in e.args:
     for m,c in collect(arg).items():out[m]=out.get(m,0)+c
    return out
   if e.is_Mul:
    out={zero:s.S.One}
    for arg in e.args:out=mul(out,collect(arg))
    return out
   if e.is_Pow and e.exp.is_Integer and e.exp>=0:
    out={zero:s.S.One}
    for _ in range(int(e.exp)):out=mul(out,collect(e.base))
    return out
   raise ValueError(e)
  terms=sorted(collect(residual).items(),key=lambda kv: -sum(kv[0]))
  print(tag,len(terms),'atom coefficients',flush=True)
  bad=[]
  for k,(mon,c) in enumerate(terms):
   c=s.cancel(c)
   if c!=0:c=s.cancel(c,extension=rt)
   if c!=0:bad.append((mon,c));print(tag,'NONZERO',mon,str(c)[:350],flush=True)
   elif k%10==0:print(tag,k,'zero',flush=True)
  pickle.dump((atoms,bad),open(Path(__file__).with_name('residual_'+tag+'.pkl'),'wb'))
  assert not bad,(tag,len(bad))
  print(tag,'EXACT ZERO',flush=True)
