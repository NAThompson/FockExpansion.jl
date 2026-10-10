"""Exact residual arithmetic over Q(sqrt(2))(T,U,pi,G,log(2))."""
import pickle,sys,time
from pathlib import Path
import sympy as s
from flint import fmpq_mpoly_ctx,fmpq
from functools import lru_cache
ctx=fmpq_mpoly_ctx.get(('T','U','R','P','G','L','J'))
T,U,R,PI,G,L,J=[ctx.gen(i) for i in range(7)];mod=R*R-2;modI=J*J+1
one=ctx.constant(1);zero=ctx.constant(0)
def red(p):return (p%mod)%modI
class Rat:
 def __init__(self,n,d=one):
  n=red(n);d=red(d);assert d
  if not n:self.n,self.d=zero,one;return
  g=n.gcd(d);self.n,self.d=n//g,d//g
 def __add__(a,b):
  g=a.d.gcd(b.d);ad=a.d//g;bd=b.d//g
  return Rat(a.n*bd+b.n*ad,ad*b.d)
 def __mul__(a,b):return Rat(a.n*b.n,a.d*b.d)
 def __pow__(a,k):return Rat(a.n**k,a.d**k) if k>=0 else Rat(a.d**(-k),a.n**(-k))
vars={'T':T,'U':U,'G':G}
@lru_cache(None)
def conv(e):
 if e==s.I:return Rat(J)
 if e==s.pi:return Rat(PI)
 if e==s.log(2):return Rat(L)
 if e==s.sqrt(2):return Rat(R)
 if e.is_Rational:return Rat(ctx.constant(fmpq(int(e.p),int(e.q))))
 if e.is_Symbol:return Rat(vars[str(e)])
 if e.is_Add:
  v=Rat(zero)
  for x in e.args:v=v+conv(x)
  return v
 if e.is_Mul:
  v=Rat(one)
  for x in e.args:v=v*conv(x)
  return v
 if e.is_Pow and e.exp.is_Integer:return conv(e.base)**int(e.exp)
 raise ValueError(e)
def getterms(e):
 atoms=sorted([x for x in e.free_symbols if str(x).startswith(('L_','C_')) or str(x) in ('AL','PS')],key=str)
 aset=set(atoms);zero_m=(0,)*len(atoms)
 def mul(p,q):
  out={}
  for m,c in p.items():
   for n,d in q.items():
    key=tuple(i+j for i,j in zip(m,n));out[key]=out.get(key,0)+c*d
  return out
 @lru_cache(None)
 def collect(e):
  if not e.free_symbols.intersection(aset):return {zero_m:e}
  if e in aset:
   m=list(zero_m);m[atoms.index(e)]=1;return {tuple(m):s.S.One}
  if e.is_Add:
   out={}
   for arg in e.args:
    for m,c in collect(arg).items():out[m]=out.get(m,0)+c
   return out
  if e.is_Mul:
   out={zero_m:s.S.One}
   for arg in e.args:out=mul(out,collect(arg))
   return out
  if e.is_Pow and e.exp.is_Integer and e.exp>=0:
   out={zero_m:s.S.One}
   for _ in range(int(e.exp)):out=mul(out,collect(e.base))
   return out
  raise ValueError(e)
 return atoms,sorted(collect(e).items(),key=lambda kv:-sum(kv[0]))
if __name__=='__main__':
 for tag in sys.argv[1:]:
  start=time.time();e=pickle.load(open(Path(__file__).with_name('raw_'+tag+'.pkl'),'rb'))
  atoms,terms=getterms(e);print(tag,len(terms),'coefficients',flush=True);bad=[]
  for i,(mon,c) in enumerate(terms):
   v=conv(c)
   if v.n:bad.append((mon,str(v.n),str(v.d)))
   print(tag,i,'ZERO' if not v.n else 'NONZERO',round(time.time()-start,2),flush=True)
  pickle.dump((atoms,bad),open(Path(__file__).with_name('flint_'+tag+'.pkl'),'wb'))
  assert not bad,(tag,len(bad))
  print(tag,'EXACT ZERO',flush=True)
