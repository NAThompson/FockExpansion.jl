"""Exact rational 55-state panel connection; serial Hermite reduction.
Only writes phase8 artifacts when invoked. No imports from Claude files.
"""
import sympy as S
import json
from pathlib import Path
from functools import lru_cache
s,H,k=S.symbols('s H k'); T=s**4+k*s**2+H; U=s*s+1; V=s*s+H
Ds=[T,U,V]; C=lambda x:S.cancel(x)
kernels=[S.Integer(1),s]+[s**j/T for j in range(4)]+[1/U,s/U,1/V,s/V]
knames=['J','K']+[f'Q{j}' for j in range(4)]+['U0','U1','V0','V1']
jets=[('A',0),('A',1),('A',2),('B',1),('B',2)]
names=[]; spec=[]
for j,nn in enumerate(knames):
 for fam,n in jets:
  if nn=='K' and n==2:continue
  names.append(f'{nn}_{fam}{n}');spec.append((j,fam,n))
names+=['FP_B1']+[f'E_{fam}{n}' for fam,n in jets]+['F0']
assert len(names)==55
idx={v:i for i,v in enumerate(names)}
def atom(j,fam,n):return idx[f'{knames[j]}_{"A" if n==0 else fam}{n}']
def endpoint(fam,n):return idx[f'E_{"A" if n==0 else fam}{n}']
def add(out,i,v):
 if v!=0:out[i]=C(out.get(i,0)+v)
def rem(a,D):return S.rem(a,D,s)
def modfrac(r,D):
 a,b=S.fraction(C(r));return rem(a*S.invert(b,D,s),D)
@lru_cache(None)
def base_reduce(r):
 """r = q' - T'/2T*q + residual with squarefree denominators."""
 r=C(r);q=S.Integer(0)
 for D in Ds:
  # Only double poles occur. modfraction removes the other denominators.
  n=modfrac(C(r*D**2),D)
  if n==0:continue
  mult=S.Rational(3,2) if D==T else S.Integer(1)
  p=rem(-n*S.invert(mult*S.diff(D,s),D,s),D)
  term=C(p/D);q+=term
  r=C(r-S.diff(term,s)+S.diff(T,s)*term/(2*T))
 coeff={}
 for di,D in enumerate(Ds):
  p=modfrac(C(r*D),D)
  if p==0:continue
  r=C(r-p/D)
  for j in range(S.degree(D,s)):
   v=S.expand(p).coeff(s,j)
   if v:coeff[2+j if di==0 else 6+(di-1)*2+j]=v
 p=S.Poly(r,s)
 assert p.degree()<=1,(r,p.degree())
 for j in (0,1):
  if p.nth(j):coeff[j]=p.nth(j)
 return C(q),coeff
@lru_cache(None)
def reduce(r,fam,n):
 q,co=base_reduce(C(r));out={}
 for j,v in co.items():add(out,atom(j,fam,n),v)
 if q!=0:
  add(out,endpoint(fam,n),q.subs(s,1))
  if n==0:add(out,idx['F0'],-q.subs(s,0))
  else:
   angleprime=1/U if fam=='A' else H/V
   for i,v in reduce(C(-q*angleprime),fam,n-1).items():add(out,i,v)
 return {i:v for i,v in out.items() if v!=0}
def build(axis):
 M={}; l=(1/(2*H) if axis==H else 0)-S.diff(T,axis)/(2*T)
 for row,(j,fam,n) in enumerate(spec):
  out=dict(reduce(C(S.diff(kernels[j],axis)+kernels[j]*l),fam,n))
  if fam=='B' and axis==H:
   add(out,row,S.Rational(n,2)/H)
   for i,v in reduce(C(-kernels[j]*s/(2*V)),fam,n-1).items():add(out,i,v)
  for col,v in out.items():
   if v:M[row,col]=v
 row=idx['FP_B1'];fpco=C(l.subs(s,0))
 out=dict(reduce(C((l-fpco)/s),'B',1));add(out,row,fpco)
 if axis==H:
  add(out,row,1/(2*H))
  for i,v in reduce(-1/(2*V),'A',0).items():add(out,i,v)
 for col,v in out.items():
  if v:M[row,col]=v
 el=(1/(2*H) if axis==H else 0)-1/(2*(1+k+H))
 for fam,n in jets:
  row=endpoint(fam,n);M[row,row]=el
  if fam=='B' and axis==H:
   M[row,row]=C(el+S.Rational(n,2)/H)
   M[row,endpoint(fam,n-1)]=-1/(2*(H+1))
 return M

def encode(M,vars):
 rows=[];clears={}
 for (i,j),v in sorted(M.items()):
  v=C(v)
  if not v:continue
  a,b=S.fraction(v);dg=[int(S.Poly(p,*vars).total_degree()) for p in (a,b)]
  rows.append([i,j,str(v),dg]);clears[i]=S.lcm(clears.get(i,S.Integer(1)),b)
 return dict(entries=rows,nnz=len(rows),max_degrees=[max(r[3][n] for r in rows) for n in (0,1)],row_clearing={str(i):str(S.factor(v)) for i,v in clears.items()})
def main():
 mats={}
 for axis in (H,k):
  M=build(axis);mats[str(axis)]=encode(M,[H,k]);print(axis,mats[str(axis)]['nnz'],mats[str(axis)]['max_degrees'],flush=True)
 Path(__file__).with_name('rational_connection.json').write_text(json.dumps(dict(atoms=names,convention='dX=M X; fixed endpoint 1; density sqrt(H)/sqrt(s^4+k*s^2+H); B=sqrt(H)*atan(s/sqrt(H)); divided factorial jets',matrices=mats),indent=2)+'\n')
if __name__=='__main__':main()
