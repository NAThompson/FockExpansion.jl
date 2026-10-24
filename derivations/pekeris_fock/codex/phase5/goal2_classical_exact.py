"""Exact finite-field Taylor cyclic-row witnesses; no numerical rank tolerance.
pi is specialized to 101; i and sqrt(2) lie in F_998244353. Witnesses concern
this explicit 18-generator representation, not independence of actual functions.
"""
from pathlib import Path
import json
P=998244353; N=8
I=pow(3,(P-1)//4,P)
# sqrt(2)=(zeta_8+zeta_8^-1), zeta_8=primitive eighth root.
z=pow(3,(P-1)//8,P); RT=(z+pow(z,-1,P))%P
assert I*I%P==P-1 and RT*RT%P==2
IJ=[(a,b) for a in range(N+1) for b in range(N+1-a)]
IX={v:k for k,v in enumerate(IJ)}
PAIRS=[[(IX[(a,b)],IX[(i-a,j-b)]) for a in range(i+1) for b in range(j+1)] for i,j in IJ]
class J:
 def __init__(self,x=0):
  assert not isinstance(x,float),x
  self.v=x if isinstance(x,list) else [int(x)%P]+[0]*(len(IJ)-1)
 def __add__(s,o):
  o=cast(o);return J([(a+b)%P for a,b in zip(s.v,o.v)])
 __radd__=__add__
 def __neg__(s):return J([(-a)%P for a in s.v])
 def __sub__(s,o):return s+-cast(o)
 def __rsub__(s,o):return cast(o)+-s
 def __mul__(s,o):
  if not isinstance(o,J):
   assert not isinstance(o,float),o
   return J([a*int(o)%P for a in s.v])
  return J([sum(s.v[a]*o.v[b] for a,b in pairs)%P for pairs in PAIRS])
 __rmul__=__mul__
 def inv(s):
  v=[0]*len(IJ);v[0]=pow(s.v[0],-1,P)
  for total in range(1,N+1):
   for k,(i,j) in enumerate(IJ):
    if i+j==total:v[k]=-v[0]*sum(s.v[a]*v[b] for a,b in PAIRS[k] if a!=0)%P
  return J(v)
 def __truediv__(s,o):return s*(o.inv() if isinstance(o,J) else pow(int(o)%P,-1,P))
 def __rtruediv__(s,o):return cast(o)*s.inv()
 def __pow__(s,n):
  if n<0:return s.inv()**(-n)
  out=J(1)
  while n:
   if n&1:out=out*s
   s=s*s;n//=2
  return out
 def diff(s,axis):
  v=[]
  for i,j in IJ:
   nxt=(i+1,j) if axis==0 else (i,j+1)
   v.append(((i+1 if axis==0 else j+1)*s.v[IX[nxt]])%P if nxt in IX else 0)
  return J(v)
def cast(v):return v if isinstance(v,J) else J(v)
def setup(t0,q0):
 t=J(t0);t.v[IX[(1,0)]]=1;q=J(q0);q.v[IX[(0,1)]]=1
 # t=exp(i*a/2), q=exp(i*p/2).
 ca=(t+1/t)/2;sa=(t-1/t)/(2*I)
 r1=(ca-sa)/RT;r2=(ca+sa)/RT;sg=r1+r2;d=r1-r2
 xi=RT*(q+1/q)/2;eta0=-RT*(q-1/q)/(2*I)
 pi=101;Z=2
 coeff=[J() for _ in range(18)]; Cs=[]
 arguments=[];idx=0
 for sign in (1,-1):
  eta=sign*eta0; ex=t*q**sign
  for mm in (-1,0,1,2):
   dd,ss=d,sg
   for _ in range(mm%4):dd,ss=-ss,dd
   kz=(dd*(-9*ss**2+6*ss*xi+8*xi**2-8)-eta*(5-2*xi**2))/(144*pi)
   kz2=-5*(xi**2-1)*(dd+eta)/(32*pi)
   nu=xi**2-sg**2 if mm in (0,2) else 9*sg**2+xi**2-2+(-8 if mm==1 else 8)*sg*eta
   nv=d*nu/(36*pi)
   coeff[idx]=(-4 if mm==2 else 4)*(Z*kz+Z**2*(kz2+nv));idx+=1
   arguments.append((ex*pow(I,mm%4,P) if mm!=2 else ex,1,sign,0 if mm==2 else mm*pi,mm==2,4))
  el=Z*(7*sg*d*eta+9*sg**3+8*sg*xi**2-26*sg+10*xi**3-12*xi)/72+Z**2*sg*(137*d*eta+64*sg**2-154*xi**2+26)/288
  coeff[11 if sign==1 else 12]=el
 c9=-eta0*(Z*(2*xi**2-5)/(9*pi)-5*Z**2*(xi**2-1)/(2*pi));coeff[8]=c9/2
 coeff[9]=J(2*Z**2)/(9*pi)*r2*(8*r2**2+xi**2-5)
 coeff[10]=J(2*Z**2)/(9*pi)*r1*(8*r1**2+xi**2-5)
 coeff[15]=Z*xi*(2*xi**2+1)/(144*pi)-5*Z**2*xi*(xi**2-1)/(32*pi)
 coeff[13]=7*Z*xi*d*sg/144+Z**2*(xi*d*sg/64+d*(xi**2-1)/(9*pi))
 coeff[16]=Z**2*sg*(4*sg**2+xi**2-9)/(36*pi)
 coeff[14]=-Z*eta0*(5*xi**2-4)/(36*pi)+Z**2*sg*xi*eta0/(18*pi)
 arguments.extend([(q*q,0,1,0,False,2),(I*t*t,1,0,pi/2,False,2),(I/(t*t),-1,0,pi/2,False,2)])
 # pi/2 above must be exact in field, not a float.
 arguments[-2]=(I*t*t,1,0,pi*pow(2,-1,P)%P,False,2)
 arguments[-1]=(I/(t*t),-1,0,pi*pow(2,-1,P)%P,False,2)
 for axis in (0,1):
  C={}; fac=2/(I*(t if axis==0 else q))
  for k,(ex,aa,pp,cc,is_tan,den) in enumerate(arguments):
   trig=(ex-1)/(I*(ex+1)) if is_tan else I*(ex+1)/(ex-1)
   f=fac*trig*(aa if axis==0 else pp)/den
   for col,val in ((13,aa),(14,pp),(17,cc)):
    if val:C[k,col]=f*val/den
  for k,sgn in ((11,1),(12,-1)):
   ex=t*q**sgn;C[k,17]=-fac*(ex-1)/(I*(ex+1))*(1 if axis==0 else sgn)/4
  C[13 if axis==0 else 14,17]=fac
  C[15 if axis==0 else 16,13 if axis==0 else 14]=2*fac
  Cs.append(C)
 return coeff,Cs

def rankadd(basis,row):
 row=[v%P for v in row]
 for p,v in basis:
  f=row[p];row=[(a-f*b)%P for a,b in zip(row,v)]
 pivot=next((i for i,v in enumerate(row) if v),None)
 if pivot is None:return False
 row=[v*pow(row[pivot],-1,P)%P for v in row]
 basis.append((pivot,row));return True

def run(t,q):
 l,C=setup(t,q);rows={(0,0):l};basis=[];selected=[];history=[]
 for order in range(N+1):
  for a in range(order+1):
   b=order-a
   if order:
    axis=0 if a else 1;parent=(a-1,b) if a else (a,b-1)
    old=rows[parent];out=[v.diff(axis) for v in old]
    for (i,j),c in C[axis].items():out[j]=out[j]+old[i]*c
    rows[a,b]=out
   if rankadd(basis,[v.v[0] for v in rows[a,b]]):selected.append((a,b))
  history.append(len(basis));print((t,q),order,len(basis),flush=True)
 result=dict(point=[t,q],ranks=history,selected=selected)
 return result
if __name__=='__main__':
 results=[run(2,3),run(5,7)]
 Path(__file__).with_name('goal2_classical_exact.json').write_text(json.dumps(dict(prime=P,i=I,sqrt2=RT,pi=101,N=N,results=results),indent=2)+'\n')
