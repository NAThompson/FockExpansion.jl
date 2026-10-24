"""Exact cyclic-row experiments on an algebraic finite-field chart.
Sparse represented states; only adjacent derivative orders retained.
Includes the known zero-jet and endpoint quotient. No numerical rank threshold.
"""
import sys,json
from pathlib import Path
from types import SimpleNamespace
import sympy as s
import goal2_classical_exact as A
# One extra jet order makes differentiated pullback coefficients valid.
A.N=15
A.IJ=[(i,k-i) for k in range(A.N+1) for i in range(k+1)]
A.IX={v:k for k,v in enumerate(A.IJ)}
A.PAIRS=[[(A.IX[(a,b)],A.IX[(i-a,j-b)]) for a in range(i+1) for b in range(j+1)] for i,j in A.IJ]
P=A.P;I=A.I;RT=A.RT;J=A.J;WORK=A.N
oldmul=J.__mul__
def mul(x,y):
 if isinstance(y,Lin):return NotImplemented
 if not isinstance(y,J):return oldmul(x,y)
 return J([sum(x.v[a]*y.v[b] for a,b in pairs)%P for pairs in A.PAIRS[:(WORK+1)*(WORK+2)//2]])
J.__mul__=mul;J.__rmul__=mul
oldadd=J.__add__
def addj(x,y):
 if isinstance(y,Lin):return NotImplemented
 return oldadd(x,y)
J.__add__=addj;J.__radd__=addj

def diffj(x,axis):
 out=[]
 for i,j in A.IJ[:(WORK+1)*(WORK+2)//2]:
  nx=(i+1,j) if axis==0 else (i,j+1);idx=A.IX.get(nx,len(x.v))
  out.append((i+1 if axis==0 else j+1)*x.v[idx]%P if idx<len(x.v) else 0)
 return J(out)
J.diff=diffj

def sqrtj(x):
 roots=s.sqrt_mod(x.v[0],P,all_roots=True)
 if not roots:raise ValueError('nonsplit algebraic chart')
 root=int(roots[0]);delta=x/x.v[0]-1;power=J(1);out=J(1);coef=1
 for k in range(1,A.N+1):
  power=power*delta;coef=coef*(3-2*k)*pow(2*k,-1,P)%P;out+=coef*power
 return root*out

class Lin:
 def __init__(self,d=None):self.d=d or {}
 def __add__(self,o):
  if o==0:return self
  assert isinstance(o,Lin),type(o)
  d=dict(self.d)
  for k,v in o.d.items():d[k]=d.get(k,J())+v
  return Lin(d)
 __radd__=__add__
 def __neg__(self):return Lin({k:-v for k,v in self.d.items()})
 def __sub__(self,o):return self+-o
 def __mul__(self,o):
  o=A.cast(o);return Lin({k:v*o for k,v in self.d.items()})
 __rmul__=__mul__
 def __truediv__(self,o):return self*A.cast(o).inv()
 def __eq__(self,o):return o==0 and not self.d if not isinstance(o,Lin) else self.d==o.d

CONTEXT=None
def geometry(b,y):return CONTEXT[:3]
def data(b,y,axis,family):
 h,z,u,side,f0=CONTEXT
 hb=h/(1-b*b);dz=[]
 for q in z[:4]:
  tq=4*(1-b)*q**3+2*(4*y*y-2)*q
  dz.append(-(1-q**4 if axis==0 else 8*y*q*q)/tq)
 dz += [J(),J(),I*hb if axis==0 else J(),-I*hb if axis==0 else J()]
 du=hb if axis==0 and side<0 else J()
 aa=[[-J(1)/2,J(),J()] for _ in range(4)]
 tj=[J(),J(1),J()] if family=='A' else [J(),J(),J()]
 sj=[J(),J(1),J()] if family=='D' else [J(),J(),J()]
 aa += [[v/(2*I) for v in tj],[-v/(2*I) for v in tj],[v/(2*I) for v in sj],[-v/(2*I) for v in sj]]
 l0=-1/(2*(1+b)) if axis==0 else J();kap=[l0,J(),J()]
 for ai,zi,vi in zip(aa,z,dz):kap=[v-w*(vi/zi) for v,w in zip(kap,ai)]
 return h,z,u,dz,du,aa,l0,kap
src=(Path(__file__).resolve().parents[1]/'phase3/panel_connection.py').read_text()
helpers=src[src.index('def add'):src.index('def geometry')].replace('m.mpf(0)','J(0)')
body=src[src.index('def rhs'):src.index('def direct_moments')].replace('m.pi','J(101)')
ns=dict(J=J,m=SimpleNamespace(re=lambda x:x),geometry=geometry,data=data)
exec(helpers+body,ns)


def context(t,q):
 ca=(t+1/t)/2;sa=(t-1/t)/(2*I)
 r=(ca-sa)/RT;rr=(ca+sa)/RT;c=RT*(q+1/q)/2;pp=-RT*(q-1/q)/(2*I)
 D=(1-c*c)**2-4*r*r*rr*rr;disc=sqrtj(D)
 z1=sqrtj((1-c*c+disc)/(2*rr*rr));z3=sqrtj((rr*rr-r*r+disc)/(pp*pp))
 zs=[z1,rr/r*z1,z3,I*c/(pp*z3)];hs=[r/rr,rr/r,c/pp,c/pp]
 sin_a=(t*t-t**-2)/(2*I);cos_a=(t*t+t**-2)/2
 sin_p=(q*q-q**-2)/(2*I);cos_p=(q*q+q**-2)/2
 maps=[(-sin_a,-(q+1/q)/2),(sin_a,-(q+1/q)/2),(cos_p,-r),(cos_p,-rr)]
 contexts=[]
 for j in range(4):
  h=hs[j];z=zs[j];side=-1 if j==1 else 1;u=h if side<0 else J(1)
  roots=[z,-z,h/z,-h/z,J(I),J(-I),I*h,-I*h]
  contexts.append((h,roots,u,side,1/(RT*(r if j==0 else rr)) if j<2 else 1/c))
 return r,rr,c,pp,sin_a,cos_a,sin_p,cos_p,maps,contexts


def build(t0,q0,include_pieces=False):
 global CONTEXT,WORK
 WORK=A.N
 coeff,Cclassic=A.setup(t0,q0)
 t=J(t0);t.v[A.IX[(1,0)]]=1;q=J(q0);q.v[A.IX[(0,1)]]=1
 r,rr,c,pp,sa,ca,sp,cp,maps,contexts=context(t,q)
 pieces=[dict(name="classical",rows=[{j:v for j,v in enumerate(coeff)}])]
 C=[dict(cc) for cc in Cclassic];L=[{17:J(1)},{330:-2*J(101)*(-cp),14:-2*c*pp,331:-sa},{j:v for j,v in enumerate(coeff)}]
 proj={i:{i:J(1)} for i in range(18)}
 X=[[[Lin({i*3+k:J(1)}) for k in range(3)] for i in range(13)] for _ in range(2)]
 for j,((b,y),ctx) in enumerate(zip(maps,contexts)):
  CONTEXT=ctx;off=18+j*78;coff=18+j*55
  for fam in range(2):
   for axis in range(2):
    out=ns['rhs'](b,y,axis,'A' if fam==0 else 'D',X[fam])
    for i,row in enumerate(out):
     for k,lin in enumerate(row):
      srcidx=off+fam*39+i*3+k
      for col,v in lin.d.items():
       for ax in range(2):
        chain=b.diff(ax) if axis==0 else y.diff(ax)
        key=(srcidx,off+fam*39+col)
        C[ax][key]=C[ax].get(key,J())+v*chain
  XA=[[[Lin({off+fam*39+i*3+k:J(1)}) for k in range(3)] for i in range(13)] for fam in range(2)]
  n0,n2=ns['moments_from_state'](b,y,*XA)
  pieces.append(dict(name=f'panel{j+1}',rows=[n0.d,n2.d]))
  Z=2
  if j<2:
   c0=-J(RT)*Z*(9*Z*b*b-18*Z-8*b*b+8)/576
   c2=-J(RT)*Z*(1-b)*(90*Z*y*y-45*Z-18*b-32*y*y+16)/288
  else:c0=Z*Z*(1-b*b)/18;c2=-2*Z*Z*(1-b)*(2*y*y-1)/9
  for k,v in (c0*n0+c2*n2).d.items():L[2][k]=v
  # Actual zero jets and physical endpoint relations; target dimension 240.
  for fam in range(2):
   for i in range(11):
    for k in range(3):
     col=coff+i*5+(0 if k==0 else (k if fam==0 else k+2))
     proj[off+fam*39+i*3+k]={col:J(1)}
  h,z,u,side,f0=ctx;e0=1/(-2*y*u)
  if j<2:
   angles=[(J(101)/4,J(),J()),(J(101)/4,J(1)/2,J())] if j==0 else [(J(101)/4,J(1)/2,J()),(J(101)/4,J(),J())]
  else:angles=[(J(101)/4,J(),J()),(J(),J(),-J(1)/2)]
  for fam,(cc,aa,ap) in enumerate(angles):
   e=off+fam*39+33;f=off+fam*39+36
   proj[e]={17:e0};proj[e+1]={17:e0*cc,13:e0*aa,14:e0*ap}
   proj[e+2]={17:e0*cc*cc/2,13:e0*cc*aa,14:e0*cc*ap,15:e0*aa*aa/2,16:e0*ap*ap/2}
   proj[f]={17:f0};proj[f+1]={};proj[f+2]={}
 # Old L,U connection in t,q coordinates.
 S=r+rr;V=-cp;sn=ca;cs=-sa;g=c*pp;den=sn*sn-V*V;az=cs+S*c
 for ax in range(2):
  ds=S.diff(ax);dc=c.diff(ax);dp=pp.diff(ax);dsn=sn.diff(ax);dcs=cs.diff(ax);dv=V.diff(ax)
  kg=dc/c+dp/pp
  j0=2*(dcs+ds*c+S*dc)/az-dsn/(1+sn)-(g*g*kg-cs*dcs)/den
  tx=2*(sn*dv-V*dsn)/den;tb=2*c*c*(cs*kg-dcs)/den
  C[ax][330,17]=(ds+dc)/(S+c)
  C[ax][331,17]=J(101)*j0+J(101)/2*tx
  C[ax][331,13]=tx;C[ax][331,14]=-pp/c*tb
 proj[330]={238:J(1)};proj[331]={239:J(1)}
 # Remove identically zero edges; store by source for sparse row action.
 edges=[]
 for cc in C:
  ee={}
  for (i,j),v in cc.items():
   if any(v.v):ee.setdefault(i,[]).append((j,v))
  edges.append(ee)
 return (L,edges,proj,pieces) if include_pieces else (L,edges,proj)


def project(row,proj):
 out=[0]*240
 for i,v in row.items():
  for j,w in proj.get(i,{}).items():out[j]=(out[j]+v.v[0]*w.v[0])%P
 return out


def flatcheck(C):
 curvature={}
 for i,ee in C[1].items():
  for j,v in ee:curvature[i,j]=curvature.get((i,j),0)+v.v[A.IX[(1,0)]]
 for i,ee in C[0].items():
  for j,v in ee:curvature[i,j]=curvature.get((i,j),0)-v.v[A.IX[(0,1)]]
 for first,second,sgn in ((C[1],C[0],1),(C[0],C[1],-1)):
  for i,ee in first.items():
   for k,v in ee:
    for j,w in second.get(k,[]):curvature[i,j]=curvature.get((i,j),0)+sgn*v.v[0]*w.v[0]
 bad={ij:v%P for ij,v in curvature.items() if v%P}
 assert not bad,('nonflat exact port',list(bad.items())[:5])
 print('PASS exact constant curvature on full 332-state presentation',flush=True)

def run(t,q,D=14):
 global WORK
 L,C,T=build(t,q);flatcheck(C);previous={(0,0):L};basis=[];selected=[];hist=[]
 for order in range(D+1):
  WORK=A.N-order;current={}
  for a in range(order+1):
   b=order-a
   if not order:rows=L
   else:
    axis=0 if a else 1;old=previous[a-1,b] if a else previous[a,b-1];rows=[]
    for row in old:
     out={i:v.diff(axis) for i,v in row.items()}
     for i,v in row.items():
      for j,w in C[axis].get(i,[]):out[j]=out.get(j,J())+v*w
     rows.append(out)
   current[a,b]=rows
   for src,row in enumerate(rows):
    if A.rankadd(basis,project(row,T)):selected.append((src+1,a,b))
  previous=current;hist.append(len(basis));print('RANK',t,q,order,len(basis),flush=True)
 return dict(point=[t,q],ranks=hist,selected=selected)

if __name__=='__main__':
 results=[]
 for t,q in [(2,3),(5,7),(3,4),(4,7),(6,9),(8,11),(9,13),(11,17),(12,19),(13,23),(15,29),(16,31),(17,37)]+[(v,v+1) for v in range(18,100)]:
  try:result=run(t,q)
  except ValueError as e:
   print('SKIP',t,q,str(e),flush=True);continue
  results.append(result)
  Path(__file__).with_name('goal2_combined_exact.json').write_text(json.dumps(dict(prime=P,pi=101,N=A.N,D=14,results=results),indent=2)+'\n')
  if len(results)==2:break
 if len(results)<2:raise RuntimeError('Need two split generic points')
