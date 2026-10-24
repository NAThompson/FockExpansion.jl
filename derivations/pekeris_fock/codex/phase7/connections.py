"""Small serial symbolic structured connections; never writes imported inputs."""
import os
os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import sympy as S
import json, hashlib
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
b,y,h,z,u=S.symbols('b y h z u'); I=S.I
jets=[(0,0),(1,0),(2,0),(0,1),(0,2)]
names=[]
for fam in ['J','K']+['J'+str(k) for k in range(8)]+['FP','E','F0']:
 for a,c in jets:
  if fam=='K' and (a,c) in [(2,0),(0,2)]:continue
  if fam=='FP' and (a,c)!=(0,1):continue
  if fam=='F0' and (a,c)!=(0,0):continue
  names.append(f'{fam}_{a}{c}')
assert len(names)==55
class Lin:
 def __init__(self,d=None):self.d=d or {}
 def __add__(self,o):
  if o==0:return self
  assert isinstance(o,Lin),o
  d=dict(self.d)
  for k,v in o.d.items():d[k]=d.get(k,0)+v
  return Lin(d)
 __radd__=__add__
 def __neg__(self):return Lin({k:-v for k,v in self.d.items()})
 def __sub__(self,o):return self+-o
 def __mul__(self,o):return Lin({k:v*o for k,v in self.d.items()})
 __rmul__=__mul__
 def __truediv__(self,o):return self*(1/o)
 def __eq__(self,o):return o==0 and not self.d if not isinstance(o,Lin) else self.d==o.d
V=[Lin({i:S.Integer(1)}) for i in range(55)];lookup=dict(zip(names,V))
banned={}
def state(fam,jet):
 name=f'{fam}_{jet[0]}{jet[1]}'
 if fam=='F0' and jet!=(0,0):return Lin()
 if name not in lookup:
  banned.setdefault(name,Lin({name:S.Integer(1)}));return banned[name]
 return lookup[name]
def geometry(bb,yy):return h,[z,-z,h/z,-h/z,I,-I,I*h,-I*h],u
def const(c):return [S.sympify(c),S.Integer(0),S.Integer(0)]
MODE="by"
def data(bb,yy,axis,family):
 _,zs,_=geometry(b,y);hb=h/(1-b*b)
 vz=-(1-z**4 if axis==0 else 8*y*z*z)/(4*(1-b)*z**3+2*(4*y*y-2)*z)
 hv=hb if axis==0 else 0
 vs=[vz,-vz,hv/z-h*vz/z**2,-hv/z+h*vz/z**2,0,0,I*hv,-I*hv]
 if MODE=='hz':
  hv=S.Integer(1) if axis==0 else S.Integer(0)
  vz=S.Integer(0) if axis==0 else S.Integer(1)
  vs=[vz,-vz,hv/z-h*vz/z**2,-hv/z+h*vz/z**2,0,0,I*hv,-I*hv]
 du=hv if SIDE<0 else 0
 aa=[const(-S.Rational(1,2)) for _ in range(4)]
 tj=[0,1,0] if family=='A' else [0,0,0]
 sj=[0,1,0] if family=='D' else [0,0,0]
 aa += [[v/(2*I) for v in tj],[-v/(2*I) for v in tj],[v/(2*I) for v in sj],[-v/(2*I) for v in sj]]
 l0=-1/(2*(1+b)) if axis==0 else 0
 if MODE=='hz':l0=-1/(h*(1+h*h)) if axis==0 else 0
 kap=const(l0)
 for ai,zi,vi in zip(aa,zs,vs):kap=[v-w*vi/zi for v,w in zip(kap,ai)]
 return h,zs,u,vs,du,aa,l0,kap
src=(ROOT/'codex/phase3/panel_connection.py').read_text()
helpers=src[src.index('def add'):src.index('def geometry')]
body=src[src.index('def rhs'):src.index('def moments_from_state')]
ns=dict(S=S,m=SimpleNamespace(mpf=S.Rational),geometry=geometry,data=data)
exec(helpers+body,ns)
# Use exact constants throughout: source's explicit floating half is mpf('.5').
def entries(side,axis):
 global SIDE
 SIDE=side;result={}
 for family in ['A','D']:
  ja=[(0,0),(1,0),(2,0)] if family=='A' else [(0,0),(0,1),(0,2)]
  fams=['J','K']+['J'+str(k) for k in range(8)]+['FP','E','F0']
  X=[[state(f,j) for j in ja] for f in fams]
  Y=ns['rhs'](b,y,axis,family,X)
  for f,row in zip(fams,Y):
   for jet,expr in zip(ja,row):
    name=f'{f}_{jet[0]}{jet[1]}'
    if name not in lookup or (family=='D' and jet==(0,0)):continue
    ri=names.index(name)
    for ci,value in expr.d.items():
     c=S.cancel(S.together(value))
     if c==0:continue
     assert isinstance(ci,int),(name,ci)
     result[ri,ci]=c
 return result

def encode(M,variables):
 data=[];mx=[0,0]
 for (i,j),c in sorted(M.items()):
  n,d=S.fraction(c)
  deg=[int(S.Poly(a,*variables).total_degree()) for a in (n,d)]
  mx=[max(a,k) for a,k in zip(mx,deg)]
  data.append([i,j,str(c),deg])
 return dict(nnz=len(data),max_degrees=mx,entries=data)
def classical():
 t,q=S.symbols('t q');pi=S.pi
 M=[{},{}];args=[]
 for sign in (1,-1):
  ex=t*q**sign
  for mm in (-1,0,1,2):args.append((ex*I**(mm%4) if mm!=2 else ex,1,sign,0 if mm==2 else mm*pi,mm==2,4))
 args += [(q*q,0,1,0,False,2),(I*t*t,1,0,pi/2,False,2),(I/(t*t),-1,0,pi/2,False,2)]
 for axis in (0,1):
  fac=2/(I*(t if axis==0 else q))
  for k,(ex,aa,pp,cc,tan,den) in enumerate(args):
   trig=(ex-1)/(I*(ex+1)) if tan else I*(ex+1)/(ex-1)
   f=fac*trig*(aa if axis==0 else pp)/den
   for col,val in [(13,aa),(14,pp),(17,cc)]:
    if val and f:M[axis][k,col]=S.cancel(f*val/den)
  for k,sgn in [(11,1),(12,-1)]:
   ex=t*q**sgn;M[axis][k,17]=S.cancel(-fac*(ex-1)/(I*(ex+1))*(1 if axis==0 else sgn)/4)
  M[axis][13 if axis==0 else 14,17]=fac
  M[axis][15 if axis==0 else 16,13 if axis==0 else 14]=2*fac
 return dict(atoms=['L(w+ -pi/4)','L(w+)','L(w+ +pi/4)','T(w+)','L(w- -pi/4)','L(w-)','L(w- +pi/4)','T(w-)','L(p/2)','L((a+pi/2)/2)','L((pi/2-a)/2)','lc(w+)','lc(w-)','a','p','a^2','p^2','1'],dt=encode(M[0],[t,q]),dq=encode(M[1],[t,q]))
def run():
 global MODE
 result=dict(convention='dX=M X; Taylor jets divided by factorial; formal degree in b,y,h,z,u, subject to stated algebraic relations',atoms=names,panels={})
 result['classical']=classical()
 for side in [1,-1]:
  for axis in [0,1]:
   key=f'{"plus" if side>0 else "minus"}_{"b" if axis==0 else "y"}'
   result['panels'][key]=encode(entries(side,axis),[b,y,h,z,u])
   print(key,result['panels'][key]['nnz'],result['panels'][key]['max_degrees'],flush=True)
 MODE='hz';result['cover']={}
 for side in [1,-1]:
  for axis in [0,1]:
   key=f'{"plus" if side>0 else "minus"}_{"h" if axis==0 else "z"}'
   M={k:S.cancel(v.subs(u,1 if side>0 else h)) for k,v in entries(side,axis).items()}
   result['cover'][key]=encode({k:v for k,v in M.items() if v!=0},[h,z])
   print(key,result['cover'][key]['nnz'],result['cover'][key]['max_degrees'],flush=True)
 (OUT/'connections.json').write_text(json.dumps(result,indent=2)+'\n')
 print('classical',[(ax,result['classical'][ax]['nnz'],result['classical'][ax]['max_degrees']) for ax in ['dt','dq']],flush=True)
if __name__=='__main__':run()
