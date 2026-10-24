"""Independent TASK12 exact specialization driver. No claude6 imports.
Derivative rows include differentiated coefficient rows, not matrix powers.
"""
import sys,types,json,time,math,itertools,hashlib,os
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
PRIME=1004535809
import sympy as sp

def load(pi=103,si=1,sr=1,signs=(1,1,1),N=9):
 a=types.ModuleType('goal2_classical_exact');sys.modules[a.__name__]=a
 src=(ROOT/'codex/phase5/goal2_classical_exact.py').read_text()
 src=src.replace('P=998244353; N=8',f'P={PRIME}; N=8').replace('pi=101;Z=2',f'pi={pi};Z=2')
 exec(compile(src,'codex_classical_embedded','exec'),a.__dict__)
 a.I=int(sp.sqrt_mod(-1,PRIME))*si%PRIME;a.RT=int(sp.sqrt_mod(2,PRIME))*sr%PRIME
 b=types.ModuleType('goal2x_connection');b.__file__=str(ROOT/'codex/phase5/goal2_combined_exact.py')
 src=(ROOT/'codex/phase5/goal2_combined_exact.py').read_text().replace('A.N=15',f'A.N={N}').replace('J(101)',f'J({pi})')
 exec(compile(src,'codex_combined_embedded','exec'),b.__dict__)
 # C only accelerates exact modular convolution; no parallelism.
 import ctypes as ct
 lib=ct.CDLL(str(OUT/'goal2x_fast.dylib')); U=ct.c_uint32; Int=ct.c_int
 offsets=[0];aa=[];bb=[]
 for pairs in a.PAIRS:
  for x,y in pairs:aa.append(x);bb.append(y)
  offsets.append(len(aa))
 off=(Int*len(offsets))(*offsets);av=(Int*len(aa))(*aa);bv=(Int*len(bb))(*bb)
 oldmul=b.J.__mul__;oldinv=b.J.inv;pure=bool(os.environ.get("GOAL2X_PURE"))
 def fastmul(x,y):
  if pure:return oldmul(x,y)
  if not isinstance(y,b.J):return oldmul(x,y)
  n=(b.WORK+1)*(b.WORK+2)//2
  for obj in (x,y):
   if not hasattr(obj,'kind'):obj.kind=2 if any(obj.v[1:]) else (1 if obj.v[0] else 0)
  if x.kind==0 or y.kind==0:return b.J([0]*n)
  if x.kind==1:return b.J([x.v[0]*v%PRIME for v in y.v[:n]])
  if y.kind==1:return b.J([y.v[0]*v%PRIME for v in x.v[:n]])
  xx=(U*len(x.v))(*x.v);yy=(U*len(y.v))(*y.v);zz=(U*n)()
  lib.convolution(n,xx,yy,zz,off,av,bv)
  return b.J(list(zz))
 def fastinv(x):
  if pure:return oldinv(x)
  if not x.v[0]:raise ZeroDivisionError()
  n=len(x.v);xx=(U*n)(*x.v);zz=(U*n)()
  lib.inverse(n,xx,zz,off,av,bv)
  return b.J(list(zz))
 b.J.__mul__=fastmul;b.J.__rmul__=fastmul;b.J.inv=fastinv
 orig=b.sqrtj;calls=[0]
 def root(x):
  k=calls[0];calls[0]+=1
  return signs[k]*orig(x)
 b.sqrtj=root
 return a,b

BASIS=[tuple(v) for v in json.loads((ROOT/'codex/phase5/goal2_candidate_generators.json').read_text())['derivatives']]
# All missing derivatives needed to close the four order-eight rows, and order-seven rows.
TARGETS=[(3,a,8-a) for a in range(4,9)]+[(3,a,9-a) for a in range(5)]

def solve(rows,targets,p):
 basis=[]
 for k,src in enumerate(rows):
  row=src[:];comb=[int(j==k) for j in range(len(rows))]
  for pivot,v,c in basis:
   f=row[pivot]
   if f:
    row=[(x-f*y)%p for x,y in zip(row,v)];comb=[(x-f*y)%p for x,y in zip(comb,c)]
  pivot=next((j for j,v in enumerate(row) if v),None)
  if pivot is None:raise ValueError('rank deficient basis')
  inv=pow(row[pivot],-1,p);basis.append((pivot,[x*inv%p for x in row],[x*inv%p for x in comb]))
 answer=[]
 for src in targets:
  row=src[:];comb=[0]*len(rows)
  for pivot,v,c in basis:
   f=row[pivot]
   if f:
    row=[(x-f*y)%p for x,y in zip(row,v)];comb=[(x+f*y)%p for x,y in zip(comb,c)]
  assert not any(row),'target outside span'
  # Recompute all 240 coordinates independently of the elimination residual.
  assert all(sum(c*r[j] for c,r in zip(comb,rows))%p==src[j] for j in range(len(src)))
  answer.append(comb)
 return answer

def relations(t,q,pi=103,si=1,sr=1,signs=(1,1,1),targets=None):
 targets=TARGETS if targets is None else targets
 a,b=load(pi,si,sr,signs);L,C,T=b.build(t,q);prev={(0,0):L};saved={}
 want=set(BASIS+targets);chi_order=max([2]+[aa+bb for ss,aa,bb in targets if ss==2])
 for order in range(10):
  b.WORK=a.N-order;current={}
  for x in range(order+1):
   y=order-x
   if not order:rows=L
   else:
    axis=0 if x else 1;old=prev[x-1,y] if x else prev[x,y-1];rows=[]
    for src,row in enumerate(old):
     if (src==0 and order>0) or (src==1 and order>chi_order):rows.append({});continue
     out={i:v.diff(axis) for i,v in row.items()}
     for i,v in row.items():
      for j,w in C[axis].get(i,[]):out[j]=out.get(j,b.J())+v*w
     rows.append(out)
   current[x,y]=rows
   for src,row in enumerate(rows):
    key=(src+1,x,y)
    if key in want:saved[key]=b.project(row,T)
  prev=current
 return solve([saved[k] for k in BASIS],[saved[k] for k in targets],PRIME)

if __name__=='__main__':
 start=time.time()
 for t,q in [(7,12),(10,17),(14,23),(19,31),(23,37),(29,43)]:
  try:r=relations(t,q)
  except ValueError as e:print('SKIP',t,q,e,flush=True);continue
  print('PASS',t,q,'seconds',time.time()-start,flush=True)
  (OUT/'goal2x_first.json').write_text(json.dumps(dict(prime=PRIME,point=[t,q],pi=103,basis=BASIS,targets=TARGETS,relations=r),indent=2)+'\n');break
