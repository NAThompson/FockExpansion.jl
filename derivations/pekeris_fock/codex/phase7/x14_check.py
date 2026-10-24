"""Independent small exact classical comparison; read-only Claude inputs."""
import sys, os, json, ast, hashlib
sys.dont_write_bytecode=True
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'codex'))
import codex_derive as D
u,v,w=D.x; sig=u+v+w; xyz=D.x
raw=json.loads((ROOT/'claude7/classical_ops.json').read_text())
loc={'u':u,'v':v,'w':w,'epsilon':D.eps,'PI':S.pi}
vals={D.eps:1,D.E:-1,D.Z:2,D.lam:1}
parse=lambda x:S.sympify(x,locals=loc).subs(vals)
ops={ast.literal_eval(k):{ast.literal_eval(d):parse(f) for d,f in o.items()} for k,o in raw['ops'].items()}
rows={k:parse(f) for k,f in raw['D'].items()}
base={k:{d:f.subs(vals) for d,f in o.items()} for k,o in D.derive().items()}
M=S.zeros(3); B=[]
for i in range(3):
 d=tuple(int(j==i) for j in range(3));B.append(base['P'][d])
 for j in range(3):
  d=tuple(int(k==i)+int(k==j) for k in range(3)); M[i,j]=base['P'][d]/(1 if i==j else 2)
x=(2*v+w)/4; y=(2*u+w)/4; Q=x*x+y*y
# da = (-2y dx +2x dy)/Q; dL(alpha/2)=(a+pi/2)*x/(2y) da.
ga=[S.cancel((-2*y*S.diff(x,z)+2*x*S.diff(y,z))/Q) for z in xyz]
gamma={}
for i in range(3):
 gamma['g13',i]={'g17':ga[i]}
 gamma['g15',i]={'g13':2*ga[i]}
 gamma['g9',i]={'g13':S.cancel(x*ga[i]/(4*y)),'g17':S.cancel(S.pi*x*ga[i]/(8*y))}
 gamma['g17',i]={}
names=['g17','g13','g15','g9']; own={}
def add(n,t,b,d,c):
 key=n,t,b; own.setdefault(key,{})[d]=own.setdefault(key,{}).get(d,0)+c
zero=(0,0,0)
for src in names:
 for n,key in [('sigma2P','S2P'),('C1','C1'),('C2','C2')]:
  for d,f in base[key].items():add(n,src,src,d,f)
 for j in range(3):
  der=tuple(int(k==j) for k in range(3))
  for i in range(3):
   for tgt,c in gamma[src,i].items():add('sigma2P',tgt,src,der,2*sig**2*M[i,j]*c)
 for i in range(3):
  for tgt,c in gamma[src,i].items():
   add('sigma2P',tgt,src,zero,sig**2*B[i]*c)
   add('C1',tgt,src,zero,2*sig*sum(M[k,i] for k in range(3))*c)
  for j in range(3):
   for tgt,c in gamma[src,j].items():
    add('sigma2P',tgt,src,zero,sig**2*M[i,j]*S.diff(c,xyz[i]))
    for t2,c2 in gamma[tgt,i].items():add('sigma2P',t2,src,zero,sig**2*M[i,j]*c*c2)
count=0;blocks=0
for key in sorted(set(own)|{k for k in ops if k[2] in names}):
 n,t,b=key
 if b not in names:continue
 o=own.get(key,{});cl=ops.get(key,{})
 for d in set(o)|set(cl):
  diff=S.cancel(rows[t]*o.get(d,0)-cl.get(d,0))
  assert diff==0,(key,d,diff)
  count+=1
 blocks+=1
print('PASS exact cleared operators:',blocks,'blocks,',count,'derivative coefficients',flush=True)
# Saved independent eps=1 recurrence operators provide a second normalization check.
a=json.loads((OUT/'angle_recurrence.json').read_text()); row=u*v*w*Q
for target,source,n,factor in [('g17','g17','diag',1),('g17','g13','lower1',1),('g13','g15','lower1',2),('g17','g15','lower2',2)]:
 old={ast.literal_eval(k):parse(v) for k,v in a['operators'][n].items()}
 for d in set(old)|set(own['sigma2P',target,source]):
  assert S.cancel(row*own['sigma2P',target,source].get(d,0)-sig**2*factor*old.get(d,0))==0
print('PASS four independent angle_recurrence operator normalizations',flush=True)
# Stencil shards may be absent or partial; compare actual available coefficients directly.
checks=[]; shardfiles=sorted((ROOT/'claude7').glob('classical_stencils_*.json'))
for path in shardfiles:
 for keystr,entry in json.loads(path.read_text()).items():
  key=ast.literal_eval(keystr)
  if key[2] not in names:continue
  # Small low-index source projection is independent of any stencil expansion.
  op=ops[key]
  for source in [(0,0,0),(1,0,1)]:
   f=S.prod(S.laguerre(j,z) for j,z in zip(source,xyz))
   expected=D.laguerre_coeffs(D.apply(op,f));actual={}
   for ds,cf in entry['rec'].items():
    shift=ast.literal_eval(ds);target=tuple(a-b for a,b in zip(source,shift))
    if min(target)<0:continue
    expr=S.sympify(cf,locals={**loc,'l':D.l,'m':D.m,'n':D.n}).subs(vals).subs(dict(zip(D.indices,target)))
    expr=S.expand(expr)
    if expr!=0:actual[target]=expr
   assert actual==expected,(path.name,key,source)
   checks.append([path.name,key,source,len(expected)])
print('Stencil inventory:',len(shardfiles),'shards;',len(checks),'small action checks',flush=True)
result={'operator_blocks':blocks,'operator_coefficients':count,'angle_normalizations':4,'stencil_files':[p.name for p in shardfiles],'stencil_checks':checks,'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'claude7/classical_ops.json',OUT/'angle_recurrence.json']+shardfiles}}
(OUT/'x14_check.json').write_text(json.dumps(result,indent=2)+'\n')
