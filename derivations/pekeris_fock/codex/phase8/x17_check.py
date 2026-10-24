"""Small serial exact panel-3 Hamiltonian cross-check, independent product rule."""
import sys
sys.dont_write_bytecode=True
import ast,json,hashlib,resource,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'codex'))
import codex_derive as D
start=time.time()
u,v,w=D.x; sig=u+v+w; zero=(0,0,0)
r=S.symbols('r1 r2 r12'); loc=dict(zip(('r1','r2','r12'),r))
inv=dict(zip(r,((2*v+w)/4,(2*u+w)/4,(u+v)/2)))
# All retained atoms have radial degree zero, so eps cancels from the
# pullback of their first derivative connection. Keep eps in scalar P.
gdist={}; inputs=[OUT/'x16_connection.txt',ROOT/'claude8/panel3_terms.txt',ROOT/'claude8/panel3_rowfactors.json',ROOT/'codex/codex_derive.py']
for line in inputs[0].read_text().splitlines():
 p,axis,src,tgt,num,den=line.split('|')
 if p=='3' and src in ('K_A0','F0'):
  gdist[src,tgt,r.index(loc[axis])]=S.sympify(num,locals=loc)/S.sympify(den,locals=loc)
g={}
for src,tgt,_ in gdist:
 for j,z in enumerate(D.x):
  g[src,tgt,j]=S.cancel(sum(gdist.get((src,tgt,i),0).subs(inv)*S.diff(inv[r[i]],z) if hasattr(gdist.get((src,tgt,i),0),'subs') else 0 for i in range(3)))
base={k:{d:f.subs({D.E:-D.eps**2,D.Z:2,D.lam:1}) for d,f in o.items()} for k,o in D.derive().items()}
M=S.zeros(3)
for i in range(3):
 for j in range(3):
  d=tuple(int(k==i)+int(k==j) for k in range(3));M[i,j]=base['P'][d]/(1 if i==j else 2)
B=[base['P'][tuple(int(k==i) for k in range(3))] for i in range(3)]
selected=[('sigma2P','K_A0','K_A0'),('sigma2P','F0','K_A0'),('C1','F0','K_A0'),('C1','K_A0','K_A0')]
own={}
for name,tgt,src in selected:
 op=dict(base['S2P' if name=='sigma2P' else name]) if tgt==src else {}
 if name=='sigma2P':
  for j in range(3):
   d=tuple(int(k==j) for k in range(3))
   op[d]=op.get(d,0)+2*sig**2*sum(M[i,j]*g.get((src,tgt,i),0) for i in range(3))
  op[zero]=op.get(zero,0)+sig**2*(sum(B[i]*g.get((src,tgt,i),0) for i in range(3))+sum(M[i,j]*(S.diff(g.get((src,tgt,j),0),D.x[i])+sum(g.get((src,k,j),0)*g.get((k,tgt,i),0) for k in ('K_A0','F0'))) for i in range(3) for j in range(3)))
 else:
  op[zero]=op.get(zero,0)+2*sig*sum(M[i,j]*g.get((src,tgt,j),0) for i in range(3) for j in range(3))
 own[name,tgt,src]=op
 print('Derived',name,tgt,src,flush=True)
# Read Claude only after constructing the independent blocks.
rows={k:S.sympify(f,locals=dict(zip(('u','v','w'),D.x))) for k,f in json.loads(inputs[2].read_text()).items()}
cl={key:{} for key in selected};terms=0
for line in inputs[1].read_text().splitlines():
 ks,ds,ms,num,den=line.split('|');name,tgt,src=ast.literal_eval(ks);key=name,tgt.split(':')[1],src.split(':')[1]
 if key not in cl:continue
 d=tuple(map(int,ds.split()));powers=list(map(int,ms.split()))
 cf=S.Rational(int(num),int(den))*S.prod(z**p for z,p in zip((*D.x,D.eps,S.pi),powers))
 cl[key][d]=cl[key].get(d,0)+cf;terms+=1
result=[]; exports={}
for key,op in own.items():
 cleared={d:S.Poly(S.cancel(rows[key[1]]*f),*D.x,D.eps) for d,f in op.items()}
 count=0; monomials=0
 for d in set(cleared)|set(cl[key]):
  a=cleared.get(d,S.Poly(0,*D.x,D.eps)); b=S.Poly(cl[key].get(d,0),*D.x,D.eps)
  assert a==b,(key,d,S.factor(a.as_expr()-b.as_expr()))
  count+=1;monomials+=len(a.terms()) if not a.is_zero else 0
 exports[str(key)]={str(d):str(S.factor(f.as_expr())) for d,f in cleared.items() if not f.is_zero}
 result.append({'block':key,'derivative_coefficients':count,'monomials':monomials})
 print('PASS exact',key,count,'derivative coefficients;',monomials,'monomials',flush=True)
# macOS ru_maxrss is bytes; Linux is KiB. Resource guard is checked on peak.
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
rssbytes=rss if sys.platform=='darwin' else rss*1024
assert rssbytes<2_000_000_000,rssbytes
summary={'blocks':result,'claude_selected_term_lines':terms,'peak_rss_bytes':rssbytes,'elapsed_seconds':time.time()-start,'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
(OUT/'x17_check.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'x17_operators.json').write_text(json.dumps(exports,indent=2)+'\n')
print('PASS all comparisons; peak RSS bytes',rssbytes,flush=True)
