"""Independent Hamiltonian subblock {1,a,a^2}; small fixed-scale exact stencils."""
from pathlib import Path
import sys,json
import sympy as S
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import codex_derive as D
OUT=Path(__file__).resolve().parent
u,v,w=D.x;x,y,c=D.a,D.b,D.c
Q=x*x+y*y;mult=-32*x*y*c
# Direct distance Hamiltonian product rule, eps=1:
# (H-E)(a^k F)=a^k(H-E)F+k*a^(k-1)B F-2k(k-1)/Q*a^(k-2)F.
B=[2*y/Q,-2*x/Q,-(x*x-y*y)*(c*c-Q)/(Q*x*y*c)]
B0=-2*(x*x-y*y)/(x*y*Q)+2*(x-y)/Q
J=S.Matrix([[-1,1,1],[1,-1,1],[2,2,-2]])
inv={x:(2*v+w)/4,y:(2*u+w)/4,c:(u+v)/2}
base=D.derive()['P'];vals={D.eps:1,D.E:-1,D.Z:2,D.lam:1}
qp=S.expand(Q.subs(inv));row=u*v*w
ops={'diag':{d:S.expand(row*qp*f.subs(vals)) for d,f in base.items()}}
bvec=J*S.Matrix(B)
op={tuple(int(i==j) for i in range(3)):S.cancel((mult*Q*bvec[j]).subs(inv))*row for j in range(3)}
op[0,0,0]=S.cancel((mult*Q*B0).subs(inv))*row
ops['lower1']={d:S.expand(f) for d,f in op.items()}
ops['lower2']={(0,0,0):S.expand(row*(-2*mult).subs(inv))}

def run():
 recs={};stats={}
 for name,op in ops.items():
  print('start',name,flush=True)
  rec=D.recurrence(op,bound=12);recs[name]=rec
  stats[name]=dict(count=len(rec),shell_range=[min(map(sum,rec)),max(map(sum,rec))])
  print(name,stats[name],flush=True)
  for src in [(0,0,0),(1,0,2),(2,1,1)]:
   f=S.prod(S.laguerre(n,z) for n,z in zip(src,D.x))
   expected=D.laguerre_coeffs(D.apply(op,f));actual={}
   for shift,coef in rec.items():
    target=tuple(a-b for a,b in zip(src,shift))
    if min(target)<0:continue
    value=S.expand(coef.subs(dict(zip(D.indices,target))))
    if value!=0:actual[target]=value
   assert actual==expected,(name,src)
 # Check product rule independently with actual angle and ungauged H.
 angle=2*S.atan(y/x)-S.pi/2
 G=S.Matrix([[1,0,(x*x+c*c-y*y)/(2*x*c)],[0,1,(y*y+c*c-x*x)/(2*y*c)],[(x*x+c*c-y*y)/(2*x*c),(y*y+c*c-x*x)/(2*y*c),2]])
 axes=[x,y,c];drift=[2/x,2/y,4/c]
 def lap(f):return sum(G[i,j]*S.diff(f,axes[i],axes[j]) for i in range(3) for j in range(3))+sum(d*S.diff(f,t) for d,t in zip(drift,axes))
 f=1+x*y+c*c
 for k in [1,2]:
  # Conjugated Laplacian includes -2 grad(s).grad and |grad(s)|^2-lap(s).
  def hg(g):return -lap(g)/2+sum((G[i,0]+G[i,1])*S.diff(g,axes[i]) for i in range(3))
  expected=k*angle**(k-1)*(sum(bb*S.diff(f,t) for bb,t in zip(B,axes))+B0*f)
  if k==2:expected-=4*f/Q
  assert S.simplify(S.expand(hg(angle**k*f)-angle**k*hg(f)-expected))==0
 (OUT/'angle_recurrence.json').write_text(json.dumps(dict(parameters={'eps':1,'E':-1,'Z':2,'lambda':1},row_factor='u*v*w*(r1^2+r2^2)',blocks={name:{str(k):str(v) for k,v in rec.items()} for name,rec in recs.items()},operators={name:{str(k):str(v) for k,v in op.items()} for name,op in ops.items()},stats=stats),indent=2)+'\n')
 print('PASS: 9 exact Laguerre actions, derivative-tail certificates, 2 actual-angle product rules',flush=True)
if __name__=='__main__':run()
