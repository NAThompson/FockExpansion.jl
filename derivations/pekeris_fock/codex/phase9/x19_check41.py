"""Independent exact distance-connection audit, only x19_* writes."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import sympy as S
import pickle,json,hashlib,time,resource
root=Path(__file__).resolve().parents[2]; out=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'codex/phase2'))
import order2 as O
start=time.time(); x,y,c=S.symbols('r1 r2 r12',positive=True); xyz=(x,y,c); Q=x*x+y*y
raw=pickle.loads((root/'claude9/closed41.pkl').read_bytes())
loc={str(v):v for v in xyz}; loc['pi']=S.pi
T={k:v.subs({s:loc.get(str(s),s) for s in v.free_symbols}) for k,v in raw.items()}
ah,bh=S.symbols('a_h b_h')
z=S.Matrix([T['1'],T['L'],T['alpha'],c*T['V'],T['U'],0])
C=O.connection(xyz)
D=lambda v,i:(v.diff(xyz[i])+C[i]*v).applyfunc(S.cancel)
z1=[D(z,i) for i in range(3)]
z2={(i,j):D(z1[j],i) for i in range(3) for j in range(i,3)}
lap=z2[0,0]+z2[1,1]+2*z2[2,2]+2/x*z1[0]+2/y*z1[1]+4/c*z1[2]+(x*x-y*y+c*c)/(x*c)*z2[0,2]+(y*y-x*x+c*c)/(y*c)*z2[1,2]
B=(S.pi-2)/(3*S.pi)
source=B*((x+y)/(x*y)*(5*c**3/S.Integer(6)-c*Q)+(x+y)*(c-Q/c))
psi42R2=(S.pi-2)*(5*S.pi-14)/(180*S.pi**2)*(Q-2*(Q-(x*x-y*y)**2/Q-(Q-c*c)**2/Q))
print('Connection and covariant derivatives constructed',flush=True)
res=[S.factor(S.cancel(S.together(v+(source+24*psi42R2 if i==0 else 0)))) for i,v in enumerate(lap)]
assert all(v==0 for v in res),res
hom=[S.factor(sum(t*S.diff(v,t) for t in xyz)-4*v) for v in z];assert all(v==0 for v in hom)
print('PASS exact homogeneous Laplacian in six-state Codex connection; all six residuals zero',flush=True)
print('PASS coefficient degree four, no Slog; free a_h,b_h harmonic additions cancel exactly',flush=True)
# Reconstruct from the identified integer relations independently, including the undetermined harmonic polynomial.
import re
R=S.symbols('R',positive=True);rec={}
for line in (root/'claude9/coef5_ident.txt').read_text().splitlines():
 lab,val,rel=map(str.strip,line.split('|'))
 if rel=='NONE':continue
 nums=list(map(int,re.findall(r'-?\d+',rel))); a=nums
 coeff=-(S.Integer(a[1])+a[2]*S.pi+a[3]*S.pi**2+a[4]/S.pi+a[5]/S.pi**2)/a[0]
 atom,mono=lab.split(' ',1); i,j,k=map(int,re.findall(r'\^(\d+)',mono));d=4-i-j-k-(atom=='V')
 rec[atom]=rec.get(atom,0)+coeff*x**i*y**j*c**k*Q**(S.Rational(d,2))
rec['1']+=(5-3/S.pi-14/S.pi**2)/360*Q**2+ah*(4*(x*x-y*y)**2-Q**2)+bh*(S.Rational(3,2)*(Q-c*c)**2-2*x*x*y*y)
assert all(S.cancel(T[k]-rec[k])==0 for k in T)
print('PASS all identified integer relations reconstruct closed41.pkl exactly',flush=True)
# Export the coefficient row with zero harmonic additions, for independent actual-function and purity tests.
lines=[]
for a,v in enumerate(z[:5]):
 for powers,coeff in S.Poly(v.subs({ah:0,bh:0}),x,y,c).terms():
  lines.append(' '.join(map(str,(a+1,*powers,S.N(coeff,78)))))
(out/'x19_closed41_coefficients.txt').write_text('\n'.join(lines)+'\n')
# Baseline constants supplied by Claude's spectral fit, used only as a comparison.
co={l.split('|')[0].strip():S.Float(l.split('|')[1].strip(),78) for l in (root/'claude9/coef5.txt').read_text().splitlines()}
bv=co['1 r1^0 r2^0 r12^4']/S.Rational(3,2); av=(co['1 r1^4 r2^0 r12^0']-2*bv)/16
(out/'x19_claude_harmonics.txt').write_text(str(S.N(av,70))+' '+str(S.N(bv,70))+'\n')
result=dict(exact_PDE_residuals=list(map(str,res)),homogeneity_residuals=list(map(str,hom)),integer_relations_match=True,rank_upper_bound=5,qualification='Membership and upper bound; this check does not assert the single function generates rank exactly five.',elapsed_seconds=time.time()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,inputs={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'claude9/closed41.pkl',root/'claude9/coef5_ident.txt',root/'claude9/coef5.txt',root/'codex/phase2/order2.py']})
(out/'x19_check41.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
