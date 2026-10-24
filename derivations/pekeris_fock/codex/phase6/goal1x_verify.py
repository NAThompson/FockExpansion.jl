"""Small serial checks for TASK_codex_11; no external writes/imports."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
import mpmath as m
import sympy as S
m.mp.dps=60
pi=m.pi
cl=lambda z:m.im(m.polylog(2,m.exp(1j*z)))
def parts(alpha,beta):
 a,b=alpha-beta,alpha+beta
 A=m.log(m.tan(b/2)*m.tan(a/2))
 B=m.log(m.tan(b/2)/m.tan(a/2))
 delta=cl(a)+cl(pi-a)-cl(pi-b)-cl(b)
 J=m.log((1+1/m.cos(a/2))*m.tan(pi/4-b/4))
 return beta*A+alpha*B-delta+pi*J,J

def fa(z):return -z/m.sin(z)+pi*m.sin(z/2)/(2*m.cos(z/2)*(1+m.cos(z/2)))
def fb(z):return z/m.sin(z)-pi/(2*m.cos(z/2))
maxerr=m.mpf(0)
for alpha,beta in [(m.mpf('.8'),m.mpf('.2')),(m.mpf('1.4'),m.mpf('-.9')),(m.mpf('.04'),m.mpf('.03'))]:
 a,b=alpha-beta,alpha+beta
 sn,cs=m.sin(alpha),m.cos(alpha); v,g=m.sin(beta),m.cos(beta)
 Jold=m.log((cs+m.sqrt((1+sn)*(1-v)))**2/((1+sn)*(g+cs)))
 U,J=parts(alpha,beta)
 errors=[abs(J-Jold),abs(m.diff(lambda aa:parts(aa,beta)[0],alpha)-fa(a)-fb(b)),abs(m.diff(lambda bb:parts(alpha,bb)[0],beta)+fa(a)-fb(b))]
 maxerr=max(maxerr,*errors)
assert maxerr<m.mpf('1e-54')
print('PASS U angle J identity and both differential components; max error',m.nstr(maxerr,5))
for beta in [m.mpf('-.9'),m.mpf(0),m.mpf('.9')]:
 assert abs(parts(pi/2,beta)[0])<m.mpf('1e-55')
print('PASS U zero exchange trace, three beta values')
z=S.symbols('z',positive=True)
f=-z/S.sin(z)+S.pi*S.sin(z/2)/(2*S.cos(z/2)*(1+S.cos(z/2)))
g=z/S.sin(z)-S.pi/(2*S.cos(z/2))
limits=[S.limit(h,z,e,dir='+' if e==0 else '-') for h in (f,g) for e in (0,S.pi)]
assert limits==[-1,1-S.pi/2,1-S.pi/2,-1]
print('PASS exact removable endpoint limits [Fa(0),Fa(pi),Fb(0),Fb(pi)] =',limits)
# Exact 3x3 bordered/normalized-tail identity, independent of a helium quadrature.
K=S.Matrix([[5,1,2],[1,4,1],[2,1,6]]); inv=K.inv(); u=inv[:,0]/inv[0,0]
assert u[0]==1 and (K*u)[1:]==[0,0] and (u.T*K*u)[0]==1/inv[0,0]
print('PASS exact normalized-tail / reciprocal inverse / defect identity')
# Algebraic minimal sequence can coexist with a geometric dominant sequence.
n=S.symbols('n',positive=True,integer=True)
h=lambda k:1/(k+1)
d=lambda k:2**k
A,B=S.symbols('A B')
sol=S.solve([h(n+1)+A*h(n)+B*h(n-1),d(n+1)+A*d(n)+B*d(n-1)],(A,B))
assert S.simplify(sol[A]+3)==0 or S.limit(sol[A],n,S.oo)==-3
assert S.limit(sol[B],n,S.oo)==2
print('PASS algebraic minimal / geometric dominant counterexample:',sol,'; roots 1,2')
import json
root=Path(__file__).resolve().parents[2]
selected=set(map(tuple,json.loads((root/'codex/phase5/goal2_combined_exact.json').read_text())['results'][0]['selected']))
targets=[(src,a,order-a) for order in (8,9) for a in range(order+1) for src in (1,2,3) if (src,a,order-a) not in selected]
assert targets[:6]==[(1,0,8),(2,0,8),(1,1,7),(2,1,7),(1,2,6),(2,2,6)]
assert next(k for k in targets if k[0]==3)==(3,4,4)
print('PASS target-selection audit: first six',targets[:6],'; first new f target (3,4,4)')
print('ALL CHECKS PASS (no shell solve, no large matrices)')
