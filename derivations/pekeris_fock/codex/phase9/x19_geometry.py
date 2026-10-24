"""Serial exact eliminations for the line-charge/Feynman radical surface."""
import sys
sys.dont_write_bytecode=True
import sympy as S
from pathlib import Path
import json,time,resource
out=Path(__file__).resolve().parent
q,z,t,s,X,Y,u,v=S.symbols('q z t s X Y u v');D=1+q**4;P=S.Matrix([q**4-1,2*q*q]);x=S.Matrix([u,v]);z2=1-u*u-v*v
checks={};facts={};start=time.time()
zero=lambda p:S.factor(p)==0
checks['unit_circle']=zero(P.dot(P)-D*D)
# Direct radicands from all four physical maps (reflections handled explicitly).
for name,pp,b,Y2 in [('panel1',P,u,(1-v)/2),('panel2',S.Matrix([-P[0],P[1]]),-u,(1-v)/2),('panel3',S.Matrix([-P[1],-P[0]]),-v,(1+u)/2),('panel4',S.Matrix([P[1],-P[0]]),-v,(1-u)/2)]:
 raw=(1-q*q)*((1+b)-(1-b)*q*q)+4*Y2*q*q
 checks[name+'_line_charge']=zero(raw-(D-pp.dot(x)))
for name,a in [('v0_e2',S.Matrix([0,1])),('v1_e1',S.Matrix([1,0])),('v1_me1',S.Matrix([-1,0]))]:
 ax=a.dot(x);px=P.dot(x);ap=a.dot(P);cross=a[0]*P[1]-a[1]*P[0]
 A=1-ax*ax;B=ap-ax*px;C=D*D-px*px;N=D-ap
 R=t*D*a+(1-t)*P;Q=S.expand(R.dot(R)-R.dot(x)**2)
 checks[name+'_Q_norm']=zero(Q-(D*D-2*t*(1-t)*D*N-R.dot(x)**2))
 checks[name+'_Q_discriminant']=zero(S.discriminant(Q,t)+4*D*D*cross*cross*z2)
 checks[name+'_gram_elimination']=zero(B*B-A*C+cross*cross*z2)
 checks[name+'_quartic_discriminant']=zero(S.discriminant(t*(1-t)*Q,t)+4*D**6*cross*cross*z2*A*A*C*C)
 # Keep the outer 1/sqrt(D) density: its inclusion changes the twist.
 Qs=S.cancel(Q.subs(t,s/(1+s))*(1+s)**2)
 checks[name+'_cubic_change']=zero(S.cancel((D*s*Qs).subs(s,X/D))-X*(A*X*X+2*B*X+C))
 cubic=X**3+2*B*X*X+A*C*X
 checks[name+'_Weierstrass_discriminant']=zero(16*S.discriminant(cubic,X)+64*A*A*C*C*z2*cross*cross)
 j=64*(4*B*B-3*A*C)**3/(A*A*C*C*(B*B-A*C))
 samples=[S.factor(j.subs({u:S.Rational(1,3),v:S.Rational(1,4),q:qq})) for qq in [S.Rational(1,2),S.Rational(2,3)]]
 checks[name+'_nonisotrivial']=samples[0]!=samples[1]
 # Generic minimality/smoothness at a fresh rational target; exact fibre factor checks.
 sub={u:S.Rational(1,3),v:S.Rational(1,4)}
 cf=S.Poly(C.subs(sub),q);cr=S.Poly(cross,q)
 checks[name+'_eight_simple_C_roots']=cf.degree()==8 and S.gcd(cf,cf.diff()).degree()==0
 checks[name+'_no_C_cross_collision']=S.gcd(cf,cr).degree()==0
 facts[name]=dict(A=str(A),B=str(S.factor(B)),C=str(S.factor(C)),N=str(S.factor(N)),cross=str(cross),Q=str(S.factor(Q)),Weierstrass=str(cubic),j_samples=list(map(str,samples)),q_degrees_a2_a4=[S.degree(2*B,q),S.degree(A*C,q)],z_degrees_a2_a4=[S.degree((2*B).subs(q*q,z),z),S.degree((A*C).subs(q*q,z),z)])
# The nonconstant density cannot be pulled through the inverse. Give an exact rational two-parameter representation.
b,r,w=S.symbols('b r w');H=(1+b)/(1-b)
checks['density_double_parameter_kernel']=zero(H/((H+r*r)*(H+w*w))-(1-b*b)/((1+b+(1-b)*r*r)*(1+b+(1-b)*w*w)))
# Endpoint/presentation-denominator control is genuinely a multiple of every x16 first-order row clearing.
r1,r2,c=S.symbols('r1 r2 r12');Q2=r1*r1+r2*r2;delta=4*r1*r1*r2*r2-(Q2-c*c)**2
common=r1**3*r2**3*c**3*Q2*(2*Q2-c*c)**2*delta**2
nclear=0
for line in (out.parent/'phase8/x16_clearing.txt').read_text().splitlines():
 den=S.sympify(line.split('|')[-1]);num,de=S.fraction(S.cancel(common/den));assert not de.free_symbols;(nclear:=nclear+1)
checks['all_92_row_clearings_divide_control']=nclear==92
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
result=dict(checks=checks,surfaces=facts,old_target_curve=str(D-P.dot(x)),normalized_measure_surface='Y0^2=D*t*(1-t)*Q(q,t)',cubic='Y^2=X^3+2*B(q)*X^2+A*C(q)*X',geometry_qualification='Generic auxiliary K3 elliptic surface over q; necessity of its H2 periods in the final physical answer is NOT proved. Full source-dependent weight-two density needs additional parameters.',elapsed_seconds=time.time()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(out/'x19_geometry.json').write_text(json.dumps(result,indent=2,default=str)+'\n');print(json.dumps(result,indent=2,default=str))
