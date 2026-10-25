"""Exact identities and a source-versus-variation counterexample, not a B_A evaluator."""
from pathlib import Path
import sympy as s,json,resource
x,y,u,t=s.symbols('x y u t',real=True)
f=s.Function('f')(x);h=s.Function('h')(x);C=s.symbols('C')
M=lambda z:s.cos(x)*s.diff(z,x,2)-s.sin(x)*s.diff(z,x)+2*s.cos(x)*z
checks={}
checks['lagrange_identity']=s.simplify(s.diff(s.cos(x)*(h*s.diff(f,x)-s.diff(h,x)*f),x)-(h*M(f)-f*M(h)))==0
checks['mode_zero']=s.trigsimp(M(s.sin(x)))==0
for sign in [-1,1]:
 hh=1+sign*s.sqrt(2)*s.cos(x)
 checks['mode_'+str(sign)]=s.trigsimp(M(hh)-sign*s.sqrt(2)*hh)==0
# An analytic order-three test input has N_b phi3=0, but a nonzero first
# log coefficient after multiplication by wv0 and the resonant inverse.
# phi3=sin(alpha). r=wv0*phi3/8. Extract the sin(Y) mode.
a=(x+y)/2;b=(y-x)/2
r=s.sqrt(2)/4*s.sin(a)**2*s.cos(b/2)
# Expand into complex exponentials: integrate each term exactly, avoiding
# large nested symbolic trigonometric antiderivatives.
expr=s.expand((s.sin(x)*s.sin(y)*r).rewrite(s.exp))
lo=s.pi/2;hi=3*s.pi/2
J=0
for term in s.Add.make_args(expr):
 exps=term.atoms(s.exp);ee=sum(e.args[0] for e in exps)
 coeff=term/s.prod(exps)
 kx=s.simplify(s.diff(ee,x)/s.I);ky=s.simplify(s.diff(ee,y)/s.I)
 def int_exp(k):
  return hi-lo if k==0 else (s.exp(s.I*k*hi)-s.exp(s.I*k*lo))/(s.I*k)
 J+=coeff*int_exp(kx)*int_exp(ky)
J=s.simplify(s.expand_complex(J))
B=s.simplify(-2*J/s.pi) # norm sin(Y)=pi/2; h(Xb)=-1
checks['analytic_source_nonzero_flux']=B!=0
# sqrt(u) twist: if p satisfies Lp=0, z=sqrt(u)*p.
A=4*u*u*(u*u-1);BB=4*u*(3*u*u-1);CC=4*u*u+1
q=-1/(2*u)
checks['sqrt_twist_D']=s.simplify(BB+2*A*q-8*u**3)==0
checks['sqrt_twist_0']=s.simplify(CC+BB*q+A*(s.diff(q,u)+q*q))-u*u==0
# Round exactly the same doubled coefficients as fac_A3.sage, without
# claiming that the rounded operator annihilates the actual physical g.
rows=[line.split() for line in Path('claude11/op_A3.txt').read_text().splitlines()]
cs={(int(k),int(j)):int(round(2*s.Rational(v))) for k,j,v in rows}
P=[sum(cs[k,j]*t**j for j in range(13)) for k in range(3)]
facts={'analytic_input_flux_sin_mode':str(B),'flux_100_digits':str(s.N(B,102)),
 'rounded_A3_coefficients':[str(s.factor(p)) for p in P],
 'sqrt_u_twist':'4*u^2*(u^2-1)*D^2 + 8*u^3*D + u^2',
 'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
assert all(checks.values()),checks
out={'checks':checks,'facts':facts}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
