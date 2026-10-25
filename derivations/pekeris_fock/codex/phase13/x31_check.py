"""Exact algebra for TASK 31. No external imports with write side effects."""
import sympy as S
import resource
from itertools import combinations
p,q,h,z,c=S.symbols('p q h z c')
vs=(p,q,h); N=1+p*q+p*h+q*h; D2=S.prod(1+x*x for x in vs)
si=lambda x:2*x/(1+x*x)
co=lambda x:(1-x*x)/(1+x*x)
d=lambda f,x:S.diff(f,x)*(1+x*x)/2
V={x:S.expand(d(N,x)-x*N/2) for x in vs}
W={x:S.expand(d(V[x],x)-x*V[x]/2) for x in vs}
def check(label,ex):
 assert S.cancel(ex)==0,(label,S.factor(ex))
 print('PASS',label)
for x,y in combinations(vs,2):
 delta=si(x)-si(y)
 check(f'{x},{y} Fsecond coefficient',(-si(x)*V[x]**2+si(y)*V[y]**2)-delta*(N*N-D2)/4)
 check(f'{x},{y} Fprime coefficient',(-si(x)*W[x]-co(x)*V[x]+si(y)*W[y]+co(y)*V[y])-3*delta*N/4)
 check(f'{x},{y} potential coefficient',-c*delta/16-delta*(-c/4)/4)
F=(4*z**3-3*z)/S.sqrt(z*z-1)
check('T3 scalar ODE',(z*z-1)*S.diff(F,z,2)+3*z*S.diff(F,z)-8*F)
n=S.symbols('n');T=S.Function('T')(z)
# General substitution identity: RHS reduces to the Chebyshev equation.
G=T/S.sqrt(z*z-1)
check('general Chebyshev substitution', ((z*z-1)*S.diff(G,z,2)+3*z*S.diff(G,z)+(1-n*n)*G)*S.sqrt(z*z-1)-((z*z-1)*S.diff(T,z,2)+z*S.diff(T,z)-n*n*T))
e=S.symbols('e');g=(4*(1+e)**3-3*(1+e))/S.sqrt(2+e)
print('characteristic expansion sqrt(e)*F =',S.series(g,e,0,3))
# Algebraic cap in rational half-angle coordinates.
cap=S.Poly(S.expand(N*N-D2),q)
print('cap F =',S.factor(cap.as_expr()))
print('cap degree in each variable:',[S.degree(cap.as_expr(),v) for v in vs])
print('cap discriminant in q =',S.factor(S.discriminant(cap.as_expr(),q)))
check('cap discriminant',S.discriminant(cap.as_expr(),q)-16*p*h*(1+p*p)*(1+h*h))
# Panel Z^1 output coefficients and chart degeneracy H=1.
b,Y,Z=S.symbols('b Y Z');sq=S.sqrt(2)
c0=-sq*Z*(9*Z*b*b-18*Z-8*b*b+8)/576
c2=-sq*Z*(1-b)*(90*Z*Y-45*Z-18*b-32*Y+16)/288
check('Z1 c0',S.expand(c0).coeff(Z,1)-sq*(b*b-1)/72)
check('Z1 c2',S.expand(c2).coeff(Z,1)-sq*(1-b)*(9*b+16*Y-8)/144)
# A generic logarithmic source is integrable in the flux, not an endpoint delta.
x=S.symbols('x',positive=True)
for m in range(5):
 assert S.limit(x*S.log(x)**m,x,0,dir='+')==0
print('PASS all five polylog flux limits x log(x)^m -> 0')
print('ru_maxrss bytes',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
print('ALL EXACT CHECKS PASSED')
