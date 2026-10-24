"""Exact reduction of the four quartic-root adjunctions to two."""
import sympy as s
r,t,c,A=s.symbols('r t c A');Q=r*r+t*t;P2=2*Q-c*c;D=(Q-c*c)**2-4*r*r*t*t
v=[(Q-c*c+A)/(2*t*t),(Q-c*c+A)/(2*r*r),(t*t-r*r+A)/P2,(r*r-t*t+A)/P2]
def red(f):
 num,den=s.fraction(s.cancel(f));return s.factor(s.rem(s.Poly(num,A),s.Poly(A*A-D,A)).as_expr()),s.factor(den)
assert red(v[1]-t*t/r**2*v[0])[0]==0
assert red(v[2]*v[3]+c*c/P2)[0]==0
print('PASS z2^2=(r2/r1)^2 z1^2')
print('PASS z3^2 z4^2=-c^2/P^2; choose z4=i*c/(P*z3)')
print('Sufficient field K(R,P,A,z1,z3), degree <= 32, not a minimality assertion')
print('Previous combined rational closure bound improves to 216*32 =',216*32)
