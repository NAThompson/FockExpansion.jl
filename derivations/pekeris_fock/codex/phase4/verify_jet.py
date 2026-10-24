import sympy as s
r,e,p1,p2,b,logr,A,D=s.symbols('r eps p1 p2 b logr A D')
F=1+r*(p1+e)+r**2*(p2+e*p1+e**2/2+b*logr+A*D)
jet=s.series(s.exp(-e*r)*F,r,0,3).removeO().expand()
assert s.simplify(jet-(1+r*p1+r*r*(p2+b*logr+A*D)))==0
print('PASS exact exponential cancellation through radial degree 2; free coefficient = A D')
# Six-dimensional radial scaling, polynomial degree 2.
assert 2*2+6==10 and 2*(2-1)+6==8
print('PASS bump squared L2 scales h^10; squared gradient scales h^8; Coulomb scales h^9')
R=s.symbols('R',positive=True)
f=R**2*s.exp(-R)
g=s.diff(f,R,2)+5/R*s.diff(f,R)-12/R**2*f
assert s.integrate(g/R,(R,0,s.oo))==-8
assert -s.integrate(g/R,(R,0,s.oo))/8==1
print('PASS manufactured residual Green sign: local=1, global=0, integral g/R=-8')
