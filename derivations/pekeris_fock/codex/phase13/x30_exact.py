"""Exact polarization, coefficient linearity, and endpoint counterexample.
No imports with file-writing side effects. stdout only.
"""
import sympy as s
from itertools import product
Z=s.symbols('Z')
# An ordered word records successive Coulomb choices at steps 1,2,3.
words={''.join(w):s.prod(-Z if a=='1' else 1 for a in w) for w in product('01',repeat=3)}
assert {w:s.expand(c).coeff(Z,1) for w,c in words.items() if s.expand(c).coeff(Z,1)}=={'001':-1,'010':-1,'100':-1}
c0,c2,q,a,d,h,P,y=s.symbols('c0 c2 q a d h P y',nonzero=True)
integrand=2/s.pi*(c0*(a-d)*(a+d)+c2*((h*d/q)**2-(q*a)**2))/s.sqrt(P+4*y*y*q*q)
assert s.simplify(integrand-c0*integrand.subs({c0:1,c2:0})-c2*integrand.subs({c0:0,c2:1}))==0
v=s.cos(s.pi/4)/s.sin(s.pi/12)
assert s.simplify(v-(1+s.sqrt(3)))==0
a0,a1,a2,a3=s.symbols('a0 a1 a2 a3')
f=a0+a1*Z+a2*Z**2+a3*Z**3
assert s.expand((f.subs(Z,1)-f.subs(Z,-1))/2)==a1+a3
assert s.expand((8*(f.subs(Z,1)-f.subs(Z,-1))-(f.subs(Z,2)-f.subs(Z,-2)))/12)==a1
print('PASS exact cubic odd-extraction contamination and corrected extraction')
print('PASS exact ordered-word polarization: coefficient Z = -(A+B+C)')
print('PASS moment integrand linearity in parameter coefficients c0,c2')
print('PASS raw endpoint angular integral = 1+sqrt(3), off diagonal')
