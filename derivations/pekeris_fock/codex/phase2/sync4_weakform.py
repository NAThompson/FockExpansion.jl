#!/usr/bin/env python3
import sys
sys.dont_write_bytecode=True
import sympy as S
import mpmath as mp
s,x,y=S.symbols('s x y',positive=True)
a,b,d=S.symbols('a b d',positive=True)
r=S.Matrix([s*(1+x)/2,s*(1-x)/2,s*y]);J=r.jacobian((s,x,y))
weight=S.factor((a*b*d).subs(dict(zip((a,b,d),r)))*(-J.det()))
assert S.simplify(weight-s**5*(1-x*x)*y/8)==0
# Gram of the six-dimensional gradients of r1,r2,r12.
g=S.Matrix([[1,0,(a*a+d*d-b*b)/(2*a*d)],
            [0,1,(b*b+d*d-a*a)/(2*b*d)],
            [(a*a+d*d-b*b)/(2*a*d),(b*b+d*d-a*a)/(2*b*d),2]])
g=g.subs(dict(zip((a,b,d),r)));w0=weight/s**5
expected=S.Matrix([[w0,0,(1-x)*(x+y*y)/8],
                   [0,w0,(1+x)*(y*y-x)/8],
                   [(1-x)*(x+y*y)/8,(1+x)*(y*y-x)/8,2*w0]])
assert (g*w0-expected).applyfunc(S.simplify)==S.zeros(3)
V=-2/a-2/b+1/d
assert S.simplify((s*V).subs(dict(zip((a,b,d),r)))*w0-(-y/1+(1-x*x)/8))==0
inverse=J.inv()
assert inverse.applyfunc(S.simplify)==S.Matrix([[1,1,0],[(1-x)/s,-(1+x)/s,0],[-y/s,-y/s,1/s]])
print('PASS exact measure, distance-gradient metric/7 kinetic weights, potential and coordinate chain rule')
mp.mp.dps=60
for m in (3,9):
 for p in range(5):
  e=mp.mpf('3.408');h=m+1;z=[mp.digamma(h)-mp.log(e)]+[mp.polygamma(k,h) for k in (1,2,3)]
  bells=[1,z[0],z[0]**2+z[1],z[0]**3+3*z[0]*z[1]+z[2],z[0]**4+6*z[0]**2*z[1]+4*z[0]*z[2]+3*z[1]**2+z[3]]
  direct=mp.diff(lambda t:mp.gamma(h+t)/e**(h+t),0,p)
  assert abs(direct-mp.gamma(h)/e**h*bells[p])<mp.mpf('1e-52')
print('PASS 10 independent Gamma-derivative radial log moments through p=4')
# Manufactured projected f(R)=R^2 log R exp(-R) has b=1,a=0.
gam=mp.euler;constant=(1-gam)+9*gam-8*gam-2
assert abs(constant+1)<mp.mpf('1e-55')
assert abs(-constant/8-mp.mpf(1)/8)<mp.mpf('1e-55')
print('PASS manufactured radial finite part=-1, recovered a=0')
