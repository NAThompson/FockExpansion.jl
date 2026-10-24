"""Exact distance homogenization of the polynomial and E/a21 parts."""
import sympy as s
R,S,c,r1,r2,Z,E,a,G,l=s.symbols('R S c r1 r2 Z E a G l',nonzero=True)
p=s.pi;ell=l+4*G/p;x=c/R;sg=S/R
raw0=x*(3-2*x*x)/72
raw1=-x*(5*x*x-6)*ell/72-p*x*(1-x*x)/36+x*(91*x*x-106)/(144*p)+(4*sg**3+4*sg**2*x-32*sg*x*x+4*sg-27*x**3+26*x)/288
raw2=-sg*(1-x*x)*ell/12+53*sg*(1-x*x)/(72*p)+p*(16*sg**3-4*sg*x*x-28*sg-45*x**3+45*x)/576-(503*sg**3-375*sg**2*x+384*sg*x*x-1212*sg-656*x**3+960*x)/864
raw3=-sg*(5*sg**2-3)/36
Q=R*R-c*c
P0=c*(3*R*R-2*c*c)/72
P1=-c*(5*c*c-6*R*R)*ell/72-p*c*Q/36+c*(91*c*c-106*R*R)/(144*p)+(4*S**3+4*S*S*c-32*S*c*c+4*S*R*R-27*c**3+26*c*R*R)/288
P2=-S*Q*ell/12+53*S*Q/(72*p)+p*(16*S**3-4*S*c*c-28*S*R*R-45*c**3+45*c*R*R)/576-(503*S**3-375*S*S*c+384*S*c*c-1212*S*R*R-656*c**3+960*c*R*R)/864
P3=-S*(5*S*S-3*R*R)/36
for i,(raw,phy) in enumerate(zip((raw0,raw1,raw2,raw3),(P0,P1,P2,P3))):
 assert s.cancel(R**3*raw-phy)==0
 print('EXACT polynomial Z power',i,'PASS')
state=E*(Z*sg*(2+r1*r2/R**2)/18-(6-x*x)*x/72)-a*(Z*sg*(1-x*x)/2-(6-5*x*x)*x/12)
ps=E*(Z*S*(2*R*R+r1*r2)/18-c*(6*R*R-c*c)/72)-a*(Z*S*Q/2-c*(6*R*R-5*c*c)/12)
assert s.cancel(R**3*state-ps)==0
print('EXACT E/a21 state polynomial PASS')
