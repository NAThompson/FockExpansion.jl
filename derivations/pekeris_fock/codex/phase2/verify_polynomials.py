#!/usr/bin/env python3
"""Exact homogenization checks against the explicit repo angular formulas."""
import sys
sys.dont_write_bytecode=True
import sympy as s
a,b,d=s.symbols('a b d',positive=True)
Z,E,A=s.symbols('Z E a21')
q=a*a+b*b; R=s.sqrt(q); S=a+b; D=q-d*d
xi=d/R; sig=S/R; v=D/q; sa=2*a*b/q; ca=(a*a-b*b)/q
pi=s.pi
cases=[]
def check(name,angular,degree,poly):
    assert s.simplify(s.expand(R**degree*angular-poly))==0,name
    assert s.Poly(poly,a,b,d).total_degree()==degree,name
    cases.append(name)
check('psi10',xi/2-Z*sig,1,d/2-Z*S)
check('psi21',-Z*(pi-2)*v/(3*pi),2,-Z*(pi-2)*D/(3*pi))
check('psi20 exposed polynomial',(1-2*E)/12-Z*sig*xi/3+Z**2*(s.Rational(1,3)+sa/2)+(A-Z*(pi+4)/(9*pi))*v,2,
    (1-2*E)*q/12-Z*S*d/3+Z**2*(q/3+a*b)+(A-Z*(pi+4)/(9*pi))*D)
check('psi31',Z*(pi-2)*xi*(5*xi**2-6)/(36*pi)+Z**2*(pi-2)*sig*v/(6*pi),3,
    Z*(pi-2)*d*(5*d*d-6*q)/(36*pi)+Z**2*(pi-2)*S*D/(6*pi))
H40=4*(a*a-b*b)**2-q*q; H42=(3*D*D-4*a*a*b*b)/2
W=q*q-8*a*a*b*b+2*D*D
check('Y40',4*ca*ca-1,4,H40)
check('Y42',(3*v*v-sa*sa)/2,4,H42)
check('psi42',Z**2*(pi-2)*(5*pi-14)*(1-2*(sa*sa-v*v))/(180*pi*pi),4,
      Z**2*(pi-2)*(5*pi-14)*W/(180*pi*pi))
check('psi41 Z1',(pi-2)*(3*(32*E-15)-8*(12*E-5)*xi*xi)/(2880*pi),4,
      (pi-2)*q*(3*(32*E-15)*q-8*(12*E-5)*d*d)/(2880*pi))
check('psi41 Z3',-(pi-2)*(4+5*sa)*v/(120*pi),4,-(pi-2)*(4*q+10*a*b)*D/(120*pi))
check('psi30 state',E*(Z*sig*(2+a*b/q)/18-(6-xi*xi)*xi/72)-A*(Z*sig*v/2-(6-5*xi*xi)*xi/12),3,
      E*(Z*S*(2*q+a*b)/18-(6*q-d*d)*d/72)-A*(Z*S*D/2-(6*q-5*d*d)*d/12))
check('psi40 exposed Z0',(12*E*E-11*E+6*A+1)/1152-(72*E*A+E-30*A-2)*v/720,4,
      (12*E*E-11*E+6*A+1)*q*q/1152-(72*E*A+E-30*A-2)*q*D/720)
check('psi40 exposed Z4',sa/36+s.Rational(7,288)+2*(4*ca*ca-1)/(135*pi),4,
      a*b*q/18+7*q*q/288+2*H40/(135*pi))
print('PASS:',len(cases),'exact source-formula homogenizations:',', '.join(cases))
