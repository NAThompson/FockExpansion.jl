# Symbolic check of Liverts-Barnea (2015) Table I entries for psi_{4,1}^{(1)}, psi_{4,1}^{(3)}
# using Lambda^2 in the coordinates a = cos(alpha), v = sin(alpha) cos(theta).
from sympy import *
a, v, E = symbols('a v E', real=True)
s = sqrt(1-a**2)                  # sin(alpha), alpha in [0, pi]
xi = sqrt(1-v)                    # r12/R
vs = sqrt(1+s)                    # varsigma = cos(a/2)+sin(a/2)
B = (pi-2)/(3*pi)
def Lam2(f):
    return -4*((1-a**2)*diff(f,a,2)+(1-v**2)*diff(f,v,2)-2*a*v*diff(f,a,v)-3*a*diff(f,a)-3*v*diff(f,v))
V0, V1 = 1/xi, 2*vs/s
p21_1 = -B*v
p31_1 = B*xi*(5*xi**2-6)/12
p31_2 = B*vs*v/2
p42_2 = (pi-2)*(5*pi-14)/(180*pi**2)*(1-2*(s**2-v**2))   # sin^2 a sin^2 t = s^2 - v^2
h1 = -2*V0*p31_1+2*E*p21_1
h3 = 2*V1*p31_2
p41_1 = (pi-2)/(2880*pi)*(3*(32*E-15)-8*(12*E-5)*xi**2)
p41_3 = -(pi-2)/(120*pi)*(4+5*s)*v
print("j=1 residual:", simplify(Lam2(p41_1)-32*p41_1-h1))
print("j=3 residual:", simplify(Lam2(p41_3)-32*p41_3-h3))
print("psi42 harmonic:", simplify(Lam2(p42_2)-32*p42_2))
