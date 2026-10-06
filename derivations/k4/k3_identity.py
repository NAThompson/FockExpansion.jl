# Numerical check of the ψ₃₀ panel identity: for Q = 1, t², the Â-coefficient R₁(t) of
# (Λ²-21)(QÂ/√P) should equal ∂_t(H/√P) with H a cubic in t. Fit H at fixed (a, v).
from sympy import *
import mpmath as mp
a, v, t, X = symbols('a v t X', real=True)
def Lam2(f):
    return -4*((1-a**2)*diff(f,a,2)+(1-v**2)*diff(f,v,2)-2*a*v*diff(f,a,v)-3*a*diff(f,a)-3*v*diff(f,v))
P = 1+a+2*v*t**2+(1-a)*t**4
m = sqrt(2*(1-v))*t
Ahat = 2/pi*atan(sqrt(P)/m)
mp.mp.dps = 30
import sys; KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
av, vv = mp.mpf('0.3'), mp.mpf('-0.2')
b = a; y2 = (1-v)/2
Qz1 = sqrt(2)*(b**2-1)/72 + sqrt(2)*(b-1)*(-18*b-32*y2+16)/288*t**2      # Z¹ part of c₀ + c₂t²
Qz2 = -sqrt(2)*(9*b**2-18)/576 + sqrt(2)*(b-1)*(90*y2-45)/288*t**2         # Z² part
for Q in (Qz1, Qz2):
    expr = Lam2(Q*Ahat/sqrt(P)) - 21*Q*Ahat/sqrt(P)
    expr = expand(expr.subs(atan(sqrt(P)/m), pi*X/2))
    R1 = lambdify((a, v, t), expr.coeff(X, 1), 'mpmath')
    R0 = lambdify((a, v, t), expr.coeff(X, 0), 'mpmath')
    ts = [mp.mpf(k)/40 for k in range(1, 40)]
    Pf = lambda tt: 1+av+2*vv*tt**2+(1-av)*tt**4
    dPf = lambda tt: 4*vv*tt+4*(1-av)*tt**3
    rows = [[(k*tt**(k-1) if k > 0 else 0)/mp.sqrt(Pf(tt)) - tt**k*dPf(tt)/(2*Pf(tt)**mp.mpf(1.5)) for k in range(KMAX)] for tt in ts]
    A = mp.matrix(rows); b = mp.matrix([R1(av, vv, tt) for tt in ts])
    h = mp.lu_solve(A.T*A, A.T*b)
    r = A*h - b
    res = max(abs(r[i]) for i in range(len(ts)))
    print(f"weight {'Z1' if Q==Qz1 else 'Z2'}: H coefficients {[mp.nstr(x, 10) for x in h]}, max residual {mp.nstr(res, 3)}")
    print("   Â-free part at t = 0.5:", mp.nstr(R0(av, vv, mp.mpf('0.5')), 12))
