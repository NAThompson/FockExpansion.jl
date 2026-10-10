"""Compact classical part of the general-angle psi_{3,0} (round 8).  E = a21 = 0 part; see FORMULA_classical_compact.txt.
Variables: r1 = cos(alpha/2), r2 = sin(alpha/2), sigma = r1 + r2, delta = r1 - r2, xi = r12/R = sqrt2 sin(psi/2), eta = sqrt(2 - xi^2),
a = alpha - pi/2, p = psi - pi, w_s = (a + s p)/4.  Functions: L(u) = int_0^u x cot x dx, T(w) = int_0^w x tan x dx."""
import mpmath as mp

def L(u): return u*mp.log(abs(2*mp.sin(u))) + mp.clsin(2, 2*u)/2
def T(w): return -w*mp.log(abs(2*mp.cos(w))) + mp.clsin(2, mp.pi - 2*w)/2
def lc(w): return mp.log(abs(2*mp.cos(w)))

def K(m, d, sg, xi, eta, Z):
    """coefficient of 4 L(w + m pi/4) (m = -1, 0, 1) and of -4 T(w) (m = 2); (d, sg) = (delta, sigma)."""
    pi = mp.pi
    for _ in range(m % 4): d, sg = -sg, d                         # S: (delta, sigma) -> (-sigma, delta), i.e. alpha -> alpha + pi
    kz = (d*(-9*sg**2 + 6*sg*xi + 8*xi**2 - 8) - eta*(5 - 2*xi**2))/(144*pi)
    kz2 = -5*(xi**2 - 1)*(d + eta)/(32*pi)
    return Z*kz + Z**2*kz2

def N(m, d, sg, xi, eta):
    """non-covariant part of the Z^2 coefficients"""
    nu = {0: xi**2 - sg**2, 2: xi**2 - sg**2, 1: 9*sg**2 + xi**2 - 2 - 8*sg*eta, -1: 9*sg**2 + xi**2 - 2 + 8*sg*eta}[m]
    return d*nu/(36*mp.pi)

def e_lam(d, sg, xi, eta, Z):
    return Z*(7*sg*d*eta + 9*sg**3 + 8*sg*xi**2 - 26*sg + 10*xi**3 - 12*xi)/72 + Z**2*sg*(137*d*eta + 64*sg**2 - 154*xi**2 + 26)/288

def CL_compact(al, ps, Z):
    al, ps, Z = mp.mpf(al), mp.mpf(ps), mp.mpf(Z); pi = mp.pi
    r1, r2 = mp.cos(al/2), mp.sin(al/2); sg, d = r1 + r2, r1 - r2
    xi, eta0 = mp.sqrt(2)*mp.sin(ps/2), mp.sqrt(2)*mp.cos(ps/2)
    a, p = al - pi/2, ps - pi
    G, l2 = mp.catalan, mp.log(2)
    v = 0
    for s in (1, -1):
        w = (a + s*p)/4; eta = s*eta0
        kk = {m: K(m, d, sg, xi, eta, Z) + Z**2*N(m, d, sg, xi, eta) for m in (-1, 0, 1, 2)}
        v += sum(4*kk[m]*L(w + m*pi/4) for m in (-1, 0, 1)) - 4*kk[2]*T(w) + e_lam(d, sg, xi, eta, Z)*lc(w)
    # one-variable functions
    c9 = -eta0*(Z*(2*xi**2 - 5)/(9*pi) - 5*Z**2*(xi**2 - 1)/(2*pi))
    v += c9*L(p/2)/2
    v += 2*Z**2/(9*pi)*(r2*(8*r2**2 + xi**2 - 5)*L(al/2) + r1*(8*r1**2 + xi**2 - 5)*L((pi - al)/2))
    # angle polynomial (centred)
    qa2 = Z*xi*(2*xi**2 + 1)/(144*pi) - 5*Z**2*xi*(xi**2 - 1)/(32*pi)
    qa = Z*7*xi*d*sg/144 + Z**2*(xi*d*sg/64 + d*(xi**2 - 1)/(9*pi))
    qp2 = Z**2*sg*(4*sg**2 + xi**2 - 9)/(36*pi)
    qp = -Z*eta0*(5*xi**2 - 4)/(36*pi) + Z**2*sg*xi*eta0/(18*pi)
    v += qa2*a**2 + qa*a + qp2*p**2 + qp*p
    # algebraic remainder E0' (polynomial in sigma, xi)
    v += E0p(sg, xi, Z, G, l2)
    return v

def E0p(sg, xi, Z, G, l2):
    pi = mp.pi
    e0 = xi*(3 - 2*xi**2)/72
    e1 = (-xi*(5*xi**2 - 6)*(l2 + 4*G/pi)/72 - pi*xi*(1 - xi**2)/36 + xi*(91*xi**2 - 106)/(144*pi)
          + (4*sg**3 + 4*sg**2*xi - 32*sg*xi**2 + 4*sg - 27*xi**3 + 26*xi)/288)
    e2 = (-sg*(1 - xi**2)*(l2 + 4*G/pi)/12 + 53*sg*(1 - xi**2)/(72*pi)
          + pi*(16*sg**3 - 4*sg*xi**2 - 28*sg - 45*xi**3 + 45*xi)/576
          - (503*sg**3 - 375*sg**2*xi + 384*sg*xi**2 - 1212*sg - 656*xi**3 + 960*xi)/864)
    e3 = -sg*(5*sg**2 - 3)/36
    return e0 + Z*e1 + Z**2*e2 + Z**3*e3
