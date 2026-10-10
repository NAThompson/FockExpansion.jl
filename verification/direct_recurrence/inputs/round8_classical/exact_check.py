"""Exact (sympy) proof that the compact classical part equals the round-2 classical part term by term.
Old: round2/classical_formula_coeffs.pkl + round3/E0_polynomial.pkl.  New: compact_classical.py (transcribed symbolically here)."""
from load import *
pi = sp.pi
dl_ = r1 - r2; sg_ = r1 + r2
def Sm(f, m):
    for _ in range(m % 4): f = f.subs({dl: -s, s: dl}, simultaneous=True)
    return f
kz = (dl*(-9*s**2 + 6*s*xi + 8*xi**2 - 8) - eta*(5 - 2*xi**2))/(144*pi)
kz2 = -5*(xi**2 - 1)*(dl + eta)/(32*pi)
nu = {0: (dl - eta)*(dl + eta), 2: (dl - eta)*(dl + eta), 1: (9*s + eta)*(s - eta), -1: (9*s - eta)*(s + eta)}
def Knew(m, sgn):
    f = Z*Sm(kz, m) + Z**2*(Sm(kz2, m) + dl*nu[m]/(36*pi))
    return f.subs(eta, sgn*eta)
enew = lambda sgn: (Z*(7*s*dl*eta + 9*s**3 + 8*s*xi**2 - 26*s + 10*xi**3 - 12*xi)/72 + Z**2*s*(137*dl*eta + 64*s**2 - 154*xi**2 + 26)/288).subs(eta, sgn*eta)
def F(e):   # to Fock normal form in r1, r2, xi, eta
    return reduce_fock(sp.expand(sp.sympify(e).subs({dl: dl_, s: sg_}, simultaneous=True)))
def oldF(c): return sp.expand(sum(to_fock(c.coeff(Z, k))*Z**k for k in (1, 2)))
ok = True
def chk(name, old, new):
    global ok
    dif = sp.expand(F(old) - F(new)); good = (dif == 0); ok &= good
    print('%-28s %s' % (name, 'exact 0' if good else dif))
pairs = {(1, 5): 0, (1, 3): -1, (1, 7): 1, (1, 1): 2, (-1, 1): 0, (-1, 3): 1, (-1, 5): 2, (-1, 7): -1}
for (sg, k), m in pairs.items():
    chk('c^%s_%d vs K_%d(%s eta)' % ('+' if sg > 0 else '-', k, m, '+' if sg > 0 else '-'), oldF(cpm(sg, k)), Knew(m, sg))
chk('c_lam^+ = e(+) - 2pi K_2(+)', oldF(coef("('L_1_1_1', 'None')")), enew(1) - 2*pi*Knew(2, 1))
chk('c_lam^- = e(-) - 2pi K_2(-)', oldF(coef("('L_1_-1_5', 'None')")), enew(-1) - 2*pi*Knew(2, -1))
chk('c_9', oldF(coef('Lam((ps-pi)/2)/2')), -eta*(Z*(2*xi**2 - 5)/(9*pi) - 5*Z**2*(xi**2 - 1)/(2*pi)))
chk('c_a', oldF(coef('2Lam(al/2)')), Z**2*r2*(8*r2**2 + xi**2 - 5)/(9*pi))
chk('-c_b', -oldF(coef('2Lam((al-pi)/2)')), Z**2*r1*(8*r1**2 + xi**2 - 5)/(9*pi))
# angle polynomial: old poly(alpha, psi) == q_a2 a^2 + q_a a + q_p2 p^2 + q_p p + q0, with q0 moved into E0
cA, cP, cAA, cPP = [coef("('%s', '')" % n) for n in ('AL', 'PS', 'AL^2', 'PS^2')]
qa2 = Z*xi*(2*xi**2 + 1)/(144*pi) - 5*Z**2*xi*(xi**2 - 1)/(32*pi)
qa = Z*7*xi*dl*s/144 + Z**2*(xi*dl*s/64 + dl*(xi**2 - 1)/(9*pi))
qp2 = Z**2*s*(4*s**2 + xi**2 - 9)/(36*pi)
qp = -Z*eta*(5*xi**2 - 4)/(36*pi) + Z**2*s*xi*eta/(18*pi)
chk('q_a2 = c_alal', oldF(cAA), qa2)
chk('q_p2 = c_psps', oldF(cPP), qp2)
chk('q_a = c_al + pi c_alal', oldF(sp.expand(cA + pi*cAA)), qa)
chk('q_p = c_ps + 2pi c_psps', oldF(sp.expand(cP + 2*pi*cPP)), qp)
q0 = oldF(sp.expand(pi/2*cA + pi*cP + pi**2/4*cAA + pi**2*cPP))
E, sym = E0_fock(); Gs, L2s = sym['G'], sym['log2']
E0p_new = (xi*(3 - 2*xi**2)/72
    + Z*(-xi*(5*xi**2 - 6)*(L2s + 4*Gs/pi)/72 + pi*xi*(xi**2 - 1)/36 + xi*(91*xi**2 - 106)/(144*pi)
         + (4*s**3 + 4*s**2*xi - 32*s*xi**2 + 4*s - 27*xi**3 + 26*xi)/288)
    + Z**2*(s*(xi**2 - 1)*(L2s + 4*Gs/pi)/12 - 53*s*(xi**2 - 1)/(72*pi)
            + pi*(16*s**3 - 4*s*xi**2 - 28*s - 45*xi**3 + 45*xi)/576
            - (503*s**3 - 375*s**2*xi + 384*s*xi**2 - 1212*s - 656*xi**3 + 960*xi)/864)
    - Z**3*s*(5*s**2 - 3)/36)
chk("E0' = E0 + q0", sp.expand(E + q0), E0p_new)
print('ALL EXACT' if ok else 'MISMATCH')
