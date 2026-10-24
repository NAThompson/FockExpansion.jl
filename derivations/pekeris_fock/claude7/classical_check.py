"""Check the distance-coordinate classical Gamma table against finite differences of the actual 18 atoms (mpmath)."""
import pickle, mpmath as m, sympy as sp
from classical_gamma import ATOMS, r1, r2, c, R, P
m.mp.dps = 50
G = pickle.load(open('classical_gamma.pkl', 'rb'))
cl2 = lambda x: m.clsin(2, x)
L = lambda z: z * m.log(abs(2 * m.sin(z))) + cl2(2 * z) / 2
T = lambda z: -z * m.log(abs(2 * m.cos(z))) + cl2(m.pi - 2 * z) / 2
lc = lambda z: m.log(abs(2 * m.cos(z)))
def atoms(x, y, z):
    Rr = m.sqrt(x * x + y * y); Pp = m.sqrt(2 * Rr * Rr - z * z)
    al = 2 * m.atan(y / x); p0 = 2 * m.atan(z / Pp)
    a = al - m.pi / 2; p = p0 - m.pi; wp, wm = (a + p) / 4, (a - p) / 4
    return [L(wp - m.pi / 4), L(wp), L(wp + m.pi / 4), T(wp), L(wm - m.pi / 4), L(wm), L(wm + m.pi / 4), T(wm),
            L(p / 2), L((a + m.pi / 2) / 2), L((m.pi / 2 - a) / 2), lc(wp), lc(wm), a, p, a * a, p * p, m.mpf(1)]
pt = [m.mpf('0.61'), m.mpf('0.37'), m.mpf('0.52')]
A0 = atoms(*pt); Rr = m.sqrt(pt[0]**2 + pt[1]**2); Pp = m.sqrt(2 * Rr**2 - pt[2]**2)
sub = {r1: pt[0], r2: pt[1], c: pt[2]}
worst = 0
for k in range(3):
    h = m.mpf('1e-20'); pp = list(pt); pm = list(pt); pp[k] += h; pm[k] -= h
    Ap, Am = atoms(*pp), atoms(*pm)
    for i in range(18):
        fd = (Ap[i] - Am[i]) / (2 * h)
        an = 0
        for (ii, j, rad), v in G[k].items():
            if ii == i:
                an += m.mpf(str(sp.N(v.subs(sub), 40))) * (Pp if rad == (0, 1) else 1) * A0[j]
        worst = max(worst, abs(fd - an))
print('max |finite difference - Gamma| over 18 atoms x 3 directions:', m.nstr(worst, 5))
