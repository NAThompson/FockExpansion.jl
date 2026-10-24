"""Check the reconstructed rational symmetric panel connection against ACTUAL integrals (mpmath quadrature + finite differences)
at a physical panel-1 point. States per family/jet: J, K, S0..S3 (= sum z_i^k J_i over quartic roots), J5..J8, Js, e, f0."""
import json
import mpmath as m
from fractions import Fraction
from panel1_reconstruct import monos
m.mp.dps = 40
con = json.load(open('panel1_connection.json'))

def states(h, Y2):
    b = (h*h - 1) / (h*h + 1)
    c4, c2, c0 = 1 - b, 4*Y2 - 2, 1 + b
    disc = m.sqrt(c2*c2 - 4*c4*c0 + 0j)
    Z2 = (-c2 + disc) / (2*c4); z = m.sqrt(Z2)
    zs = [z, -z, h/z, -h/z, 1j, -1j, 1j*h, -1j*h]
    W = lambda q: m.sqrt(c4*q**4 + c2*q*q + c0)
    out = {}
    for fam in ('A', 'D'):
        for n in range(3):
            ang = (lambda q: m.atan(q)) if fam == 'A' else (lambda q: m.atan(q/h))
            f = lambda q, n=n: ang(q)**n / (m.factorial(n) * W(q))
            J = m.quad(f, [0, 1]); K = m.quad(lambda q: q*f(q), [0, 1])
            Ji = [m.quad(lambda q, zi=zi: f(q)/(q - zi), [0, 1]) for zi in zs]
            base = 1/m.sqrt(c0) if n == 0 else 0
            Js = m.quad(lambda q: (f(q) - base)/q, [0, 1])          # + base*log(1) = 0
            e = f(1); f0 = base
            S = [sum(zs[i]**k * Ji[i] for i in range(4)) for k in range(4)]
            vec = [J, K] + S + Ji[4:] + [Js, e, f0]
            for st, val in enumerate(vec): out[(fam, st, n)] = val
    return out

def coef(key, h, Y2):
    r = con.get(key)
    if r is None: return 0
    mn, md = monos(r['dn']), monos(r['dd'])
    N = sum(Fraction(c) * h**a * Y2**bb for c, (a, bb) in zip(r['num'], mn))
    D = sum(Fraction(c) * h**a * Y2**bb for c, (a, bb) in zip(r['den'], md))
    return N / D

h0, Y0 = m.mpf('1.7'), m.mpf('0.21')       # physical: h = r1/r2 > 0, Y2 = r12^2/(2R^2) in (0,1)
X0 = states(h0, Y0)
worst = 0
for ax, (dh, dY) in ((0, (m.mpf('1e-12'), 0)), (1, (0, m.mpf('1e-12')))):
    Xp, Xm = states(h0 + dh, Y0 + dY), states(h0 - dh, Y0 - dY)
    step = 2 * (dh if ax == 0 else dY)
    for fam in ('A', 'D'):
        for st in range(13):
            for n in range(3):
                fd = (Xp[(fam, st, n)] - Xm[(fam, st, n)]) / step
                r = st * 3 + n
                an = 0
                for c in range(39):
                    re = coef(f'{ax}{fam}|0|{r}|{c}', h0, Y0); im = coef(f'{ax}{fam}|1|{r}|{c}', h0, Y0)
                    if re == 0 and im == 0: continue
                    an += (m.mpf(re) + 1j * m.mpf(im)) * X0[(fam, c // 3, c % 3)]
                err = abs(fd - an) / (1 + abs(fd))
                worst = max(worst, err)
print('max relative |finite difference - connection*state| over 2 axes x 2 families x 39 states:', m.nstr(worst, 4))
