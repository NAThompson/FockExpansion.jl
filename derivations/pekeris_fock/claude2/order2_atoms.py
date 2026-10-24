"""Order-2 closed atom set and its gradient table Gamma (rational coefficients in r1, r2, r12).

Atoms (P = sqrt(2R^2 - r12^2) = |x1+x2|, R^2 = r1^2+r2^2, s_h = r1+r2, d = r1-r2):
  weight 2: Om = logs - DL  (chi = c*Om/(6 pi) + lower weight, c = s_h d / R^2)
  weight 1: al = 2 atan(r2/r1), be = asin((R^2-r12^2)/R^2), beP = P*be,
            lR = log(R/s_h), la = log((s_h+r12)/R), lb = log((d+r12)^2/(r12 P + s_h d))
  weight 0: one = 1, P
grad(atom) = sum_{atom'} Gamma[atom][atom'] * atom', Gamma rational in r.
"""
import sympy as sp

r1, r2, r12, Ps = sp.symbols('r1 r2 r12 P', positive=True)
RV = (r1, r2, r12)
sh, d = r1 + r2, r1 - r2
R2 = r1**2 + r2**2
P2 = 2 * R2 - r12**2
ATOMS = ['Om', 'al', 'beP', 'lR', 'la', 'lb', 'one', 'P']
WEIGHT = {'Om': 2, 'al': 1, 'beP': 1, 'lR': 1, 'la': 1, 'lb': 1, 'one': 0, 'P': 0}


def dP(k):  # dP/dr_k = P * d(P^2)/dr_k / (2 P^2)
    return Ps * sp.diff(P2, RV[k]) / (2 * P2)


def tot(expr, k):
    return sp.diff(expr, RV[k]) + sp.diff(expr, Ps) * dP(k)


def reduceP(expr):
    """Reduce a polynomial in P modulo P^2 = P2: returns (a, b) with expr = a + b P."""
    poly = sp.Poly(sp.expand(expr), Ps)
    a = b = 0
    for (n,), c in poly.terms():
        if n % 2 == 0:
            a += c * P2**(n // 2)
        else:
            b += c * P2**(n // 2)
    return sp.expand(a), sp.expand(b)


def split(expr):
    """expr (rational in r and P) -> (a, b), expr = a + b P, a, b rational in r."""
    num, den = sp.fraction(sp.together(sp.expand(expr)))
    A, B = reduceP(den)
    num = sp.expand(num * (A - B * Ps))
    den = sp.expand(A**2 - B**2 * P2)
    a, b = reduceP(num)
    return sp.cancel(a / den), sp.cancel(b / den)


v = (R2 - r12**2) / R2
L5 = r12 * Ps + sh * d


def omega_beta():
    """beta-coefficients of grad Omega (from atom_closure/atom_reduce: all log letters cancel)."""
    import atom_closure as ac
    from atom_reduce import reduce, idx
    out = []
    for k in range(3):
        e = sp.sympify(reduce(ac.allrows[k])[idx['be']])
        out.append(e.subs({ac.R: sp.sqrt(R2), ac.P: Ps, ac.r1: r1, ac.r2: r2, ac.r12: r12}))
        # sanity: alpha coefficient equals dlog Qa
        ea = sp.sympify(reduce(ac.allrows[k])[idx['al']]).subs({ac.R: sp.sqrt(R2), ac.P: Ps, ac.r1: r1, ac.r2: r2, ac.r12: r12})
        assert sp.simplify(ea - sp.diff(sp.log((sh**2 - r12**2) / (r12**2 - d**2)), RV[k])) == 0
        # no other weight-1 content
        for n in ac.W1:
            if n not in ('al', 'be'):
                assert sp.simplify(sp.sympify(reduce(ac.allrows[k])[idx[n]])) == 0, n
    return out


OMEGA_BETA = None


def gamma_table():
    global OMEGA_BETA
    if OMEGA_BETA is None:
        OMEGA_BETA = omega_beta()
    G = {a: [dict() for _ in range(3)] for a in ATOMS}
    for k in range(3):
        def add(atom, target_pair, expr):
            a, b = split(expr)
            t0, t1 = target_pair
            if a != 0:
                G[atom][k][t0] = sp.cancel(G[atom][k].get(t0, 0) + a)
            if b != 0:
                G[atom][k][t1] = sp.cancel(G[atom][k].get(t1, 0) + b)
        dbe = tot(v, k) * R2 / (r12 * Ps)                 # d beta = dv / g
        add('al', ('one', 'P'), sp.diff(2 * sp.atan(r2 / r1), RV[k]))
        add('beP', ('one', 'P'), Ps * dbe)                 # P d(beta)  (P*P -> P^2)
        add('beP', ('beP', 'beP'), dP(k) / Ps)             # beta dP = (beta P) * dP/P
        add('lR', ('one', 'P'), sp.diff(sp.log(sp.sqrt(R2) / sh), RV[k]))
        add('la', ('one', 'P'), sp.diff(sp.log((sh + r12) / sp.sqrt(R2)), RV[k]))
        add('lb', ('one', 'P'), 2 * sp.diff(d + r12, RV[k]) / (d + r12) - tot(L5, k) / L5)
        add('P', ('P', 'P'), dP(k) / Ps)                   # dP = P * (dP/P)
        # Omega: alpha * dlog((s_h^2 - r12^2)/(r12^2 - d^2)) + beta * B_k  (B_k from atom_closure, reduced)
        Qa = (sh**2 - r12**2) / (r12**2 - d**2)
        add('Om', ('al', 'al'), sp.diff(sp.log(Qa), RV[k]))
        add('Om', ('be', 'beP'), OMEGA_BETA[k])
    return G


# ---------------- numerical atom values (mpmath) ----------------
def atom_values(a, b, c):
    import mpmath as mp
    a, b, c = mp.mpf(a), mp.mpf(b), mp.mpf(c)
    R = mp.sqrt(a**2 + b**2); P = mp.sqrt(2 * R**2 - c**2); S, D = a + b, a - b
    al = 2 * mp.atan(b / a); be = mp.asin((R**2 - c**2) / R**2)
    cl2 = lambda t: mp.clsin(2, t)
    L5 = c * P + S * D
    logs = (be * (2 * mp.log(2 * a * b) - 2 * mp.log(L5))
            + (be + al) * (mp.log(S - c) + mp.log(S + c) - mp.log(2 * a * b))
            + (be - al) * (mp.log(c - D) + mp.log(c + D) - mp.log(2 * a * b)))
    DL = cl2(al - be) + cl2(mp.pi + be - al) - cl2(mp.pi - al - be) - cl2(al + be)
    return {'Om': logs - DL, 'al': al, 'beP': be * P, 'lR': mp.log(R / S), 'la': mp.log((S + c) / R),
            'lb': mp.log((D + c)**2 / L5), 'one': mp.mpf(1), 'P': P}


def check_gamma(G, pt=('0.62', '0.41', '0.53')):
    import mpmath as mp
    mp.mp.dps = 50
    h = mp.mpf('1e-18')
    worst = 0
    for k in range(3):
        pp = [mp.mpf(x) for x in pt]; pm = list(pp); pp[k] += h; pm[k] -= h
        vp, vm, v0 = atom_values(*pp), atom_values(*pm), atom_values(*[mp.mpf(x) for x in pt])
        sub = {r1: mp.mpf(pt[0]), r2: mp.mpf(pt[1]), r12: mp.mpf(pt[2])}
        for a in ATOMS:
            fd = (vp[a] - vm[a]) / (2 * h)
            an = sum(mp.mpf(str(sp.N(cf.subs(sub), 40))) * v0[t] for t, cf in G[a][k].items())
            worst = max(worst, abs(fd - an))
    return worst


if __name__ == '__main__':
    G = gamma_table()
    print('max |finite difference - Gamma| over atoms and directions:', check_gamma(G))
    for a in ATOMS:
        print(a, [sorted(G[a][k].keys()) for k in range(3)])


