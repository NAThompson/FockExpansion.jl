"""Classical order-3 block in distance coordinates: Gamma table from Codex's exact (t,q) connection (codex/phase7/connections.json).
dX/dr_k = Gamma_k X with Gamma_k = i t (d_k alpha/2) M_t + i q (d_k p0/2) M_q, t = (1-i)(r1+i r2)/(sqrt2 R), q = (c - i P)/(sqrt2 R).
Entries are split as a + b R + c P + d R P with a..d rational in (r1, r2, r12); imaginary parts must vanish (atoms are real)."""
import json, sys
import sympy as sp
r1, r2, c, R, P, I = sp.symbols('r1 r2 r12 R P', positive=True) + (sp.I,)
RV = (r1, r2, c)
Q = r1**2 + r2**2; P2 = 2 * Q - c**2
con = json.load(open('../codex/phase7/connections.json'))['classical']
ATOMS = con['atoms']; n = len(ATOMS)
tS, qS, piS = sp.symbols('t q pi')
def parse(s): return sp.sympify(s, locals={'t': tS, 'q': qS, 'pi': sp.pi, 'I': sp.I})
Mt = {(e[0], e[1]): parse(e[2]) for e in con['dt']['entries']}
Mq = {(e[0], e[1]): parse(e[2]) for e in con['dq']['entries']}
tval = (1 - I) * (r1 + I * r2) / (sp.sqrt(2) * R)
qval = (c - I * P) / (sp.sqrt(2) * R)
dalpha2 = [-r2 / Q, r1 / Q, 0]                     # d(alpha/2) = d atan(r2/r1)
dp02 = [(P * 0 - c * (2 * r1) / P) / (2 * Q), (P * 0 - c * (2 * r2) / P) / (2 * Q), (P - c * (-c) / P) / (2 * Q)]
# d atan(c/P) = (P dc - c dP)/(P^2 + c^2), dP = (2 r1 dr1 + 2 r2 dr2 - c dc)/P, P^2 + c^2 = 2Q

def reduce_rad(expr):
    """polynomial in R, P -> dict {(i,j): coeff} for R^i P^j, i,j in {0,1}, using R^2 = Q, P^2 = P2."""
    poly = sp.Poly(sp.expand(expr), R, P)
    out = {}
    for (i, j), cf in poly.terms():
        k = (i % 2, j % 2)
        out[k] = sp.expand(out.get(k, 0) + cf * Q**(i // 2) * P2**(j // 2))
    return out

def split4(expr):
    num, den = sp.fraction(sp.together(expr))
    num, den = sp.expand(num), sp.expand(den)
    # make the denominator real: multiply by its complex conjugate (I -> -I; R, P real)
    dc = den.subs(I, -I)
    num, den = sp.expand(num * dc), sp.expand(den * dc)
    for sym in (P, R):
        dd = {k: v for k, v in reduce_rad(den).items() if v != 0}
        if all((k[0] if sym == R else k[1]) == 0 for k in dd): continue
        conj = sum(cf * (R**i) * (P**j) * (-1 if (sym == R and i) or (sym == P and j) else 1) for (i, j), cf in dd.items())
        num, den = sp.expand(num * conj), sp.expand(den * conj)
    dd = {k: v for k, v in reduce_rad(den).items() if v != 0}
    assert set(dd) <= {(0, 0)}, list(dd)
    den = dd[(0, 0)]
    out = {}
    for k, cf in reduce_rad(num).items():
        v = sp.cancel(cf / den)
        if v != 0:
            assert sp.im(sp.expand(v)) == 0 or sp.expand(v).has(I) is False, ('complex', v)
            out[k] = v
    return out

def gamma():
    G = [dict() for _ in range(3)]
    for k in range(3):
        for (i, j), m in Mt.items():
            G[k][i, j] = G[k].get((i, j), 0) + I * tval * dalpha2[k] * m.subs({tS: tval, qS: qval})
        for (i, j), m in Mq.items():
            G[k][i, j] = G[k].get((i, j), 0) + I * qval * dp02[k] * m.subs({tS: tval, qS: qval})
    out = [dict() for _ in range(3)]
    for k in range(3):
        for ij, e in G[k].items():
            s4 = split4(e)
            for rad, v in s4.items():
                v = sp.nsimplify(v) if False else v
                if sp.simplify(sp.im(v.subs({r1: sp.Rational(3, 5), r2: sp.Rational(2, 7), c: sp.Rational(1, 2)}))) != 0:
                    raise ValueError(f'complex coefficient {ij} {rad}')
                out[k][ij + (rad,)] = sp.factor(v)
    return out

if __name__ == '__main__':
    G = gamma()
    import pickle; pickle.dump(G, open('classical_gamma.pkl', 'wb'))
    for k in range(3):
        rads = {}
        for (i, j, rad), v in G[k].items(): rads[rad] = rads.get(rad, 0) + 1
        degs = [max(sp.Poly(sp.fraction(v)[0], *RV).total_degree(), sp.Poly(sp.fraction(v)[1], *RV).total_degree()) for v in G[k].values()]
        print(f'd/d{RV[k]}: {len(G[k])} nonzero (atom,atom,radical) entries, radical sectors {rads}, max num/den degree in distances {max(degs)}')
