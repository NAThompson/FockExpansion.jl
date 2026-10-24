"""Machine check: apply sigma^2 * P directly (sympy differentiation, exact log) to
exp(-sigma/2) [F0 + log(sigma) F1 + log(sigma)^2 F2] with random finite Laguerre
coefficient arrays, expand each log^j coefficient in the Laguerre basis, and compare
with the recurrence (sector j: sigma2P A_j + (j+1) C1 A_{j+1} + (j+1)(j+2) C2 A_{j+2})."""
import itertools
import random
import sympy as sp
from pekeris_operator import pekeris_operator, u, v, w, eps, E, Z, lam
from recurrence import (op_scaled, compose_log_ops, operator_to_stencil,
                        stencil_to_recurrence, l, m, n)

X = (u, v, w)
sigma = u + v + w
T = sp.Symbol('T')


def to_laguerre(poly_expr):
    out = {}
    p = sp.Poly(sp.expand(poly_expr), *X)
    for (a, b, c), cf in p.terms():
        exps = []
        for k in (a, b, c):
            exps.append({j: sp.factorial(k) * (-1)**j * sp.binomial(k, j) for j in range(k + 1)})
        for (i, ci), (j, cj), (k, ck) in itertools.product(exps[0].items(), exps[1].items(), exps[2].items()):
            out[(i, j, k)] = out.get((i, j, k), 0) + cf * ci * cj * ck
    return {k: sp.Rational(c) for k, c in out.items() if c != 0}


def apply_rec(rec, A, params):
    """(rec A)(idx) = sum_off rec[off](idx) A(idx + off) for all idx where nonzero."""
    out = {}
    targets = set()
    for (i, j, k) in A:
        for off in rec:
            targets.add((i - off[0], j - off[1], k - off[2]))
    for t in targets:
        if min(t) < 0:
            continue
        s = 0
        for off, c in rec.items():
            src = (t[0] + off[0], t[1] + off[1], t[2] + off[2])
            if src in A:
                s += c.subs({l: t[0], m: t[1], n: t[2]}).subs(params) * A[src]
        if s != 0:
            out[t] = sp.Rational(s)
    return out


def main(seed=0, pekeris_relation=False):
    random.seed(seed)
    params = {eps: sp.Rational(3, 2), Z: 2, lam: 1,
              E: -sp.Rational(9, 4) if pekeris_relation else -sp.Rational(7, 3)}
    P = op_scaled()
    S2P, C1, C2 = compose_log_ops(P)
    recs = {name: stencil_to_recurrence(operator_to_stencil(Op))
            for name, Op in [('S2P', S2P), ('C1', C1), ('C2', C2)]}
    A = []
    for k in range(3):
        A.append({idx: sp.Rational(random.randint(-9, 9), random.randint(1, 5))
                  for idx in itertools.product(range(3), repeat=3) if random.random() < 0.5})
    F = [sum(c * sp.laguerre(i, u) * sp.laguerre(j, v) * sp.laguerre(kk, w) for (i, j, kk), c in Ak.items()) for Ak in A]
    G = F[0] + sp.log(sigma) * F[1] + sp.log(sigma)**2 * F[2]
    # apply P (dict of derivative orders) directly
    PG = 0
    for (a, b, c), coef in P.items():
        d = G
        for x, k in zip(X, (a, b, c)):
            if k:
                d = sp.diff(d, x, k)
        PG += coef.subs(params) * d
    expr = sp.expand(sp.expand(sigma**2 * PG).subs(sp.log(sigma), T))
    expr = sp.expand(sp.cancel(expr))
    poly = sp.Poly(expr, T)
    ok = True
    for j in range(3):
        direct = to_laguerre(poly.coeff_monomial(T**j))
        pred = {}
        terms = [(recs['S2P'], A[j], 1)]
        if j + 1 < 3:
            terms.append((recs['C1'], A[j + 1], j + 1))
        if j + 2 < 3:
            terms.append((recs['C2'], A[j + 2], (j + 1) * (j + 2)))
        for rec, Ak, w8 in terms:
            for key, val in apply_rec(rec, Ak, params).items():
                pred[key] = pred.get(key, 0) + w8 * val
        pred = {k: c for k, c in pred.items() if c != 0}
        diff = set(direct) ^ set(pred)
        bad = [k for k in set(direct) & set(pred) if sp.simplify(direct[k] - pred[k]) != 0]
        print(f'seed {seed}, eps^2=-E: {pekeris_relation}, sector {j}: {len(direct)} Laguerre coefficients;',
              'support mismatch', len(diff), 'value mismatch', len(bad))
        for k in bad[:4]: print("   ", k, direct[k], pred[k])
        ok &= not diff and not bad
    return ok


if __name__ == '__main__':
    import sys
    allok = all([main(0, False)] if 'quick' in sys.argv else [main(0, False), main(1, True), main(2, False)])
    print('ALL OK' if allok else 'FAILED')
