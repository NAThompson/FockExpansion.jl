"""Laguerre recurrences for the Pekeris operator and its log(u+v+w) sectors.

Psi = exp(-sigma/2) * sum_k log(sigma)^k F_k,  sigma = u+v+w,
F_k = sum A^{(k)}_{lmn} L_l(u) L_m(v) L_n(w).
Sector j equation (multiplied by sigma^2):
  sigma^2 P F_j + (j+1) C1 F_{j+1} + (j+1)(j+2) C2 F_{j+2} = 0,
  C1 = 2 sigma sum_j m_j d_j + (sigma sum_i b_i - sum_ij M_ij),  C2 = sum_ij M_ij.
"""
import itertools
import json
import sympy as sp
from pekeris_operator import pekeris_operator, u, v, w, eps, E, Z, lam

l, m, n = sp.symbols('l m n', integer=True)
X = (u, v, w)
IDX = (l, m, n)
sigma = u + v + w


def op_scaled():
    """P multiplied by eps (polynomial in eps)."""
    return {k: sp.expand(eps * c) for k, c in pekeris_operator().items()}


def split_operator(P):
    """P = sum_ij M_ij d_i d_j + sum_i b_i d_i + c (M symmetric)."""
    M = sp.zeros(3, 3)
    b = [0, 0, 0]
    c = 0
    for k, coef in P.items():
        if sum(k) == 2:
            nz = [i for i in range(3) for _ in range(k[i])]
            i, j = nz
            if i == j:
                M[i, i] += coef
            else:
                M[i, j] += coef / 2
                M[j, i] += coef / 2
        elif sum(k) == 1:
            b[k.index(1)] += coef
        else:
            c += coef
    return M, b, c


def compose_log_ops(P):
    M, b, c = split_operator(P)
    msum = [sp.expand(sum(M[i, j] for i in range(3))) for j in range(3)]
    trM = sp.expand(sum(M))
    C1 = {}
    for j in range(3):
        key = tuple(1 if t == j else 0 for t in range(3))
        C1[key] = sp.expand(2 * sigma * msum[j])
    C1[(0, 0, 0)] = sp.expand(sigma * sum(b) - trM)
    C2 = {(0, 0, 0): trM}
    S2P = {k: sp.expand(sigma**2 * cf) for k, cf in P.items()}
    return S2P, C1, C2


# ---------- one-dimensional Laguerre actions ----------
def xmul(d, nn):
    out = {}
    for (typ, j), c in d.items():
        q = nn + j
        if typ == 'L':
            terms = [(('L', j + 1), -(q + 1)), (('L', j), 2 * q + 1), (('L', j - 1), -q)]
        else:  # x L'_q = q L_q - q L_{q-1}
            terms = [(('L', j), q), (('L', j - 1), -q)]
        for key, cc in terms:
            out[key] = sp.expand(out.get(key, 0) + c * cc)
    return {k: c for k, c in out.items() if c != 0}


def action_1d(a, k, nn):
    """x^a d^k L_nn as {('L'|'D', shift): coeff}."""
    if k == 0:
        d = {('L', 0): sp.Integer(1)}
    elif k == 1:
        if a >= 1:
            d = {('L', 0): nn, ('L', -1): -nn}
            a -= 1
        else:
            d = {('D', 0): sp.Integer(1)}
    elif k == 2:
        if a >= 2:
            xl1 = {('L', 0): nn, ('L', -1): -nn}  # x L'
            t = xmul(xl1, nn)
            for key, c in xl1.items():
                t[key] = sp.expand(t.get(key, 0) - c)
            t2 = xmul({('L', 0): nn}, nn)
            for key, c in t2.items():
                t[key] = sp.expand(t.get(key, 0) - c)
            d = {kk: c for kk, c in t.items() if c != 0}
            a -= 2
        elif a == 1:
            d = {('L', -1): -nn, ('D', 0): sp.Integer(-1)}
            a -= 1
        else:
            raise ValueError('bare second derivative')
    else:
        raise ValueError
    for _ in range(a):
        d = xmul(d, nn)
    return d


def operator_to_stencil(Op):
    """Map an operator dict to the action on L_l L_m L_n:
    {((t1,s1),(t2,s2),(t3,s3)): poly(l,m,n)}; then check D terms cancel."""
    out = {}
    for der, coef in Op.items():
        poly = sp.Poly(coef, *X)
        for (a1, a2, a3), cc in poly.terms():
            acts = [action_1d(a, k, nn) for a, k, nn in zip((a1, a2, a3), der, IDX)]
            for (k1, c1), (k2, c2), (k3, c3) in itertools.product(*[a.items() for a in acts]):
                key = (k1, k2, k3)
                out[key] = out.get(key, 0) + cc * c1 * c2 * c3
    out = {k: sp.expand(c) for k, c in out.items()}
    out = {k: c for k, c in out.items() if c != 0}
    bad = [k for k in out if any(t == 'D' for t, _ in k)]
    assert not bad, f'non-cancelling bare derivative terms: {bad[:5]}'
    return out


def stencil_to_recurrence(st):
    """Coefficient of L_{l}L_{m}L_{n} in Op[sum A]: sum_delta c(l-delta) A(l-delta).
    Returns {offset (as A index shift = -delta): coeff poly in (l,m,n)}."""
    rec = {}
    for key, c in st.items():
        delta = tuple(s for _, s in key)
        sub = {l: l - delta[0], m: m - delta[1], n: n - delta[2]}
        off = tuple(-d for d in delta)
        rec[off] = sp.expand(rec.get(off, 0) + c.subs(sub, simultaneous=True))
    return {k: c for k, c in rec.items() if c != 0}


def summarize(rec, name):
    offs = sorted(rec)
    N = [sum(o) for o in offs]
    print(f'{name}: {len(offs)} offsets; shell offsets dN in [{min(N)}, {max(N)}];',
          'max|offset_i| =', max(max(abs(x) for x in o) for o in offs))


if __name__ == '__main__':
    P = op_scaled()
    st = operator_to_stencil(P)
    rec_general = stencil_to_recurrence(st)
    summarize(rec_general, 'P (general eps)')
    # Pekeris choice eps^2 = -E
    rec_pek = {k: sp.expand(c.subs(E, -eps**2)) for k, c in rec_general.items()}
    rec_pek = {k: c for k, c in rec_pek.items() if c != 0}
    summarize(rec_pek, 'P (eps^2=-E)')
    S2P, C1, C2 = compose_log_ops(P)
    print('sum M =', sp.factor(C2[(0, 0, 0)]))
    out = {}
    for name, Op in [('sigma2P', S2P), ('C1', C1), ('C2', C2)]:
        r = stencil_to_recurrence(operator_to_stencil(Op))
        summarize(r, name + ' (general eps)')
        rp = {k: sp.expand(c.subs(E, -eps**2)) for k, c in r.items()}
        rp = {k: c for k, c in rp.items() if c != 0}
        summarize(rp, name + ' (eps^2=-E)')
        out[name] = {str(k): str(c) for k, c in rp.items()}
        out[name + '_general'] = {str(k): str(c) for k, c in r.items()}
    out['P'] = {str(k): str(c) for k, c in rec_pek.items()}
    out['P_general'] = {str(k): str(c) for k, c in rec_general.items()}
    with open('recurrence.json', 'w') as f:
        json.dump(out, f, indent=1)
