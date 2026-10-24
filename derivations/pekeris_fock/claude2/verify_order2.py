"""Machine check of the order-2 atom recurrence (order2_recurrence.json) by direct operator application.

Psi = exp(-sigma/2) sum_{k=0,1} log(sigma)^k sum_b b(r) F_{k,b}(u,v,w), random small Laguerre arrays.
LHS(pt) = sigma^2 * [-(2v+w)(2u+w)(u+v)] * exp(sigma/2) (H - E) Psi, computed with mpmath numerical derivatives
of the actual functions (Clausen, atan, asin, logs, sqrt) in Hylleraas coordinates.
RHS(pt) = sum_{k,b'} log(sigma)^k b'(pt) G_{k,b'}(pt) / D_{b'}(pt), G from the stencils.
Parameters: eps = 3/2, E = -eps^2, Z = 2, lambda = 1 (as in the derivation)."""
import json
import random
import mpmath as mp
import sympy as sp
import order2_atoms as oa

mp.mp.dps = 60
EPS = mp.mpf(3) / 2
rec = json.load(open('order2_recurrence.json'))
ATOMS = rec['atoms']
l, m, n, eps = sp.symbols('l m n epsilon')
u, v, w = sp.symbols('u v w')


def compile_block(blk):
    out = []
    for off, cf in blk.items():
        f = sp.lambdify((l, m, n), sp.sympify(cf).subs(eps, sp.Rational(3, 2)), 'mpmath')
        out.append((eval(off), f))
    return out


BLK = {name: {eval(k): compile_block(b) for k, b in rec[name].items()} for name in ('sigma2P', 'C1', 'C2')}
DEN = {b: sp.lambdify((u, v, w), sp.sympify(rec['D'][b]).subs(sp.Symbol('epsilon'), sp.Rational(3, 2)), 'mpmath') for b in ATOMS}


def lag(nn, x):
    return mp.laguerre(nn, 0, x)


def psi(r1, r2, r12, A):
    U, V, W = EPS * (r2 + r12 - r1), EPS * (r1 + r12 - r2), 2 * EPS * (r1 + r2 - r12)
    sig = U + V + W
    at = oa.atom_values(r1, r2, r12)
    tot = 0
    for (k, b), arr in A.items():
        F = sum(c * lag(i, U) * lag(j, V) * lag(q, W) for (i, j, q), c in arr.items())
        tot += mp.log(sig)**k * at[b] * F
    return mp.exp(-sig / 2) * tot


def H_minus_E(f, r1, r2, r12):
    E = -EPS**2
    d = lambda *o: mp.diff(f, (r1, r2, r12), o)
    lap = (d(2, 0, 0) + 2 / r1 * d(1, 0, 0) + d(0, 2, 0) + 2 / r2 * d(0, 1, 0) + 2 * d(0, 0, 2) + 4 / r12 * d(0, 0, 1)
           + (r1**2 - r2**2 + r12**2) / (r1 * r12) * d(1, 0, 1) + (r2**2 - r1**2 + r12**2) / (r2 * r12) * d(0, 1, 1))
    V = -2 / r1 - 2 / r2 + 1 / r12
    return -lap / 2 + (V - E) * f(r1, r2, r12)


def rhs(r1, r2, r12, A):
    U, V, W = EPS * (r2 + r12 - r1), EPS * (r1 + r12 - r2), 2 * EPS * (r1 + r2 - r12)
    sig = U + V + W
    at = oa.atom_values(r1, r2, r12)
    total = 0
    for k in (0, 1):
        for bp in ATOMS:
            G = 0
            for name, kk, w8 in (('sigma2P', k, 1), ('C1', k + 1, k + 1), ('C2', k + 2, (k + 1) * (k + 2))):
                for (t, b), stencil in BLK[name].items():
                    if t != bp or (kk, b) not in A:
                        continue
                    arr = A[(kk, b)]
                    # equation index e: sum_off coeff(e) A(e+off); evaluate G(pt) = sum_e (...) L_e(pt)
                    targets = {(i - o[0], j - o[1], q - o[2]) for (i, j, q) in arr for o, _ in stencil}
                    for e in targets:
                        if min(e) < 0:
                            continue
                        val = 0
                        for o, f in stencil:
                            src = (e[0] + o[0], e[1] + o[1], e[2] + o[2])
                            if src in arr:
                                val += f(*e) * arr[src]
                        if val != 0:
                            G += w8 * val * lag(e[0], U) * lag(e[1], V) * lag(e[2], W)
            total += mp.log(sig)**k * at[bp] * G / DEN[bp](U, V, W)
    return total


if __name__ == '__main__':
    random.seed(3)
    A = {}
    for k in (0, 1):
        for b in ATOMS:
            A[(k, b)] = {(i, j, q): mp.mpf(random.randint(-5, 5)) / 7 for i in range(2) for j in range(2) for q in range(2) if random.random() < 0.4}
    for pt in [('0.62', '0.41', '0.53'), ('0.35', '0.8', '0.6'), ('0.9', '0.5', '1.1')]:
        r1, r2, r12 = map(mp.mpf, pt)
        Uq, Vq, Wq = EPS * (r2 + r12 - r1), EPS * (r1 + r12 - r2), 2 * EPS * (r1 + r2 - r12)
        sig = Uq + Vq + Wq
        f = lambda a, b, c: psi(a, b, c, A)
        lhs = sig**2 * (-(2 * Vq + Wq) * (2 * Uq + Wq) * (Uq + Vq)) * mp.exp(sig / 2) * H_minus_E(f, r1, r2, r12)
        r = rhs(r1, r2, r12, A)
        print(pt, 'LHS', mp.nstr(lhs, 20), ' RHS', mp.nstr(r, 20), ' rel.diff', mp.nstr(abs(lhs - r) / abs(lhs), 3), flush=True)
